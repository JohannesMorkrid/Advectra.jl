# ------------------------------------------------------------------------------------------
#                                      CFL Diagnostic                                       
# ------------------------------------------------------------------------------------------

# ---------------------------------- Velocity Functions ------------------------------------

"""
    compute_velocity(state, prob, time, velocity::Symbol=:ExB)

  Compute the velocity fields, i.e. (v_x, v_y), based on `velocity` method.

  ### `velocity` options:
  - `:ExB`: uses (-∂ϕ∂y, ∂ϕ∂x)
  - `:burger`: uses the state as the velocity magnitude.

  ### Return
  Tuple (v_x, v_y).
"""
function compute_velocity(state, prob, time; velocity::Symbol=:ExB)
    _compute_velocity(state, prob, time, Val(velocity))
end

# Calculate velocity U_ExB = ̂z×∇Φ   
function _compute_velocity(state_hat, prob, time, velocity::Val{:ExB})
    @unpack domain, operators = prob
    @unpack solve_phi, diff_x, diff_y = operators
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    ϕ_hat = solve_phi(n_hat, Ω_hat)
    return (get_bwd(domain) * -diff_y(ϕ_hat), get_bwd(domain) * diff_x(ϕ_hat))
end

function _compute_velocity(state, prob, time, velocity::Val{:burger})
    (selectdim(state, ndims(prob.domain) + 1, 1),)
end

# ---------------------------------- Component Functions -----------------------------------

"""
    compute_cfl(velocities::Tuple, prob, time; component::Symbol=:x)

  Compute CFL-criterion for each velocity, based on which `component` is requested.

  ### `component` options
  - `:x` or `:radial`: return v_x*dt/dx.
  - `:y` or `:poloidal`: return v_y*dt/dy.
  - `:both`: return both :x and :y components.
  - `:magnitude`: return hypot(v_x, v_y, ...)*dt/hypot(dx, dy, ...).
"""
function compute_cfl(velocities::Tuple, prob, time; component::Symbol=:x)
    compute_cfl(velocities, prob, time, Val(component))
end

function compute_cfl(velocities::Tuple, prob, time, component::Union{Val{:x},Val{:radial}})
    return (first(velocities) .* (prob.dt / prob.domain.dx),)
end

function compute_cfl(velocities::Tuple, prob, time,
                     component::Union{Val{:y},Val{:poloidal}})
    return (last(velocities) .* (prob.dt / prob.domain.dy),)
end

function compute_cfl(velocities::Tuple, prob, time, component::Val{:both})
    return (first(velocities) .* (prob.dt / prob.domain.dx),
            last(velocities) .* (prob.dt / prob.domain.dy))
end

function compute_cfl(velocities::Tuple, prob, time, component::Val{:magnitude})
    Δ = hypot(differential_elements(prob.domain)...)
    return (hypot.(velocities...) .* (prob.dt / Δ),)
end

function compute_cfl(velocities::Tuple, prob, time, component::Val{T}) where {T}
    error("Unknown component: :$T")
end

# ----------------------------------- Printing Methods -------------------------------------

"""
    print_cfl(maximas::Tuple, prob, time, component::Val{Symbol})

  Print the max CFL number alongside index and time in a human readable format.
"""
function print_cfl(maximas::Tuple, prob, time, component::Union{Val{:x},Val{:radial}})
    maxima = first(maximas)
    println("CFL max (x): $(round(first(maxima), sigdigits=3))",
            " @ $(last(maxima).I), t=$(round(time, digits=3))")
end

function print_cfl(maximas::Tuple, prob, time, component::Union{Val{:y},Val{:poloidal}})
    maxima = first(maximas)
    println("CFL max (y): $(round(first(maxima), sigdigits=3))",
            " @ $(last(maxima).I), t=$(round(time, digits=3))")
end

function print_cfl(maximas::Tuple, prob, time, component::Val{:both})
    x_maxima, y_maxima = maximas
    println("CFL max (x): $(round(first(x_maxima), sigdigits=3)) @ $(last(x_maxima).I)",
            ", max (y): $(round(first(y_maxima), sigdigits=3)) @ $(last(y_maxima).I)",
            ", t=$(round(time, digits=3))")
end

function print_cfl(maximas::Tuple, prob, time, component::Val{:magnitude})
    maxima = first(maximas)
    println("CFL max (magnitude): $(round(first(maxima), sigdigits=3))",
            " @ $(last(maxima).I), t=$(round(time, digits=3))")
end

# -------------------------------------- Main Method ---------------------------------------

"""
    cfl(state, prob, time, velocity=Val(:ExB), component=Val(:x); silent=false)

  Compute CFL based on which `velocity` and what `component`, if `silent=false` then the 
  results are printed. 
  
  See [`compute_velocity`](@ref) and [`compute_cfl`](@ref) for available options for 
  `velocity` and `component` respectively.  

  ### Return
  `Array` where row represent the component and the first column is the max cfl, with the
  remaining columns being indices.
""" # TODO consider a better output format (not worthy of an issue)
function cfl(state, prob, time, velocity=Val(:ExB), component=Val(:x); silent::Bool=false)
    velocities = _compute_velocity(state, prob, time, velocity)
    CFLs = compute_cfl(velocities, prob, time, component)

    maximas = findmax.(CFLs)
    results = stack([[first(maxima), last(maxima).I...] for maxima in maximas]; dims=1)

    !silent ? print_cfl(maximas, prob, time, component) : nothing # Catch breakdown?

    return results
end

function cfl(state::AbstractArray, prob, time; velocity::Symbol=:ExB, component::Symbol=:x,
             silent::Bool=false)
    cfl(state, prob, time, Val(velocity), Val(component); silent=silent)
end
# ------------------------------------- Build Related --------------------------------------

# :x, :radial, :y, :poloidal, :both, :magnitude
function cfl_metadata(context::String, component::Symbol)
    components = (component == :both) ? [:x, :y] : [component]
    rows = join(components, "-component; ") * "-component"
    indices = "y-index; x-index"
    return "$context Courant-Friedrichs-Lewy number, rows: $rows, columns: max CFL; $indices."
end

assumes_spectral(::Val{:burger}) = false
# Default
assumes_spectral(::Val{T}) where {T} = true

function requires_operator(::Val{cfl}; velocity_method, kwargs...)
    if velocity_method == :ExB
        return [OperatorRecipe(:diff_x), OperatorRecipe(:diff_y),
                OperatorRecipe(:solve_phi)]
    else
        return []
    end
end

function build_diagnostic(::Val{:cfl}; velocity=:ExB, component=:magnitude,
                          silent=false, kwargs...)
    context = isuppercase(string(velocity)[1]) ? string(velocity) :
              titlecase(string(velocity))
    Diagnostic(; name="$context CFL",
               method=cfl,
               metadata=cfl_metadata(context, component),
               assumes_spectral_state=assumes_spectral(Val(velocity)),
               args=(Val(velocity), Val(component)),
               kwargs=(; silent=silent))
end