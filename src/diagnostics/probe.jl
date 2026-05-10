# ------------------------------------------------------------------------------------------
#                                     Probe Diagnostics                                     
# ------------------------------------------------------------------------------------------

# ------------------------------------ Probing Methods -------------------------------------

# ------------------------------------- One Position ---------------------------------------

"""
    get_index(position, domain)

  Return `Tuple` of indices, the closest index for each dimension.
"""
function get_index(position::NTuple{2,Number}, domain::AbstractDomain)
    return argmin(abs.(domain.y .- position[2])), argmin(abs.(domain.x .- position[1]))
end

"""
    probe_field(field::AbstractArray, domain::AbstractDomain, position::Tuple)
    probe_field(field::AbstractInterpolation, domain::AbstractDomain, position::Tuple)

  Probe the `field` at one `position` using the `domain` to determine the index.
"""
function probe_field(field::AbstractArray, domain::AbstractDomain, position)
    @allowscalar field[get_index(position, domain)...]
end

function probe_field(field::AbstractInterpolation, domain::AbstractDomain,
                     position::NTuple{2,Number})
    field(position[2], position[1])
end

# ---------------------------------- Multiple Positions ------------------------------------

"""
    probe_field(field::AbstractArray, domain::AbstractDomain, positions::AbstractArray,
                     interpolation::Nothing=nothing)
    probe_field(field::AbstractArray, domain::AbstractDomain, positions::AbstractArray,
                     interpolation::AbstractInterpolation)

  Probe the `field` at multiple `positions` using the `domain` to determine the indices. 
  The method is dispatched on the `interpolation` method.
"""
function probe_field(field::AbstractArray, domain::AbstractDomain,
                     positions::AbstractArray{<:Tuple},
                     interpolation::Nothing=nothing)
    data = [probe_field(field, domain, position) for position in positions]

    # Return either the one point, or the array
    length(data) == 1 ? data[1] : data
end

# ------------------------------------- Indices Based --------------------------------------

# Uses pre-computed indices
function probe_field(field::AbstractArray, domain::AbstractDomain,
                     indices::AbstractArray{<:Integer},
                     interpolation::Nothing=nothing)
    data = field[indices]
    # Return either the one point, or the array
    length(data) == 1 ? only(Array(data)) : data
end

# GPU optimized in-place without interpolation
function probe_field!(out::T, field::T, indices::T) where {T<:AbstractGPUArray}
    i = threadIdx().x
    if i <= length(indices)
        out[i] = field[indices[i]]
    end
    return
end

# TODO fix SubArray issue [#26](https://github.com/JohannesMorkrid/Advectra.jl/issues/26)
# GPU "optimized" out-of place without interpolation
# function probe_field(field::AbstractGPUArray, domain::AbstractDomain,
#                      indices::AbstractGPUArray{<:Integer},
#                      interpolation::Nothing=nothing)
#     data = zeros(size(indices)) |> memory_type(domain)
#     probe_field!(data, field, indices)

#     # Return either the one point, or the array
#     length(data) == 1 ? data[1] : data
# end

# --------------------------------- Interpolated Probing -----------------------------------

function probe_field(field::AbstractArray, domain::AbstractDomain, positions::AbstractArray,
                     interpolation::AbstractInterpolation)
    interpolated = interpolation((domain.y, domain.x), field)
    data = [probe_field(interpolated, domain, position) for position in positions]
    # Return either the one point, or the array
    length(data) == 1 ? data[1] : data
end

# --------------------------------- Construction Related -----------------------------------

"""
    ensure_in_domain(position, domain::Domain)

  Checks if the `position` is inside the domain. Throws error if not.
"""
function ensure_in_domain(position::NTuple{2,Number}, domain::Domain)
    if (first(domain.x) <= first(position)) && (last(domain.x) >= first(position)) &&
       (first(domain.y) <= last(position)) && (last(domain.y) >= last(position))
        nothing
    else
        error("Point at position $position is outside the domain bounds \
        x∈[$(first(domain.x)), $(last(domain.x))], y∈[$(first(domain.y)), $(last(domain.y))].")
    end
end

"""
    validate_positions(positions, domain)

  Check that all points have the correct length and are within the bounds of the domain.
"""
function validate_positions(positions, domain)
    points = get_points(domain)
    ndim = length(points)
    for (i, position) in enumerate(positions)
        if length(position) != ndim
            error("Point $i has an invalid length. Expected a tuple of length " *
                  string(ndim) * ", but got length $(length(position)).")
        end

        # May throw an error
        ensure_in_domain(position, domain)
    end
end

"""
    prepare_positions(positions, domain)

  Convert positions into an `Array{<:Tuples}`, and check if the positions are in the domain. 
"""
function prepare_positions(positions, domain)
    # Check if the user sent in tuple of points or single point
    if isa(positions, Tuple) && isa(positions[1], Number)
        positions = [positions]
    end

    # Convert all positions to the desired precision type
    PrecisionType = get_precision(domain)
    positions = [PrecisionType.(position) for position in positions]

    # Bound and length checks the positions
    validate_positions(positions, domain)

    return positions
end

"""
    prepare_indices(positions, domain::Domain)

  Prepare an `Array` of `LinearIndices` related to the fields on the domain, for quicker
  parallel access.
"""
function prepare_indices(positions, domain::Domain)
    Ind = LinearIndices((size(domain)))
    [Ind[get_index(position, domain)...] for position in positions] |> domain.MemoryType
end

"""
    prepare_metadata(; positions::AbstractArray{<:Tuple}, quantities::String)
    prepare_metadata(; positions::AbstractArray{<:Tuple}, quantities::AbstractArray{<:String})

  Prepare a human readable overview of positions of probes and what is measured.
"""
function prepare_metadata(positions::AbstractArray{<:Tuple}, quantities::AbstractArray)
    return string("Probe at positions: ", join(positions, "; "), ", measuring: ",
                  join(quantities, ", "), ".")
end

function prepare_metadata(positions::AbstractArray{<:Tuple}, quantities::String)
    return string("Probe at positions: ", join(positions, "; "), ", measuring: ",
                  quantities, ".")
end

"""
    build_probe_diagnostic(; name, method, positions, domain, quantities, 
    assumes_spectral_state=false, interpolation)
  
  General build method for probe diagnostics which prepares the positions and metadata.
"""
function build_probe_diagnostic(; name, method, positions, domain, quantities,
                                assumes_spectral_state=false, interpolation=nothing)
    # Bound and type checking
    positions = prepare_positions(positions, domain)
    # Can pre-compute indices if not using interpolation
    if isnothing(interpolation)
        args = (prepare_indices(positions, domain),)
    else
        args = (positions,)
    end
    Diagnostic(; name=name,
               method=method,
               metadata=prepare_metadata(positions, quantities),
               assumes_spectral_state=assumes_spectral_state,
               args=args,
               kwargs=(; interpolation=interpolation,))
end

# ------------------------------------------------------------------------------------------
#                                      Specific Probes                                      
# ------------------------------------------------------------------------------------------

# ---------------------------------------- Density -----------------------------------------

"""
    probe_density(state, prob, time; positions, interpolation=nothing)
  
  Probe the density field, n, at the given `positions`. 
"""
function probe_density(state, prob, time, positions; interpolation=nothing)
    n = selectdim(state, ndims(prob.domain) + 1, 1)
    probe_field(n, prob.domain, positions, interpolation)
end

function build_diagnostic(::Val{:probe_density}; domain, positions,
                          interpolation=nothing, kwargs...)
    build_probe_diagnostic(; name="Density probe",
                           method=probe_density,
                           positions=positions,
                           domain=domain,
                           quantities="density",
                           interpolation=interpolation)
end

# --------------------------------------- Vorticity ----------------------------------------

"""
    probe_vorticity(state, prob, time; positions, interpolation=nothing)
  
  Probe the vorticity field, Ω, at the given `positions`. 
"""
function probe_vorticity(state, prob, time, positions; interpolation=nothing)
    Ω = selectdim(state, ndims(prob.domain) + 1, 2)
    probe_field(Ω, prob.domain, positions, interpolation)
end

function build_diagnostic(::Val{:probe_vorticity}; domain, positions,
                          interpolation=nothing, kwargs...)
    build_probe_diagnostic(; name="Vorticity probe",
                           method=probe_vorticity,
                           positions=positions,
                           domain=domain,
                           quantities="vorticity",
                           interpolation=interpolation)
end

# --------------------------------------- Temperature ----------------------------------------

"""
    probe_temperature(state, prob, time; positions, interpolation=nothing)
  
  Probe the temperature field, Ω, at the given `positions`. 
"""
function probe_temperature(state, prob, time, positions; interpolation=nothing)
    Ω = selectdim(state, ndims(prob.domain) + 1, 3)
    probe_field(Ω, prob.domain, positions, interpolation)
end

function build_diagnostic(::Val{:probe_temperature}; domain, positions,
                          interpolation=nothing, kwargs...)
    build_probe_diagnostic(; name="Temperature probe",
                           method=probe_temperature,
                           positions=positions,
                           domain=domain,
                           quantities="temperature",
                           interpolation=interpolation)
end

# --------------------------------------- Potential ----------------------------------------

"""
    probe_potential(state, prob, time; positions, interpolation=nothing)
  
  Probe the potential field, ϕ, at the given `positions`. 
"""
function probe_potential(state_hat, prob, time, positions; interpolation=nothing)
    @unpack domain, operators = prob
    @unpack solve_phi = operators
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    ϕ = get_bwd(domain) * solve_phi(n_hat, Ω_hat)
    probe_field(ϕ, domain, positions, interpolation)
end

requires_operator(::Val{:probe_potential}; kwargs...) = [OperatorRecipe(:solve_phi)]

function build_diagnostic(::Val{:probe_potential}; domain, positions,
                          interpolation=nothing, kwargs...)
    build_probe_diagnostic(; name="Phi probe",
                           method=probe_potential,
                           positions=positions,
                           domain=domain,
                           quantities="potential",
                           assumes_spectral_state=true,
                           interpolation=interpolation)
end

# ------------------------------------ Radial Velocity -------------------------------------

"""
    probe_radial_velocity(state, prob, time; positions, interpolation=nothing)
  
  Probe the radial velocity field, vᵣ, at the given `positions`. 
"""
function probe_radial_velocity(state, prob, time, positions; interpolation=nothing)
    @unpack domain, operators = prob
    @unpack solve_phi, diff_y = operators
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    ϕ_hat = solve_phi(n_hat, Ω_hat)
    v_x_hat = -diff_y(ϕ_hat)
    v_x = get_bwd(domain) * v_x_hat
    probe_field(v_x, domain, positions, interpolation)
end

function requires_operator(::Val{:probe_radial_velocity}; kwargs...)
    [OperatorRecipe(:solve_phi), OperatorRecipe(:diff_y)]
end

function build_diagnostic(::Val{:probe_radial_velocity}; domain, positions,
                          interpolation=nothing, kwargs...)
    build_probe_diagnostic(; name="Radial velocity probe",
                           method=probe_radial_velocity,
                           positions=positions,
                           domain=domain,
                           quantities="radial velocity",
                           assumes_spectral_state=true,
                           interpolation=interpolation)
end

# -------------------------------------- All Fields ----------------------------------------

"""
    probe_all(state, prob, time; positions, interpolation=nothing)
  
  Probe the density (n), vorticity (Ω), potential (ϕ), radial velocity field (vᵣ) and 
  radial flux (Γ), at the given `positions` all "at once". 
"""
function probe_all(state_hat, prob, time, positions; interpolation=nothing)
    @unpack domain, operators = prob
    @unpack solve_phi, diff_y = operators

    # Calculate spectral fields
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    ϕ_hat = solve_phi(n_hat, Ω_hat)
    v_x_hat = -diff_y(ϕ_hat)

    # Cache for transformation
    cache = zeros(size(domain)) |> memory_type(prob.domain)

    # Transform to physical space and probe fields
    n = mul!(cache, get_bwd(domain), n_hat)
    n_p = probe_field(n, domain, positions, interpolation)
    Ω = mul!(cache, get_bwd(domain), Ω_hat)
    Ω_p = probe_field(Ω, domain, positions, interpolation)
    ϕ = mul!(cache, get_bwd(domain), ϕ_hat)
    ϕ_p = probe_field(ϕ, domain, positions, interpolation)
    v_x = mul!(cache, get_bwd(domain), v_x_hat)
    v_x_p = probe_field(v_x, domain, positions, interpolation)

    #Combine fields for output (The last field is the flux Γ=nvₓ)
    cat(n_p, Ω_p, ϕ_p, v_x_p, n_p .* v_x_p; dims=3)
end

function requires_operator(::Val{:probe_all}; kwargs...)
    [OperatorRecipe(:solve_phi), OperatorRecipe(:diff_y)]
end

function build_diagnostic(::Val{:probe_all}; domain, positions, interpolation=nothing,
                          kwargs...)
    build_probe_diagnostic(; name="All probe",
                           method=probe_all,
                           positions=positions,
                           domain=domain,
                           quantities=["density",
                               "vorticity",
                               "potential",
                               "radial velocity",
                               "radial flux"],
                           assumes_spectral_state=true,
                           interpolation=interpolation)
end