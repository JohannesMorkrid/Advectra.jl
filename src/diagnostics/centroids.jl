# ------------------------------------------------------------------------------------------
#                                      COM Diagnostic                                       
# ------------------------------------------------------------------------------------------

# ---------------------------------------- Helper ------------------------------------------

"""
    radial_centroid(field, prob, time)

  Compute radial centroid position for field (u) `X_C = ∑ux/∑u`, excluding the boundary.
"""
function radial_centroid(field, prob, time)
    # 2:end is because the boundaries are periodic and thus should not contribute
    f = @view field[2:end, 2:end]
    x = prob.domain.x[2:end]'
    return mapreduce((fi, xi) -> fi * xi, +, f, x) / sum(f)
end

"""
    poloidal_centroid(field, prob, time)

  Compute poloidal centroid position for field (u) `Y_C = ∑uy/∑u`, excluding the boundary.
"""
function poloidal_centroid(field, prob, time)
    # 2:end is because the boundaries are periodic and thus should not contribute
    f = @view field[2:end, 2:end]
    y = prob.domain.y[2:end]
    return mapreduce((fi, yi) -> fi * yi, +, f, y) / sum(f)
end

"""
    _centroid_velocity!(memory, R_C, time)
  
  Compute centroid velocity based on current `R_C` and previous `memory` positions.
  
  !!! warning
  The `memory` `Dict` is altered to store the current state as the previous for next sample.
"""
function _centroid_velocity!(memory, R_C, time)
    # Check that do not divide by zero
    if memory["previous_time"] == time
        V_C = 0.0
    else
        V_C = (R_C .- memory["previous_position"]) ./ (time .- memory["previous_time"])
    end

    # Store for next computation
    memory["previous_position"] = R_C
    memory["previous_time"] = time

    return V_C
end

# -------------------------------------- Main Method ---------------------------------------

"""
    centroid(state, prob, time, memory::Dict=Dict(); field=first, dims::Symbol=:radial)

Compute the centroid position `R_C` and velocity `V_C` for a given field and direction.

The field used for the centroid computation is specified via `field`, a function
applied to all slices of `state`. Defaults to `first` (i.e. the first slice).

### Keyword Arguments
- `field`: method applied to splatted slices get desired field, e.g. 
`(ρ, ω, T) -> ρ .* T` (default: `first`)
- `dims`: centroid component, `:radial` or `:poloidal` (default: `:radial`)

### Returns
`Array` where first entry is the position `R_C` and the second is the velocity `V_C`.

### Examples
```julia
  # Radial centroid of density/first field (default)
  centroid(state, prob, time, memory)

  # Poloidal centroid of temperature
  centroid(state, prob, time; field=(n, Ω, T) -> T, dims=:poloidal)

  # Radial centroid of pressure, with velocity tracking
  centroid(state, prob, time, memory; field=(n, Ω, T) -> n .* T)
```
"""
function centroid(state, prob, time, memory::Dict=Dict(); field = (n, args...) -> n,
                  dims::Symbol=:radial)
    slices = eachslice(state; dims=ndims(state))
    f = field(slices...)

    centroid_position = if dims === :radial
        radial_centroid
    elseif dims === :poloidal
        poloidal_centroid
    else
        throw(ArgumentError("Unknown dims: $dims. Expected :radial or :poloidal."))
    end

    R_C = centroid_position(f, prob, time)
    V_C = !isempty(memory) ? _centroid_velocity!(memory, R_C, time) : NaN

    return [R_C, V_C]
end

# ------------------------------------ Radial Centroid -------------------------------------

function build_diagnostic(::Val{:radial_COM}; kwargs...)
    args = (Dict("previous_position" => 0.0, "previous_time" => 0.0),)
    Diagnostic(; name="Radial COM",
               method=centroid,
               metadata="Radial Center-of-mass (COM) diagnostics, columns: X_COM, V_COM.",
               args=args,
               kwargs=(; field=(n, args...) -> n, dims=:radial))
end

function build_diagnostic(::Val{:radial_COT}; kwargs...)
    args = (Dict("previous_position" => 0.0, "previous_time" => 0.0),)
    Diagnostic(; name="Radial COT",
               method=centroid,
               metadata="Radial Center-of-temperature (COT) diagnostics, columns: X_COT, V_COT.",
               args=args,
               kwargs=(; field=(n, Ω, T, args...) -> T, dims=:radial))
end

function build_diagnostic(::Val{:radial_COP}; kwargs...)
    args = (Dict("previous_position" => 0.0, "previous_time" => 0.0),)
    Diagnostic(; name="Radial COP",
               method=centroid,
               metadata="Radial Center-of-pressure (COP) diagnostics, columns: X_COP, V_COP.",
               args=args,
               kwargs=(; field=(n, Ω, T, args...) -> n .* T, dims=:radial))
end

# ----------------------------------- Poloidal Centroid ------------------------------------

function build_diagnostic(::Val{:poloidal_COM}; kwargs...)
    args = (Dict("previous_position" => 0.0, "previous_time" => 0.0),)
    Diagnostic(; name="Poloidal COM",
               method=centroid,
               metadata="Poloidal Center-of-mass (COM) diagnostics, columns: Y_COM, V_COM.",
               args=args,
               kwargs=(; field=(n, args...) -> n, dims=:poloidal))
end

function build_diagnostic(::Val{:poloidal_COT}; kwargs...)
    args = (Dict("previous_position" => 0.0, "previous_time" => 0.0),)
    Diagnostic(; name="Poloidal COT",
               method=centroid,
               metadata="Poloidal Center-of-temperature (COT) diagnostics, columns: Y_COT, V_COT.",
               args=args,
               kwargs=(; field=(n, Ω, T, args...) -> T, dims=:poloidal))
end

function build_diagnostic(::Val{:poloidal_COP}; kwargs...)
    args = (Dict("previous_position" => 0.0, "previous_time" => 0.0),)
    Diagnostic(; name="Poloidal COP",
               method=centroid,
               metadata="Poloidal Center-of-pressure (COP) diagnostics, columns: X_COP, V_COP.",
               args=args,
               kwargs=(; field=(n, Ω, T, args...) -> n .* T, dims=:poloidal))
end
