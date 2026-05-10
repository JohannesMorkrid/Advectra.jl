# ------------------------------------------------------------------------------------------
#                                      COM Diagnostic                                       
# ------------------------------------------------------------------------------------------

# ---------------------------------------- Helper ------------------------------------------

"""
    compute_radial_COM_position(field, prob, time)

  Compute radial COM position for field (n) `X_COM = ∑nx/∑n`, excluding the boundary.
"""
function compute_radial_COM_position(field, prob, time)
    # 2:end is because the boundaries are periodic and thus should not contribute
    x = prob.domain.x[2:end]'
    return @views sum(field[2:end, 2:end] .* x) / sum(field[2:end, 2:end])
end

"""
    _compute_radial_COM_velocity!(memory, X_COM, time)
  
  Compute radial COM velocity based on current `X_COM` and previous `memory` positions.
  
  !!! warning
  The `memory` `Dict` is altered to store the current state as the previous for next sample.
"""
function _compute_radial_COM_velocity!(memory, X_COM, time)
    # Check that do not divide by zero
    if memory["previous_time"] == time
        V_COM = 0.0
    else
        V_COM = (X_COM .- memory["previous_position"]) ./ (time .- memory["previous_time"])
    end

    # Store for next computation
    memory["previous_position"] = X_COM
    memory["previous_time"] = time

    return V_COM
end

# -------------------------------------- Main Method ---------------------------------------

"""
    radial_COM(state, prob, time, memory=Dict(); field_idx::Int = 1)

  Compute the radial position `X_COM` and velocity `V_COM` as long as `memory` is passed.

  ### Returns
  `Array` where first column/entry is the position and the second is the velocity.
"""
function radial_COM(state, prob, time, memory::Dict=Dict(); field_idx::Int=1)
    slices = eachslice(state; dims=ndims(state))
    density = slices[field_idx]
    X_COM = compute_radial_COM_position(density, prob, time)

    if !isempty(memory)
        V_COM = _compute_radial_COM_velocity!(memory, X_COM, time)
    else
        V_COM = NaN
    end

    return [X_COM, V_COM]
end

function build_diagnostic(::Val{:radial_COM}; kwargs...)
    args = (Dict("previous_position" => 0.0, "previous_time" => 0.0),)
    Diagnostic(; name="Radial COM",
               method=radial_COM,
               metadata="Radial Center-of-mass (COM) diagnostics, columns: X_COM, V_COM.",
               args=args)
end