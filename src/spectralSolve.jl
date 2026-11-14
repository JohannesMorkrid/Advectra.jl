# ------------------------------------------------------------------------------------------
#                                      Spectral Solve                                       
# ------------------------------------------------------------------------------------------

# Assuming for now that dt is fixed.
# The first step is stored during initialization of Output.
"""
    spectral_solve(prob, scheme, output; debug)

  Integrates `prob` using `scheme`, writing output via `output` at each step.

  ## Arguments
  - `prob`: The `SpectralODEProblem` to solve.
  - `scheme`: The ODE algorithm to use (default: `MSS3()`).
  - `output`: An `Output` instance managing writers and diagnostics
    (default: `Output(prob)` with a single `HDF5Writer`).

  ## Keyword Arguments
  - `debug`: If `true`, rethrows errors instead of printing them (default: `false`).
"""
function spectral_solve(prob::SpectralODEProblem, scheme::AbstractODEAlgorithm=MSS3(),
    output::Output=Output(prob); debug::Bool=false)

    # Initialize cache and tracking from checkpoint if resume=true and one exists
    cache, step, t = initialize_solve(prob, scheme, output)

    dt = prob.dt

    # Calculate number of steps #TODO ceil or floor?
    total_steps = floor(Int, (last(prob.tspan) - first(prob.tspan)) / dt)

    # Enable CTRL+C from terminal outside of interactive mode
    Base.exit_on_sigint(false)

    try
        # This method assumes step number does not overflow!
        while step < total_steps
            perform_step!(cache, prob, t)

            step += 1
            t = first(prob.tspan) + step * dt

            handle_output!(output, step, cache.u, prob, t)
        end
    catch error
        # Interupt the error, so that the code does not halt when not in debug mode
        debug ? rethrow(error) : showerror(stdout, error)
    end

    # Store checkpoint so the simulation can be resumed later
    store_checkpoint!(output, cache, step, t)

    # Final flush
    flush!(output, step)

    return output
end

"""
    initialize_solve(prob, scheme, output)

  Initializes the integrator cache and time tracking. If any writer supports checkpointing
  and holds a stored checkpoint, the cache and time are restored from it.
"""
function initialize_solve(prob::SpectralODEProblem, scheme::AbstractODEAlgorithm,
    output::Output)

    restored = restore_checkpoint(output, prob, scheme)
    isnothing(restored) && return get_cache(prob, scheme), 0, first(prob.tspan)

    # Resume and resync schedules to checkpoint t
    cache, step, t, Δt_prev = restored
    setup_schedules!(output, prob, t; Δt_prev)
    return restored
end