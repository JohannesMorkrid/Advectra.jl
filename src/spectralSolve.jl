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
  - `resume`: If `true`, restores the integrator state from the most recent checkpoint
    stored in the first capable writer (default: `false`).
  - `debug`: If `true`, rethrows errors instead of printing them (default: `false`).
"""
function spectral_solve(prob::SpectralODEProblem, scheme::AbstractODEAlgorithm=MSS3(),
    output::Output=Output(prob);
    resume::Bool=false,
    debug::Bool=false)

    # Initialize cache and tracking — restores checkpoint if resume=true and one exists
    cache, step, t = initialize_solve(prob, scheme, output; resume)

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
    output::Output; resume::Bool=false)

    restored = resume ? restore_checkpoint(output, prob, scheme) : nothing
    isnothing(restored) ? (get_cache(prob, scheme), 0, first(prob.tspan)) : restored
end