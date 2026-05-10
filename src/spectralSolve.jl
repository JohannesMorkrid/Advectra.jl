# Assuming for now that dt is fixed
# If custom outputter is not provided, then resort to default
# First step is stored during initilization of output
function spectral_solve(prob::SOP, scheme::SA=MSS3(),
                        output::O=Output(prob; store_hdf=false);
                        debug::Bool=false,) where {SOP<:SpectralODEProblem,
                                                   SA<:AbstractODEAlgorithm,O<:Output}

    # Initialize cache and tracking
    cache, t, step = initialize_solve(prob, scheme, output)

    # Time step
    dt = prob.dt

    # Calculate number of steps #TODO ceil or floor?
    total_steps = floor(Int, (last(prob.tspan) - first(prob.tspan)) / dt)

    # Enable CTRL+C from terminal outside of interactive mode
    Base.exit_on_sigint(false)

    try
        # This method assumes step number does not overflow!
        while step < total_steps
            perform_step!(cache, prob, t)

            # Increment step and time
            step += 1
            t = first(prob.tspan) + step * dt

            handle_output!(output, step, cache.u, prob, t)
        end
    catch error
        # Interupt the error, so that the code does not halt when not in debug mode
        debug ? rethrow(error) : showerror(stdout, error)
    end

    # Store the cache to be able to resume simulations
    save_checkpoint!(output, cache, step, t)

    # TODO catch edge case

    # Write buffer to file
    output.store_hdf ? flush(output.simulation.file) : nothing

    # Returns output struct
    return output
end

function initialize_solve(prob::SOP, scheme::SA,
                          output::O) where {SOP<:SpectralODEProblem,
                                            SA<:AbstractODEAlgorithm,O<:Output}
    if output.resume && output.store_hdf && haskey(output.simulation, "checkpoint")
        cache = restore_checkpoint(output.simulation, prob, scheme)
        t = read(output.simulation, "checkpoint/time")
        step = read(output.simulation, "checkpoint/step")
    else
        cache = get_cache(prob, scheme)
        t = first(prob.tspan)
        step = 0
    end
    return cache, t, step
end
