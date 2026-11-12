# ------------------------------------------------------------------------------------------
#                                   Diagnostics Sampling                                    
# ------------------------------------------------------------------------------------------

# ------------------------------- Initialization ------------------------------------------

"""
    initialize_diagnostics(prob, state, state_hat, t0)

  Builds each `Diagnostic` from `prob.diagnostic_recipes` and samples it once at `t0`.
  Returns the vector of diagnostics and the vector of initial samples.
"""
function initialize_diagnostics(prob::SpectralODEProblem, state, state_hat, t0)
    diagnostics = Diagnostic[]
    initial_samples = []

    for recipe in prob.diagnostic_recipes
        diagnostic = build_diagnostic(recipe.method, prob; recipe.kwargs...)
        push!(diagnostics, diagnostic)

        # Determine if sampling is done in physical or spectral space
        input = diagnostic.assumes_spectral_state ? state_hat : state
        # Sample once and store output
        sample = diagnostic(input, prob, t0)
        push!(initial_samples, sample)
    end

    return diagnostics, initial_samples
end

# ----------------------------------- Runtime Handling -------------------------------------

"""
    handle_output!(output, step, state, prob, time)

  Called by the integrator after every time step. Removes aliased modes, samples all due
  diagnostics (transforming the state at most once per step regardless of how many writers
  are active), flushes all writers, and checks for NaN breakdown.
"""
function handle_output!(output::Output, step::Integer, state::AbstractArray,
    prob::SpectralODEProblem, time::Number)

    output.transformed = false

    # Remove modes after each step using user defined function
    prob.remove_modes(state, prob.domain)

    handle_diagnostics!(output, step, state, prob, time)

    flush!(output, step)

    # Check if first value is NaN, if one value is NaN the whole Array will turn NaN after FFT
    assert_no_nan(state, time)
end

"""
    handle_diagnostics!(output, step, state, prob, time)

  Iterates over all diagnostics and samples those that are due at this `step`.
"""
function handle_diagnostics!(output::Output, step::Integer, state, prob, time)
    for diagnostic in output.diagnostics
        if should_sample(diagnostic, step, time)
            sample_diagnostic!(output, diagnostic, step, state, prob, time)
            advance_schedule!(diagnostic)
        end
    end
end

"""
    sample_diagnostic!(output, diagnostic, step, state, prob, time)

  Evaluates `diagnostic` and dispatches the result to every writer. Transforms the state to
  physical space first if the diagnostic requires it and this hasn't been done yet this step.
  Each writer determines its own write position, no `idx` is passed.
"""
function sample_diagnostic!(output::Output, diagnostic, step::Integer, state, prob, time)
    # Check if diagnostic assumes physical field and transform if not yet done
    if !diagnostic.assumes_spectral_state && !output.transformed
        # Transform state (updates state_buffer)
        transform_state!(output, state, bwd_plan(prob.domain))
    end

    # Apply diagnostic to the correct input (either physical or spectral state)
    input = diagnostic.assumes_spectral_state ? state : output.state_buffer
    sample = diagnostic(input, prob, time)

    isnothing(sample) && return

    for writer in output.writers
        write_to!(writer, diagnostic, sample, step, time)
    end
end

"""
    transform_state!(output, state_hat, plan)

  Applies the inverse spectral transform and `output.physical_transform` to
  `output.state_buffer`. Sets `output.transformed = true` so this is done at most once
  per step even when multiple diagnostics need the physical state.
"""
function transform_state!(output::Output, state_hat, plan)
    spectral_transform!(output.state_buffer, plan, state_hat)
    output.physical_transform(output.state_buffer)
    output.transformed = true
end
