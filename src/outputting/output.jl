# ------------------------------------------------------------------------------------------
#                                          Output                                           
# ------------------------------------------------------------------------------------------

"""
    Output

  Coordinates simulation output: owns the diagnostics, manages the physical-space transform
  buffer, drives all writers, and handles checkpointing. The simulation loop calls only
  `handle_output!` and `save_checkpoint!` — it never touches writers directly.

  # Fields
  - `diagnostics`: Array of `Diagnostic`s to evaluate each step.
  - `state_buffer`: Pre-allocated buffer for the physical-space state.
  - `writers`: Vector of `AbstractWriter`s that receive each sample.
  - `physical_transform`: Function applied to the state after the inverse transform.
  - `transformed`: Flag — `true` if the state has already been transformed this step.

  # Constructor

    Output(prob; writers, physical_transform, storage_limit, simulation_name, resume)

  ## Keyword Arguments
  - `writers`: One or more `AbstractWriter` instances (default: `[HDF5Writer()]`).
  - `physical_transform`: Applied after the inverse spectral transform (default: `identity`).
  - `storage_limit`: Optional total storage cap as a string, e.g. `"2 GB"` (default: `""`).
  - `simulation_name`: Name for the top-level group in each writer (default: `"advectra"`).
    Supported symbols: `:timestamp`, `:parameters`, `nothing` (flat, no group), or a
    plain `String`.
  - `resume`: If `true`, appends to existing output files and restores the checkpoint.
    If `false` and a file exists, an error is raised — pass `overwrite_existing=true`
    on the writer to explicitly overwrite (default: `false`).
  - `print_report`: If `true`, prints the storage estimate report at construction time
    (default: `true`).

  # Usage examples

    # First run — file does not exist yet
    output = Output(prob; writers=[HDF5Writer(filename="run.h5")])

    # Resume after walltime — appends to existing file and restores checkpoint
    output = Output(prob; writers=[HDF5Writer(filename="run.h5")], resume=true)

    # Start fresh explicitly — overwrites existing file (always warns)
    output = Output(prob; writers=[HDF5Writer(filename="run.h5", overwrite_existing=true)])

    # Multiple writers with custom diagnostics
    output = Output(prob;
        writers            = [HDF5Writer(filename="run.h5"), LocalWriter()],
        physical_transform = my_transform,
        simulation_name    = :timestamp,
        storage_limit      = "10 GB"
    )
"""
mutable struct Output{DV<:AbstractArray{<:Diagnostic},UB<:AbstractArray,
                      WV<:AbstractVector{<:AbstractWriter},PT<:Function}
    diagnostics::DV # TODO perhaps diagnostics will be removed all together from Output? Or removed from SpectralODEProblem?
    state_buffer::UB
    writers::WV
    physical_transform::PT
    transformed::Bool
    resume::Bool
end

function Output(prob::SpectralODEProblem;
                writers::AbstractVector{<:AbstractWriter}=[HDF5Writer()],
                physical_transform::PT=identity,
                storage_limit::Union{AbstractString,Nothing}="",
                simulation_name::Union{AbstractString,Symbol,Nothing}="advectra",
                resume::Bool=false,
                print_report::Bool=true) where {PT<:Function}

    # Prepare physical-space buffer and build diagnostics. Diagnostics whose schedule is
    # driven by storage_limit are built with Every(-1, typemax(Int)) and resolved below.
    state, t0 = prepare_initial_state(prob; physical_transform=physical_transform)
    diagnostics, initial_samples = initialize_diagnostics(prob, state, prob.u0_hat, t0)

    # Resolve storage_limit-driven schedules (infer Every(n) from per-diagnostic storage_limit)
    resolve_schedules!(diagnostics, initial_samples, prob)

    # Estimate storage for all diagnostics (pure, no side effects)
    estimates = storage_estimates(diagnostics, initial_samples, prob)

    # Optionally print the report
    print_report && print_storage_report(estimates, writers)

    # Validate per-diagnostic and total limits (errors if exceeded)
    validate_storage!(estimates, prob, writers, storage_limit)

    # Resolve simulation name once — all writers share the same resolved string
    simulation_name = resolve_simulation_name(simulation_name, prob)

    # Give every writer a chance to build its file structure
    for writer in writers
        initialize!(writer, prob, diagnostics, initial_samples, t0; simulation_name, resume)
    end

    return Output(diagnostics, state, writers, physical_transform, true, resume)
end

# ------------------------------- Convenience Constructors ---------------------------------

"""
    Output(prob, writer; kwargs...)

  Convenience form accepting a single writer without wrapping it in a vector.
"""
function Output(prob::SpectralODEProblem, writer::AbstractWriter; kwargs...)
    Output(prob; writers=[writer], kwargs...)
end

# --------------------------------- Simulation Name Resolution ----------------------------

"""
    resolve_simulation_name(simulation_name, prob)

  Resolves `simulation_name` into a concrete HDF5 group name string. Called once by
  `Output` before passing the result to each writer's `initialize!`.

  ### Supported values
  * `nothing`      — `"/"` (flat structure, write directly to file root)
  * `:timestamp`   — `Dates.now()` formatted as a string (unique per run)
  * `:parameters`  — parameter names and values joined into a string
  * A plain `String` — used as-is
"""
function resolve_simulation_name(simulation_name, prob)::String
    simulation_name === nothing && return "/"
    simulation_name == :timestamp && return "$(now())"
    simulation_name == :parameters && return parameter_string(prob.p)
    simulation_name isa String && return simulation_name
    error("$simulation_name is not a valid simulation_name. " *
          "Use nothing, :timestamp, :parameters, or a plain String.")
end

# --------------------------------- Prepare Initial State ---------------------------------

"""
    prepare_initial_state(prob; physical_transform=identity)

  Copies `prob.u0`, applies `physical_transform` in-place, and returns it alongside `t0`.
"""
function prepare_initial_state(prob; physical_transform=identity)
    state = copy(prob.u0)
    physical_transform(state)
    t0 = first(prob.tspan)
    return state, t0
end

# ------------------------------- Diagnostic Initialization --------------------------------

"""
    initialize_diagnostics(prob, state, state_hat, t0)
 
  Builds each `Diagnostic` from `prob.diagnostic_recipes` and takes an initial sample at
  `t0`. Schedule alignment is left to `resolve_schedules!`, which runs immediately after.
"""
function initialize_diagnostics(prob, state, state_hat, t0)
    diagnostics = Diagnostic[]
    initial_samples = []

    for recipe in prob.diagnostic_recipes
        diagnostic = build_diagnostic(recipe.method, prob; schedule=recipe.schedule,
                                      recipe.kwargs...)

        push!(diagnostics, diagnostic)

        input = diagnostic.assumes_spectral_state ? state_hat : state
        sample = diagnostic(input, prob, t0)
        push!(initial_samples, sample)
    end

    return diagnostics, initial_samples
end

"""
    setup_schedules!(output, prob, t; Δt_prev=nothing)

  Realigns all diagnostic schedules to the current time `t`. Called at `Output`
  construction (with `t = t0`) and again in `initialize_solve` after restoring a
  checkpoint (with `t = checkpoint_t`). `Δt_prev` is the timestep from the previous
  run — when provided and different from `prob.dt`, `Every` will warn that the sampling
  interval in physical time has changed.
"""
function setup_schedules!(output::Output, prob, t; Δt_prev=nothing)
    for diagnostic in output.diagnostics
        setup!(diagnostic, prob.dt, first(prob.tspan), t; Δt_prev)
    end
end

# --------------------------------- Checkpoint Interface ----------------------------------

"""
    store_checkpoint!(output, cache, step, time)

  Fans out to all writers that support checkpointing. All supporting writers receive the
  checkpoint — redundancy is cheap since checkpointing is rare and data is small.
"""
function store_checkpoint!(output::Output, cache::AbstractCache, step::Integer,
                           time::Number)
    for writer in output.writers
        store_checkpoint!(writer, cache, step, time)
    end
end

"""
    restore_checkpoint(output, prob, scheme)

  Finds the first writer with a stored checkpoint and restores `(cache, step, t, Δt_prev)`
  from it. Returns `nothing` if no checkpoint is found, leaving the caller to handle the
  fresh-start fallback. `Δt_prev` is the timestep of the previous run, used by
  `setup_schedules!` to warn if it changed.
"""
function restore_checkpoint(output::Output, prob::SpectralODEProblem,
                            scheme::AbstractODEAlgorithm)
    output.resume || return nothing
    checkpoint_writer = findfirst(has_checkpoint, output.writers)
    isnothing(checkpoint_writer) && return nothing
    restored = restore_checkpoint(output.writers[checkpoint_writer], prob, scheme)
    cache, step, t, Δt_prev = restored
    @info "Resuming from checkpoint at t=$t (step=$step)"
    return restored
end

# ------------------------------------ Base overloads --------------------------------------

flush!(output::Output, step) = foreach(w -> flush!(w, step), output.writers)

Base.close(output::Output) =
    for writer in output.writers
        close(writer)
    end

function Base.show(io::IO, m::MIME"text/plain", output::Output)
    print(io, "Output")
    output.physical_transform !== identity &&
        print(io, " (physical_transform: ", nameof(output.physical_transform), ")")
    print(io, "\nWriters:")
    for writer in output.writers
        print(io, "\n  ")
        show(io, m, writer)
    end
    print(io, "\nDiagnostics:")
    for diagnostic in output.diagnostics
        print(io, "\n  - ")
        show(io, m, diagnostic)
    end
end

# ---------------------------------- Handling Of Output ------------------------------------

"""
    # Store final state
    # function write_final_state()
    # Get finale index
    #idx = "end"

    #     # Store in hdf if user wants
    #     output.store_hdf ? write_state(output.simulation, idx, U, t) : nothing

    #     # Store in locally if user wants
    #     output.store_locally ? write_local_state(output, idx, U, t) : nothing
    # end
"""
