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
  - `strides`: Per-diagnostic stride (steps between samples).
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
  - `resume`: Resume from an existing output file (default: `false`).

  # Usage example

    # Single writer (default)
    output = Output(prob)

    # Custom writers
    output = Output(prob;
        writers = [
            HDF5Writer(filename="run_full.h5"),
            HDF5Writer(filename="run_scalars.h5"),
            LocalWriter(),
        ],
        physical_transform = my_transform,
        simulation_name    = :timestamp,
        storage_limit      = "10 GB"
    )
"""
mutable struct Output{DV<:AbstractArray{<:Diagnostic},UB<:AbstractArray,
    WV<:AbstractVector{<:AbstractWriter},PT<:Function}
    diagnostics::DV
    strides::Vector{Int} # TODO remove this in favor of the schedules stored in diagnostics
    state_buffer::UB
    writers::WV
    physical_transform::PT
    transformed::Bool
end

function Output(prob::SpectralODEProblem;
    writers::AbstractVector{<:AbstractWriter}=[HDF5Writer()],
    physical_transform::PT=identity,
    storage_limit::AbstractString="",
    simulation_name::Union{AbstractString,Symbol,Nothing}="advectra",
    resume::Bool=false) where {PT<:Function}

    # Prepare physical-space buffer and sample each diagnostic once at t0
    state, t0 = prepare_initial_state(prob; physical_transform=physical_transform)
    diagnostics, initial_samples = initialize_diagnostics(prob, state, prob.u0_hat, t0)
    strides = determine_strides(initial_samples, prob, storage_limit)

    # Resolve simulation name once — all writers share the same resolved string
    simulation_name = resolve_simulation_name(simulation_name, prob)

    # Give every writer a chance to build its file structure
    for writer in writers
        initialize!(writer, prob, diagnostics, initial_samples, t0; simulation_name, resume)
    end

    return Output(diagnostics, strides, state, writers, physical_transform, true)
end

# ------------------------------- Convenience Constructors ---------------------------------

"""
    Output(prob, writer; kwargs...)

  Convenience form accepting a single writer without wrapping it in a vector.
"""
Output(prob::SpectralODEProblem, writer::AbstractWriter; kwargs...) =
    Output(prob; writers=[writer], kwargs...)

# --------------------------------- Simulation Name Resolution ----------------------------

"""
    resolve_simulation_name(simulation_name, prob) -> String

  Resolves `simulation_name` into a concrete HDF5 group name string. Called once by
  `Output` before passing the result to each writer's `initialize!`.

  ### Supported values
  * `nothing`      — `"/"` (flat structure, write directly to file root)
  * `:timestamp`   — `Dates.now()` formatted as a string (unique per run)
  * `:parameters`  — parameter names and values joined into a string
  * A plain `String` — used as-is
"""
function resolve_simulation_name(simulation_name, prob)
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

# --------------------------------- Checkpoint Interface -----------------------------------

"""
    store_checkpoint!(output, cache, step, time)

  Fans out to all writers that support checkpointing. All supporting writers receive the
  checkpoint — redundancy is cheap since checkpointing is rare and data is small.
"""
function store_checkpoint!(output::Output, cache::AbstractCache, step::Integer, time::Number)
    for writer in output.writers
        store_checkpoint!(writer, cache, step, time)
    end
end

"""
    restore_checkpoint(output, prob, scheme) -> Union{Tuple, Nothing}

  Finds the first writer with a stored checkpoint and restores `cache, step, t` from it.
  Returns `nothing` if no checkpoint is found, leaving the caller to handle the fresh-start
  fallback.
"""
function restore_checkpoint(output::Output, prob::SpectralODEProblem,
    scheme::AbstractODEAlgorithm)
    checkpoint_writer = findfirst(has_checkpoint, output.writers)
    isnothing(checkpoint_writer) && return nothing
    cache, step, t = restore_checkpoint(output.writers[checkpoint_writer], prob, scheme)
    @info "Resuming from checkpoint at t=$t (step=$step)"
    return cache, step, t
end

# ------------------------------------ Base overloads --------------------------------------

flush!(output::Output, step) = foreach(w -> flush!(w, step), output.writers)

function Base.close(output::Output)
    for writer in output.writers
        close(writer)
    end
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