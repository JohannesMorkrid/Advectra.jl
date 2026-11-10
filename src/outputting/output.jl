# ------------------------------------------------------------------------------------------
#                                         Output Struct
# ------------------------------------------------------------------------------------------

# TODO fix function description
"""
    Output

  Coordinates simulation output: manages diagnostics, drives the writer(s), and buffers the
  physical-space state for diagnostics that need it.

  Writers are passed explicitly, so any combination of `HDF5Writer`, `LocalWriter`, or
  future writers (NetCDF, Zarr, …) can be used together or alone. The simulation loop
  calls only `handle_output!` and never touches the writers directly.

  # Fields
  - `diagnostics`: Array of `Diagnostic`s to evaluate each step.
  - `strides`: Per-diagnostic stride (steps between samples).
  - `state_buffer`: Pre-allocated buffer for the physical-space state.
  - `writers`: Vector of `AbstractWriter`s that receive each sample.
  - `physical_transform`: Function applied to the state after the inverse transform.
  - `transformed`: Flag — `true` if the state has already been transformed this step.

  # Constructor

    Output(prob; writers, physical_transform, storage_limit)

  ## Keyword Arguments
  - `writers`: One or more `AbstractWriter` instances (default: `[HDF5Writer(prob)]`).
  - `physical_transform`: Applied after the inverse spectral transform (default: `identity`).
  - `storage_limit`: Optional total storage cap as a string, e.g. `"2 GB"` (default: `""`).
  
  - `simulation_name`: Name of the top-level HDF5 group (default: `:filename`).
    Supported symbols: `:filename`, `:timestamp`, `:parameters`, or a plain `String`.
  - `resume`: If `true`, open an existing file and validate its attributes (default: `false`).
"""
mutable struct Output{DV<:AbstractArray{<:Diagnostic},UB<:AbstractArray,
    WV<:AbstractVector{<:AbstractWriter},PT<:Function}
    diagnostics::DV
    strides::Vector{Int}
    state_buffer::UB
    simulation_name::String
    writers::WV
    physical_transform::PT
    transformed::Bool
    resume::Bool
end

# TODO add #resume::Bool=false, and #simulation_name::Union{AbstractString,Symbol}=:filename,

function Output(prob::SpectralODEProblem; simulation_name::AbstractString="",
    writers::AbstractVector{<:AbstractWriter}=[HDF5Writer(prob)],
    physical_transform::PT=identity, resume::Bool=false,
    storage_limit::AbstractString="") where {PT<:Function}

    # Prepare physical-space state buffer and sample each diagnostic once at t0
    state, t0 = prepare_initial_state(prob; physical_transform=physical_transform)
    diagnostics, initial_samples = initialize_diagnostics(prob, state, prob.u0_hat, t0)
    strides = determine_strides(initial_samples, prob, storage_limit)

    # Give every writer a chance to allocate its storage
    for writer in writers
        initialize!(writer, prob, diagnostics, initial_samples, t0; resume, simulation_name)
    end

    return Output(diagnostics, strides, state, simulation_name, writers, physical_transform, true, resume)
end

# --------------------------------- Prepare Initial State ----------------------------------

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

# ----------------------------------- Runtime Handling -------------------------------------

"""
    handle_output!(output, step, state, prob, time)

  Called by the integrator after every time step. Removes aliased modes, samples all due
  diagnostics (transforming the state at most once per step), flushes writers, and checks
  for NaN breakdown.
"""
function handle_output!(output::Output, step::Integer, state::AbstractArray,
    prob::SpectralODEProblem, time::Number)

    output.transformed = false

    # Remove modes after each step using user defined function
    prob.remove_modes(state, prob.domain)

    handle_diagnostics!(output, step, state, prob, time)

    for writer in output.writers
        flush!(writer, step)
    end

    # Check if first value is NaN, if one value is NaN the whole Array will turn NaN after FFT
    assert_no_nan(state, time)
end

# --------------------------------- Convenience Constructors ------------------------------

"""
    Output(prob, writer::AbstractWriter; kwargs...)

  Convenience form accepting a single writer without wrapping it in a vector.
"""
Output(prob::SpectralODEProblem, writer::AbstractWriter; kwargs...) =
    Output(prob; writers=[writer], kwargs...)

"""
    resolve_simulation_name(simulation_name, prob, filename)

  Converts the `simulation_name` keyword into a group name string.

  ### Supported symbols
  * `:filename`   — basename of the file without extension.
  * `:timestamp`  — `Dates.now()` formatted as a string.
  * `:parameters` — parameter names and values joined into a string.
"""
function resolve_simulation_name(simulation_name, prob, filename)
    simulation_name == :filename && return first(splitext(basename(filename)))
    simulation_name == :timestamp && return "$(now())"
    simulation_name == :parameters && return parameter_string(prob.p)
    simulation_name isa String && return simulation_name
    error("$simulation_name is not a valid input")
end

# ------------------------------------ Base overloads -------------------------------------

function Base.close(output::Output)
    for writer in output.writers
        close(writer)
    end
end

function Base.show(io::IO, m::MIME"text/plain", output::Output)
    print(io, "Output")
    if output.physical_transform !== identity
        print(io, " (physical_transform: ", nameof(output.physical_transform), ")")
    end
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