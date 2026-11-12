# ------------------------------------------------------------------------------------------
#                                        HDF5 Writer                                        
# ------------------------------------------------------------------------------------------

"""
    HDF5Writer <: AbstractWriter

  Writes diagnostic samples and checkpoints to an HDF5 file. One HDF5 group is created
  per diagnostic inside the simulation group, each holding a `"data"` dataset
  (pre-allocated to the expected number of samples) and a `"time"` dataset for time points.

  Construction is lazy — no file structure is built until `initialize!` is called by
  `Output`. The file itself is opened at construction time to hold the placeholder group.

  # Constructor

    HDF5Writer(; filename, overwrite_existing, flush_interval, h5_kwargs...)

  ## Keyword Arguments
  - `filename`: Path to the HDF5 file (default: random temp name).
  - `overwrite_existing`: If `true`, always overwrites an existing file. Always produces
    a warning to prevent being left on by accident (default: `false`).
  - `flush_interval`: Minutes between forced file flushes (default: `10`).
  - `h5_kwargs...`: Extra keyword arguments forwarded to HDF5 dataset creation
    (e.g. `blosc=3` for compression).
"""
struct HDF5Writer <: AbstractWriter
    simulation::HDF5.Group        # placeholder until initialize! is called
    overwrite_existing::Bool
    flush_interval::Int
    last_flush_time::Ref{DateTime}
    h5_kwargs::NamedTuple
end

function HDF5Writer(; filename::AbstractString=basename(tempname()) * ".h5",
    overwrite_existing::Bool=false, flush_interval::Int=10, h5_kwargs...)

    file = open_or_create_hdf5_file(filename)
    h5_kwargs = merge((blosc=3,), h5_kwargs)
    simulation = placeholder_group(file)
    return HDF5Writer(simulation, overwrite_existing, flush_interval, now(), h5_kwargs)
end

# --------------------------------------- Interface ----------------------------------------

function initialize!(writer::HDF5Writer, prob, diagnostics, initial_samples, strides, t0;
    simulation_name::String="advectra", resume::Bool=false, kwargs...)
    simulation = open_or_create_simulation_group(writer, simulation_name, prob; resume)

    # Transfer HDF5 hid_t and file properties to the placeholder
    setsimulation!(writer, simulation)

    # TODO to be updated when updating the diagnostic interface
    # N_steps = compute_number_of_steps(prob)
    # for (diagnostic, sample, stride) in zip(diagnostics, initial_samples, strides)
    #     if diagnostic.stores_data
    #         N_samples = cld(N_steps, stride) + 1
    #         setup_diagnostic_group(simulation, diagnostic, N_samples, sample, t0;
    #             writer.h5_kwargs)
    #     end
    # end
end

# TODO need to have extent growable!!!
function write_to!(writer::HDF5Writer, diagnostic, idx::Integer, sample, time)
    isinitialized(writer) ||
        error("HDF5Writer has not been initialized — call initialize! before write_to!")
    h5group = writer.simulation[diagnostic.name]
    # TODO better check on ndims
    h5group["data"][fill(:, ndims(sample))..., idx] = sample
    h5group["time"][idx] = time
end

# TODO make less ambiguos
function flush!(writer::HDF5Writer, step)
    isopen(writer) || return
    if now() - writer.last_flush_time[] >= Minute(writer.flush_interval) ||
        step % 1024 == 0
        flush(writer.simulation.file)
        writer.last_flush_time[] = now()
    end
end

Base.close(writer::HDF5Writer) = close(writer.simulation.file)
Base.isopen(writer::HDF5Writer) = HDF5.isvalid(writer.simulation.file)
isinitialized(writer::HDF5Writer) = HDF5.isvalid(writer.simulation)
supports_checkpoint(::HDF5Writer) = true
has_checkpoint(writer::HDF5Writer) = isinitialized(writer) && haskey(writer.simulation, "checkpoint")

"""
    HDF5.file(writer::HDF5Writer)

  Returns the underlying HDF5 file handle, consistent with HDF5.jl's own `HDF5.file(group)`
  convention.
"""
HDF5.file(writer::HDF5Writer) = writer.simulation.file

"""
    setsimulation!(writer::HDF5Writer, simulation::HDF5.Group)

  Transfers ownership of `simulation` into `writer`, replacing the placeholder group.
  Invalidates `simulation` by setting its id to `-1` to prevent double-close.
"""
function setsimulation!(writer::HDF5Writer, simulation::HDF5.Group)
    writer.simulation.id = simulation.id
    writer.simulation.file = simulation.file
    simulation.id = HDF5.API.hid_t(-1)
end

# --------------------------------- Checkpoint Interface -----------------------------------

"""
    store_checkpoint!(writer::HDF5Writer, cache, step, time)

  Dumps all non-`AbstractTableau` fields of `cache` into a `"checkpoint"` group inside
  the simulation group, overwriting any previously stored checkpoint.
"""
function store_checkpoint!(writer::HDF5Writer, cache::AbstractCache, step::Integer,
    time::Number)
    isinitialized(writer) || return
    checkpoint = create_or_open_group(writer.simulation, "checkpoint")
    for key in fieldnames(typeof(cache))
        val = getfield(cache, key)
        # Do not want to bother with storing Tableau, easy to recover
        isa(val, AbstractTableau) && continue
        isa(val, AbstractArray) && (val = Array(val))  # download from GPU if needed
        rewrite_dataset(checkpoint, string(key), val)
    end
    rewrite_dataset(checkpoint, "time", time)
    rewrite_dataset(checkpoint, "step", step)  # TODO: may collide with a cache field named step
end # TODO needs some rethinking

"""
    restore_checkpoint(writer::HDF5Writer, prob, scheme)

  Restores the integrator cache for `scheme` from the checkpoint stored in `writer`.
"""
function restore_checkpoint(writer::HDF5Writer, prob::SpectralODEProblem,
    scheme::AbstractODEAlgorithm)
    # TODO validate that dt and other critical parameters match
    # abs(t_c - (t0 + step_c * Δt_c)) < tol || error("Checkpoint is inconsistent")
    # t0_new == t0 || error("t0 cannot change on resume")
    cache = get_cache(prob, scheme)
    for key in fieldnames(typeof(cache))
        field = getfield(cache, key)
        # Skip restoring Tableau
        isa(field, AbstractTableau) && continue
        data = adapt(typeof(field), read(writer.simulation["checkpoint"], string(key)))
        setproperty!(cache, key, data)
    end
    t = read(writer.simulation["checkpoint"], "time")
    step = read(writer.simulation["checkpoint"], "step")
    return cache, step, t
end

# ---------------------------------- File / Group Setup ------------------------------------

"""
    open_or_create_hdf5_file(filename)

  Resolves the full path for `filename`, creates any missing directories, and opens the
  file in `"cw"` mode (read/write, create if not exists).
"""
function open_or_create_hdf5_file(filename::AbstractString)
    filename = add_h5_if_missing(filename)
    filepath = determine_filepath(filename)
    mkpath(dirname(filepath))
    return h5open(filepath, "cw")
end

"""
    placeholder_group(file::HDF5.File)

  Creates a placeholder `HDF5.Group` for `file` with an invalid id. Used to initialise
  `HDF5Writer` before `initialize!` replaces it with the real simulation group.
"""
placeholder_group(file::HDF5.File) = HDF5.Group(-1, file)

# ----------------------------- Simulation Group Flag Checks -------------------------------

check_contradictory_flags(resume, overwrite_existing) =
    resume && overwrite_existing &&
        error("resume=true and overwrite_existing=true are contradictory. " *
            "Set overwrite_existing=false to resume, or resume=false to overwrite.")

check_ambiguous_existing(resume, overwrite_existing, exists, simulation_name, filepath) =
    !resume && !overwrite_existing && exists &&
        error("Simulation \"$simulation_name\" already exists in \"$filepath\". " *
            "Pass resume=true to continue or overwrite_existing=true to overwrite.")

warn_overwrite_existing(overwrite_existing, exists, simulation_name, filepath) =
    overwrite_existing && @warn "overwrite_existing=true — " *
        (exists ? "overwriting \"$simulation_name\" in \"$filepath\"."
                                 : "no existing file found, starting fresh.")

warn_no_checkpoint(resume, exists, simulation_name, filepath) =
    resume && !exists && @warn "resume=true but \"$simulation_name\" not found " *
        "in \"$filepath\" — starting fresh."

# ------------------------------------------------------------------------------------------

"""
    open_or_create_simulation_group(writer, simulation_name, prob; resume, overwrite_existing)

  Opens or creates the top-level HDF5 group for a simulation run, implementing the
  full resume/overwrite logic. See the `Output` docstring for the full behaviour matrix.
"""
function open_or_create_simulation_group(writer::HDF5Writer, simulation_name::String, prob;
    resume::Bool=false)

    file = HDF5.file(writer)
    filepath = file.filename
    exists = haskey(file, simulation_name)
    overwrite_existing = writer.overwrite_existing

    check_contradictory_flags(resume, overwrite_existing)
    check_ambiguous_existing(resume, overwrite_existing, exists, simulation_name, filepath)
    warn_overwrite_existing(overwrite_existing, exists, simulation_name, filepath)
    warn_no_checkpoint(resume, exists, simulation_name, filepath)

    if resume && exists
        # Existing group: Open it
        simulation = open_group(file, simulation_name)
        validate_resume_attributes(simulation, prob)
        @info "Resume validation successful for simulation: $simulation_name"
        return simulation
    end

    if overwrite_existing && exists
        # Fresh (overwriten) run: delete the existing file and start clean
        close(file)
        rm(filepath)
        file = open_or_create_hdf5_file(filepath)
    end

    # Fresh start: Create group and write attributes
    simulation = create_group(file, simulation_name)
    write_attributes(simulation, prob)
    return simulation
end

# TODO check if new dim of fields, in that case probably should re-create simulation group
"""
    setup_diagnostic_group(simulation, diagnostic, N_samples, sample, t0; h5_kwargs)

  Creates an HDF5 group for `diagnostic` inside `simulation`, pre-allocating datasets for
  `N_samples` time steps and writing the initial sample. On a resume path, extends the
  existing datasets to the new `N_samples`.
"""
function setup_diagnostic_group(simulation, diagnostic, N_samples, sample, t0; h5_kwargs)
    if !haskey(simulation, diagnostic.name)
        # Create diagnostic group
        h5group = create_group(simulation, diagnostic.name)

        # Create dataset to store samples and associated time
        dset = create_dataset(h5group, "data", datatype(eltype(sample)),
            (size(sample)..., typemax(Int64));
            chunk=default_chunk(sample), h5_kwargs...)
        HDF5.set_extent_dims(dset, (size(sample)..., N_samples))
        dset = create_dataset(h5group, "time", datatype(eltype(t0)), (typemax(Int64),);
            chunk=default_chunk(t0), h5_kwargs...)
        HDF5.set_extent_dims(dset, (N_samples,))

        # Add metadata
        write_attribute(h5group, "metadata", diagnostic.metadata)

        # TODO where should this logic be?
        # Store initial sample
        h5group["data"][fill(:, ndims(sample))..., 1] = sample
        h5group["time"][1] = t0

        # Store the initial conditions
        #write_state(simulation, 1, state, t0)
    else
        # Resume path: extend datasets to new size
        h5group = open_group(simulation, diagnostic.name)
        # Extend size of arrays
        # Open dataset
        dset = open_dataset(h5group, "data")
        HDF5.set_extent_dims(dset, (size(sample)..., N_samples))
        dset = open_dataset(h5group, "time")
        HDF5.set_extent_dims(dset, (N_samples,))
    end
end

# -------------------------------------- Attributes ----------------------------------------
"""
    write_attributes(simulation, prob::SpectralODEProblem)
    write_attributes(simulation, domain::AbstractDomain)

  Writes `dt`, domain properties, and physical parameters as HDF5 attributes on the
  simulation group. Also stores the config file that created this `Output`, if available.
"""
function write_attributes(simulation, prob::SpectralODEProblem)
    write_attribute(simulation, "dt", prob.dt)
    write_attributes(simulation, prob.domain)
    for (key, val) in pairs(prob.p)
        write_attribute(simulation, string(key), val)
    end

    # Store the config file used to create the Output struct
    configfile = Base.source_path()
    !isempty(configfile) && open(configfile, "r") do file
        write_attribute(simulation, "config-file", read(file, String))
    end
end

function write_attributes(simulation, domain::AbstractDomain)
    # Construct list of attributes by removing derived attributes
    attributes = setdiff(fieldnames(typeof(domain)),
        (:x, :y, :kx, :ky, :SC, :transforms, :precision, :MemoryType))
    for attribute in attributes
        write_attribute(simulation, string(attribute), getproperty(domain, attribute))
    end
end

"""
    validate_resume_attributes(simulation::HDF5.Group, prob::SpectralODEProblem)

  Checks that critical numerical parameters (Nx, Ny, dt, real_transform) in the existing
  HDF5 group match the current `prob`.

  !!! warning "User Responsibility"
      Advectra only checks a subset of critical numerical parameters (Nx, Ny, dt,
      real_transform). It does **not** perform an exhaustive check of all physical input
      parameters. It is the responsibility of the user to ensure the simulations are
      physically compatible.
"""
function validate_resume_attributes(simulation::HDF5.Group, prob::SpectralODEProblem)
    @warn "Advectra only checks a subset of critical numerical parameters (Nx, Ny, dt, " *
        "real_transform). Advectra does not check if ALL input parameters are identical " *
        "to the run you are resuming from. It is the responsibility of the user to make " *
        "sure the simulations are compatible."

    # List of attributes to check and their corresponding values in 'prob' or 'prob.domain'
    expected_values = Dict("dt" => prob.dt,
        "Nx" => prob.domain.Nx,
        "Ny" => prob.domain.Ny,
        "real_transform" => prob.domain.real_transform)

    for (key, expected) in expected_values
        haskey(attributes(simulation), key) ||
            error("Resume Error: Attribute '$key' missing from existing simulation group.")
        actual = read_attribute(simulation, key)
        actual != expected &&
            error("Resume Mismatch: '$key' was $(actual) in file, but is $(expected) now.")
    end
    return true
end

# ------------------------------------- HDF5 Helpers ---------------------------------------

"""
    create_or_open_group(parent, path; properties...)

  Creates or opens an HDF5 group at `path` within `parent`.
"""
function create_or_open_group(parent::Union{HDF5.File,HDF5.Group}, path::AbstractString;
    properties...)
    haskey(parent, path) ? open_group(parent, path, properties...) :
    create_group(parent, path, properties...)
end

"""
    rewrite_dataset(parent, name, data; pv...)

  Deletes the dataset `name` if it exists, then writes a fresh one.
"""
function rewrite_dataset(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, data;
    pv...)
    haskey(parent, name) && delete_object(parent[name])
    write_dataset(parent, name, data; pv...)
end

"""
    default_chunk(sample; scalar_target_bytes, array_target_bytes)

  Computes a sensible HDF5 chunk shape from `sample`.
"""
function default_chunk(sample; scalar_target_bytes=64 * 1024,
    array_target_bytes=1024 * 1024)
    elsize = sizeof(eltype(sample))
    if ndims(sample) == 0
        n = max(1, round(Int, scalar_target_bytes / elsize))
        return (n,)
    else
        array_bytes = prod(size(sample)) * elsize
        time_chunk = max(1, round(Int, array_target_bytes / array_bytes))
        return (size(sample)..., time_chunk)
    end
end

"""
    add_h5_if_missing(filename)

  Appends `.h5` to `filename` if it has no extension.
"""
add_h5_if_missing(filename::AbstractString) =
    splitext(filename)[end] == "" ? splitext(filename)[1] * ".h5" : filename

# --------------------------------------- Base.show ----------------------------------------

function Base.show(io::IO, ::MIME"text/plain", writer::HDF5Writer)
    print(io, "HDF5Writer → ", writer.simulation.file.filename,
        " (group: ", writer.simulation.id == -1 ? "invalid" : HDF5.name(writer.simulation), ")")
end
