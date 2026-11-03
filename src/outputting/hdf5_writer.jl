# ------------------------------------------------------------------------------------------
#                                        HDF5 Writer
# ------------------------------------------------------------------------------------------

"""
    HDF5Writer <: AbstractWriter

  Writes diagnostic samples to an HDF5 file. One HDF5 group is created per diagnostic,
  each holding a `"data"` dataset (pre-allocated to the expected number of samples) and a
  `"t"` dataset for the corresponding time points.

  # Constructor

    HDF5Writer(; filename, flush_interval, h5_kwargs...)

  ## Keyword Arguments
  - `filename`: Path to the HDF5 file (default: random temp name).
  - `flush_interval`: Minutes between forced file flushes (default: `10`).
  - `h5_kwargs...`: Extra keyword arguments forwarded to HDF5 dataset creation (e.g. compression).
"""
mutable struct HDF5Writer <: AbstractWriter
    simulation::HDF5.Group
    h5_kwargs::NamedTuple
    flush_interval::Int
    last_flush_time::DateTime
end

function HDF5Writer(; filename::AbstractString=basename(tempname()) * ".h5",
    flush_interval::Int=10,
    h5_kwargs...)

    file = open_hdf5_file(filename)
    h5_kwargs = merge((blosc=3,), h5_kwargs)
    simulation = placeholder_group(file)
    return HDF5Writer(simulation, h5_kwargs, flush_interval, now())
end

# ------------------------------------ Interface ------------------------------------------

function initialize!(writer::HDF5Writer, prob, diagnostics, initial_samples, strides, t0)

    # Setup simulation
    simulation = open_or_create_simulation_group(filename, simulation_name, prob;
        h5_kwargs=h5_kwargs, resume=resume)

    # Transfer HDF5 hid_t and file properties to the placeholder
    writer.simulation.id = simulation.id
    writer.simulation.file = simulation.file

    N_steps = compute_number_of_steps(prob)
    for (diagnostic, sample, stride) in zip(diagnostics, initial_samples, strides)
        if diagnostic.stores_data
            N_samples = cld(N_steps, stride) + 1
            setup_diagnostic_group(writer.simulation, diagnostic, N_samples, sample, t0;
                writer.h5_kwargs)
        end
    end
end

function write!(writer::HDF5Writer, diagnostic, idx::Integer, sample, time)
    h5group = writer.simulation[diagnostic.name]
    # TODO better check on ndims
    h5group["data"][fill(:, ndims(sample))..., idx] = sample
    h5group["t"][idx] = time
end

function flush!(writer::HDF5Writer, step)
    if now() - writer.last_flush_time >= Minute(writer.flush_interval) ||
        step % 1024 == 0
        flush(writer.simulation.file)
        writer.last_flush_time = now()
    end
end

Base.close(writer::HDF5Writer) = close(writer.simulation.file)

# ----------------------------------- File / Group Setup ----------------------------------

# TODO check if new dim of fields, in that case probably should re-create simulation group
"""
    open_or_create_simulation_group(filename, simulation_name, prob; h5_kwargs, resume)

  Opens or creates the top-level HDF5 group for a simulation. If `resume=false` and the
  file already exists it is deleted and re-created to guarantee a clean start.
"""
function open_or_create_simulation_group(file, simulation_name, prob;
    h5_kwargs=(blosc=3,), resume=false)

    # TODO move to output.jl and rename variable
    #group_name = resolve_simulation_name(simulation_name, prob, filename)

    # TODO remove temporary
    group_name = simulation_name

    if !haskey(file, group_name)
        # Fresh start: Create group and write attributes
        simulation = create_group(file, group_name)
        write_attributes(simulation, prob)
    else
        # Existing group: Open it
        simulation = open_group(file, group_name)
        if resume
            validate_resume_attributes(simulation, prob)
            @info "Resume validation successful for simulation: $group_name"
        else
            # Fresh (non-resume) run: delete the existing file and start clean
            filepath = file.filename
            close(file)
            rm(filepath)
            file = h5open(filepath, "cw")
            simulation = create_group(file, group_name)
            write_attributes(simulation, prob)
        end
    end

    return simulation
end

# TODO check if new dim of fields, in that case probably should re-create simulation group
"""
    setup_diagnostic_group(simulation, diagnostic, N_samples, sample, t0; h5_kwargs)

  Creates an HDF5 group for `diagnostic` inside `simulation`, pre-allocating datasets for
  `N_samples` time steps and writing the initial sample. If the group already exists (resume
  path), the datasets are extended to the new `N_samples`.
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
        dset = create_dataset(h5group, "t", datatype(eltype(t0)), (typemax(Int64),);
            chunk=default_chunk(t0), h5_kwargs...)
        HDF5.set_extent_dims(dset, (N_samples,))

        # Add metadata
        write_attribute(h5group, "metadata", diagnostic.metadata)

        # TODO where should this logic be?
        # Store initial sample
        h5group["data"][fill(:, ndims(sample))..., 1] = sample
        h5group["t"][1] = t0

        # Store the initial conditions
        #write_state(simulation, 1, state, t0)
    else
        # Resume path: extend datasets to new size
        h5group = open_group(simulation, diagnostic.name)
        # Extend size of arrays
        # Open dataset
        dset = open_dataset(h5group, "data")
        HDF5.set_extent_dims(dset, (size(sample)..., N_samples))
        dset = open_dataset(h5group, "t")
        HDF5.set_extent_dims(dset, (N_samples,))
    end
end

# -------------------------------------- Attributes ---------------------------------------

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
        if !haskey(attributes(simulation), key)
            error("Resume Error: Attribute '$key' missing from existing simulation group.")
        end
        actual = read_attribute(simulation, key)
        if actual != expected
            error("Resume Mismatch: '$key' was $(actual) in file, but is $(expected) now.")
        end
    end
    return true
end

# --------------------------------------- HDF5 Helpers ------------------------------------

function open_hdf5_file(filename)
    filename = add_h5_if_missing(filename)
    filepath = determine_filepath(filename)
    mkpath(dirname(filepath))
    return h5open(filepath, "cw") # Open in "cw" mode (Read/Write, create if not exists)
end

placeholder_group(file::HDF5.File) = HDF5.Group(-1, file)

"""
    create_or_open_group(parent, path; properties...)

  Creates or opens an HDF5 group at `path` within `parent`.
"""
function create_or_open_group(parent::Union{HDF5.File,HDF5.Group}, path::AbstractString;
    properties...)
    haskey(parent, path) ? open_group(parent, path, properties...) :
    create_group(parent, path, properties...)
end

# TODO remove unused method?
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

# TODO move to output.jl
"""
    determine_filepath(filename)

  Resolves the full path for `filename`. When no directory is given, places the file in
  an `output/` subdirectory relative to the source file (interactive) or `pwd()` (script).
"""
function determine_filepath(filename)
    isempty(dirname(filename)) || return filename
    isinteractive() ? joinpath(Base.source_dir(), "output", filename) :
    joinpath(pwd(), "output", filename)
end

# ----------------------------------- Base.show -------------------------------------------

# TODO fix for placeholder_group
function Base.show(io::IO, ::MIME"text/plain", writer::HDF5Writer)
    #print(io, "HDF5Writer → ", writer.simulation.file.filename,
    #    " (group: ", HDF5.name(writer.simulation), ")")
end
