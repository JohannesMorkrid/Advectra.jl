# ------------------------------------------------------------------------------------------
#                                         Outputer                                          
# ------------------------------------------------------------------------------------------

"""
    Output{DV, U, UB, T, S, PT, K}

  A mutable struct for managing simulation output, diagnostics, and storage.
  
  # Fields
  - `stride::Int`: Number of steps between stored samples.
  - `diagnostics::DV`: Array of `Diagnostic`s.
  - `u::U`: Array for storing simulation states in memory.
  - `state_buffer::UB`: Buffer for storing physical space state.
  - `t::T`: Array for storing time points.
  - `simulation::S`: HDF5 group or `Nothing`, representing the simulation storage.
  - `physical_transform::PT`: Function to transform the state in physical space.
  - `store_hdf::Bool`: Whether to store output in HDF5 format.
  - `store_locally::Bool`: Whether to store output in memory.
  - `transformed::Bool`: Indicates if the state has been transformed to physical space.
  - `h5_kwargs::K`: Named tuple of keyword arguments for HDF5 storage.
  - `resume::Bool`: Whether to resume from an existing output file.
  
  # Constructor
  
    Output(prob::SOP; filename::FN, diagnostics::DV, stride::Integer, \
    physical_transform::PT, simulation_name::SN, store_hdf::Bool, store_locally::Bool, \
    storage_limit::AbstractString, h5_kwargs...)

  Creates an `Output` object for a given `SpectralODEProblem`` `prob`. Handles setup for \
  HDF5 and local storage, initializes diagnostics, and manages sampling stride.
  
  ## Keyword Arguments
  - `filename`: Name of the HDF5 file for output (default: random temporary name).
  - `diagnostics`: Array of diagnostics to compute (default: `DEFAULT_DIAGNOSTICS`).
  - `stride`: Number of steps between samples (default: -1 (lets the program decide)).
  - `physical_transform`: Function to transform state in physical space (default: `identity`).
  - `simulation_name`: Name for the simulation group in HDF5 (default: `:timestamp`).
  - `store_hdf`: Store output in HDF5 file (default: `true`).
  - `store_locally`: Store output in memory (default: `true`).
  - `storage_limit`: Limit for field storage (default: empty string).
  - `h5_kwargs...`: Additional keyword arguments for HDF5 storage (merged with defaults).
  - `resume`: Resume from existing output file (default: `false`).

  # Usage

  Create an `Output` object to manage simulation results, diagnostics, and storage options 
  for a given spectral ODE problem.
"""
mutable struct Output{DV<:AbstractArray{<:Diagnostic},UB<:AbstractArray,T<:AbstractArray,
                      S<:Union{HDF5.Group,Nothing},PT<:Function,
                      K<:NamedTuple}
    diagnostics::DV
    strides::Vector{Int}
    state_buffer::UB
    t::T # TODO replace with simulationDict?
    simulation::S
    physical_transform::PT
    store_hdf::Bool
    store_locally::Bool
    transformed::Bool
    h5_kwargs::K #Possibly also called a filter
    flush_interval::Int
    last_flush_time::DateTime
    resume::Bool

    function Output(prob::SOP; filename::FN=basename(tempname()) * ".h5",
                    physical_transform::PT=identity,
                    simulation_name::SN=:filename,
                    store_hdf::Bool=true,
                    store_locally::Bool=true,
                    storage_limit::AbstractString="",
                    flush_interval::Int=10,
                    resume::Bool=false,
                    h5_kwargs...) where {SOP<:SpectralODEProblem,
                                         FN<:AbstractString,PT<:Function,
                                         SN<:Union{AbstractString,Symbol}}

        # Prepare initial state
        state, t0 = prepare_initial_state(prob; physical_transform=physical_transform)

        # Build diagnostics and do initial sampling pass
        diagnostics, initial_samples = initialize_diagnostics(prob, state, prob.u0_hat, t0)

        #built_diagnostics = Diagnostic[] # Currently disabled, todo re-enable
        strides = determine_strides(initial_samples, prob, storage_limit)

        # Merge h5_kwargs with default kwargs
        h5_kwargs = merge((blosc=3,), h5_kwargs)

        # Setup HDF5 storage if wanted
        simulation = setup_hdf5_storage(prob, t0;
                                        filename=filename,
                                        simulation_name=simulation_name,
                                        diagnostics=diagnostics,
                                        initial_samples=initial_samples,
                                        strides=strides,
                                        store_hdf=store_hdf,
                                        h5_kwargs=h5_kwargs,
                                        resume=resume)

        # Setup local (in memory) storage if wanted
        u, t = setup_local_storage(state, t0; store_locally=false) # Currently disabled TODO re-enable

        # Create output
        new{typeof(diagnostics),typeof(state),typeof(t),typeof(simulation),
            typeof(physical_transform),typeof(h5_kwargs)}(diagnostics, strides, state, t,
                                                          simulation, physical_transform,
                                                          store_hdf, store_locally, true,
                                                          h5_kwargs, flush_interval, now(),
                                                          resume)
    end
end

# --------------------------------- Prepare Initial State ----------------------------------

"""
    prepare_initial_state(prob; physical_transform=identity)

  Prepares the initial state by applying the user defined `physical_transform`, if any, to 
  a copy of the initial condition stored in `prob`, returning it alongside the initial time. 
"""
function prepare_initial_state(prob; physical_transform=identity)
    state = copy(prob.u0)
    physical_transform(state)
    t0 = first(prob.tspan)
    return state, t0
end

# -------------------------------------- HDF5 Setup ----------------------------------------

function setup_hdf5_storage(prob, t0;
                            filename,
                            simulation_name,
                            diagnostics,
                            initial_samples,
                            strides,
                            store_hdf=true,
                            h5_kwargs=(blosc=3,),
                            resume)
    simulation = setup_simulation_group(filename, simulation_name, prob; resume=resume,
                                        store_hdf=store_hdf, h5_kwargs=h5_kwargs)

    if !resume
        # if file already exists it will be deleted and must be created again
        rm(simulation.file.filename)
        simulation = setup_simulation_group(filename, simulation_name, prob; resume=false,
                                            store_hdf=store_hdf, h5_kwargs=h5_kwargs)
    end

    if !isnothing(simulation)
        N_steps = compute_number_of_steps(prob)
        for (diagnostic, sample, stride) in zip(diagnostics, initial_samples, strides)
            if diagnostic.stores_data
                N_samples = cld(N_steps, stride) + 1
                setup_diagnostic_group(simulation, diagnostic, N_samples, sample, t0;
                                       h5_kwargs)
            end
        end
    end

    return simulation
end

# TODO check if new dim of fields, in that case probably should re-create simulation group
"""
    setup_hdf5_storage(filename, simulation_name, N_samples::Int, state, prob, t0; 
                       store_hdf=store_hdf, h5_kwargs=h5_kwargs)

  Creates a *HDF5* file, if not existing, and writes a group with `simulation_name` to it, 
  refered to as a `simulation` group. If the simulation group does not exists, the `h5_kwargs` 
  are applied to the `"fields"` and `"t"` datasets. The opened `simulation` is returned.
"""
function setup_simulation_group(filename, simulation_name, prob; store_hdf=true,
                                h5_kwargs=(blosc=3,), resume=false)
    if store_hdf
        filename = add_h5_if_missing(filename)
        filepath = determine_filepath(filename)
        mkpath(dirname(filepath))

        # Open in "cw" mode (Read/Write, create if not exists)
        file = h5open(filepath, "cw")
        group_name = handle_simulation_name(simulation_name, prob, filename)

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
            end
        end

        return simulation
    else
        return nothing
    end
end

function determine_filepath(filename)
    if !isempty(dirname(filename))
        return filename
    end

    if isinteractive()
        return joinpath(Base.source_dir(), "output", filename)
    else
        joinpath(pwd(), "output", filename)
    end
end
# ------------------------------------- HDF5 Helpers ---------------------------------------

"""
    add_h5_if_missing(filename::AbstractString)

  Makes sure the filename has an extension, if not the `.h5` extension is added. 
"""
function add_h5_if_missing(filename::AbstractString)
    splitext(filename)[end] == "" ? splitext(filename)[1] * ".h5" : filename
end

"""
    handle_simulation_name(simulation_name)

  Creates a `simulation_name` string based on the users input. 
  ### Supported symbols:
  * `:filename` creates a string with the basename of the file without the extension.
  * `:timestamp` creates a timestamp string using `Dates.now()`.
  * `:parameters` creates a string with the parameter names and values.
"""
function handle_simulation_name(simulation_name, prob, filename)
    # Handle simulation_name

    if simulation_name == :filename
        return first(splitext(basename(filename)))
    elseif simulation_name == :timestamp
        return "$(now())"
    elseif simulation_name == :parameters
        return parameter_string(prob.p)
    elseif simulation_name isa String
        return simulation_name
    else
        error("$simulation_name is not a valid input")
    end
end

"""
    write_attributes(simulation, prob::SpectralODEProblem)
    write_attributes(simulation, domain::AbstractDomain)

  Writes the esential attributes of the container to the simulation group `simulation`. The 
  `SpectralODEProblem` also writes the `domain` properties.
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

# -------------------------------- HDF5 Diagnostic Storage ---------------------------------

"""
    setup_simulation_group(file, simulation_name, N_samples, state, prob, t0; h5_kwargs)

  Creates a *HDF5* group with `simulation_name` (a "simulation"), and allocates the correct 
  sizes based on `N_samples` with the fields being chunked with additinal `h5_kwargs` applied.
  In addition the inital condition is written along with the attributes of the `prob`.
"""
function setup_diagnostic_group(simulation, diagnostic, N_samples, sample, t0; h5_kwargs)
    if !haskey(simulation, diagnostic.name)
        # Create diagnostic group
        h5group = create_group(simulation, diagnostic.name)

        # Create dataset to store samples and associated time
        dset = create_dataset(h5group, "data", datatype(eltype(sample)),
                              (size(sample)..., typemax(Int64)); chunk=(size(sample)..., 1),
                              h5_kwargs...)
        HDF5.set_extent_dims(dset, (size(sample)..., N_samples))
        dset = create_dataset(h5group, "t", datatype(eltype(t0)), (typemax(Int64),);
                              chunk=(1,), h5_kwargs...)
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
        h5group = open_group(simulation, diagnostic.name)
        # Extend size of arrays
        # Open dataset
        dset = open_dataset(h5group, "data")
        HDF5.set_extent_dims(dset, (size(sample)..., N_samples))
        dset = open_dataset(h5group, "t")
        HDF5.set_extent_dims(dset, (N_samples,))
    end
end

"""
    write_to(h5group::HDF5.Group, idx, data, time)

  Writes the `data` at time `time` to the  `h5group` (HDF5).
"""
function write_to(h5group::HDF5.Group, idx, data, time)
    # TODO better check on ndims
    h5group["data"][fill(:, ndims(data))..., idx] = data
    h5group["t"][idx] = time
end

# ---------------------------------- Setup Local Storage -----------------------------------

# TODO clean-up here later, make local storage a Dict
"""
    setup_local_storage(state, t0, N_samples; store_locally=store_locally)
  
  Allocates vectors in memory for storing the fields alongside the time if the user wants it,
  otherwise empty vectors are returned.
"""
function setup_local_storage(state, t0; store_locally=store_locally)
    if store_locally
        # Allocate local memory for fields
        u = [zero(state) for _ in 1:N_samples]
        u[1] .= state
        t = zeros(N_samples)
        t[1] = t0
    else
        u, t = [], []
    end

    return u, t
end

function setup_local_key()
    # Allocate arrays
    diagnostic.data = [zero(id) for _ in 1:N] #Vector{typeof(id)}(undef, N)
    diagnostic.t = zeros(N)

    # Store intial diagnostic
    if isa(id, AbstractArray)
        diagnostic.data[1] .= id
    else
        diagnostic.data[1] = copy(id)
    end
    diagnostic.t[1] = first(prob.tspan)
end

"""
    write_local_state(output, idx, u, t)

  Writes the state `u` at time `t` to the local storage in the `Output` struct.
"""
function write_local_state(output::Output, idx, u, t) # (diagnostic::Diagnostic, idx, data, t)
    output.u[idx] .= u
    output.t[idx] = t
    if isa(data, AbstractArray)
        output.diagnostic.data[idx] .= data
    else
        diagnostic.data[idx] = copy(data)
    end
    diagnostic.t[idx] = t
end

# ---------------------------- Stride And Storage Size Related -----------------------------

"""
    prepare_sampling(stride::Int, storage_limit, prob::SpectralODEProblem)

  Prepares the sampling strategy for storing simulation states based on the desired stride 
  and storage limit. Determines the number of samples to store and the stride between samples.
  
  #### Arguments
  - `stride::Int`: The proposed stride between samples. If set to `-1`, the function 
  will automatically recommend an appropriate stride based on the storage limit.
  - `storage_limit`: The maximum allowed storage size for the output of fields, as a 
  string (e.g., `"100 MB"`). The limit does not affect the storage size of the Diagnostics! 
  If empty, no storage constraint is applied.
  - `prob::SpectralODEProblem`: The `SpectralODEProblem`` containing the size of the fields.
    
  #### Notes
  - If both `stride` and `storage_limit` are unspecified, all steps are recorded.
  - The function validates and adjusts `stride` to ensure it is within feasible bounds.
  - Issues a warning if the last step has a different stride than the rest.
"""
function determine_sampling_strategy(sample, stride::Int, storage_limit,
                                     prob::SpectralODEProblem; context="")
    # Compute total number of simulation steps
    N_steps = compute_number_of_steps(prob)

    # Handle all the different scenarios, stride = -1 => let the program decide
    if !isempty(storage_limit)
        storage_bytes = parse_storage_limit(storage_limit)
        if stride == -1
            stride = recommend_stride(storage_bytes, N_steps, sample; context=context)
        else
            check_storage_size(storage_bytes, N_steps, stride, sample; context=context)
        end
    elseif stride == -1
        stride = 1
    end

    # Validate stride, might change stride if too large
    stride = validate_stride(N_steps, stride; context=context)

    # Compute number of data points to record and warn if last step has differnt step size
    N_samples = cld(N_steps, stride) + 1

    return N_samples, stride
end

"""
    compute_number_of_steps(prob::SpectralODEProblem)
  
  Computes the number of time steps the integrator needs to solve the problem `prob`. The 
  last step might be fractional and is therefore rounded up due to the fixed `prob.dt`.
"""
function compute_number_of_steps(prob::SpectralODEProblem)
    ceil(Int, (last(prob.tspan) - first(prob.tspan)) / prob.dt)
end

"""
    recommend_stride(storage_limit::Int, N_steps::Int, sample::AbstractArray; context="")
  
  Recommends the closest divisor to the minimum stride needed to fullfil the `storage_limit`.
  If the `storage_limit` is too strict an error is thrown, which informs the user of the 
  minimum limit.
"""
function recommend_stride(storage_limit::Int, N_steps::Int, sample::AbstractArray;
                          context="")
    field_bytes = length(sample) * sizeof(eltype(sample))
    # Determine how many fields can be fully stored
    max_samples = storage_limit ÷ field_bytes

    # At least two should be stored, start and end
    if max_samples < 2
        throw(ArgumentError(context * "The storage limit ($(format_bytes(storage_limit))) \
          is too small. The sample alone requires $(format_bytes(field_bytes)), however at \
          least two samples are required (minimum limit: $(format_bytes(2*field_bytes)))."))
    end

    # Compute minimum stride to achieve the max number of samples
    min_stride = ceil(Int, N_steps / (max_samples - 1))

    # Picks closest divisor that does not exceed storage limit, while N_steps ≥ min_stride
    recommended_stride = next_divisor(N_steps, min_stride)

    return recommended_stride
end

# Edge case
recommend_stride(storage_limit, N_steps, sample::Nothing; context="") = 1

"""
    check_storage_size(storage_limit::Int, N_steps::Int, stride::Int, sample; context="")
  
  Checks that the needed storage does not exceed the `storage_limit`, otherwise an error is 
  thrown, which recommends the minimum divisor satisfying the storage limit. In addition the
   error checks of `recommend_stride` are performed, which may trigger before the storage check.
"""
function check_storage_size(storage_limit::Int, N_steps::Int, stride::Int, sample;
                            context="")
    min_stride = recommend_stride(storage_limit, N_steps, sample; context=context)

    storage_need = compute_storage_need(N_steps, stride, sample; context=context)
    if storage_need > storage_limit
        throw(ArgumentError(context * "The total output requires \
                              $(format_bytes(storage_need)), which exceeds the storage \
                              limit of $(format_bytes(storage_limit)). Consider increasing \
                              the `stride` (minimum recommended: $min_stride) or the `storage_limit`."))
    end
end

"""
    compute_storage_need(N_steps::Int, stride::Int, sample::AbstractArray; context="")

  Computes the storage needed to store `N_steps÷stride`samples with `sizeof(sample)`. 
"""
function compute_storage_need(N_steps::Int, stride::Int, sample::AbstractArray; context="")
    stride < 1 ? throw(ArgumentError(context * "stride must be ≥ 1, got $stride")) : nothing
    (cld(N_steps, stride) + 1) * length(sample) * sizeof(eltype(sample))
end

function compute_storage_need(N_steps, stride, sample::Number; context="")
    stride < 1 ? throw(ArgumentError(context * "stride must be ≥ 1, got $stride")) : nothing
    (cld(N_steps, stride) + 1) * sizeof(eltype(sample))
end

# Edge case
compute_storage_need(N_steps, stride, sample::Nothing; context="") = 0

"""
    validate_stride(N_steps::Int, stride::Int; context="")
  
    Validates and adjusts the `stride`. If `stride`:

  - exceeds `N_steps`, it is set to `N_steps` and a warning is issued.
  - is less than 1 an `ArgumentError` is thrown.
  - does not evenly divide `N_steps`, a warning is issued and a suggested divisor is provided.

  Returns the validated (and possibly adjusted) `stride`.
"""
function validate_stride(N_steps::Int, stride::Int; context="")
    if N_steps < stride
        @warn context * "stride ($stride) exceeds total steps ($N_steps). \
                 Adjusting to stride = N_steps ($N_steps)."
        stride = N_steps
    end

    stride < 1 ? throw(ArgumentError(context * "stride must be ≥ 1, got $stride")) : nothing

    if N_steps % stride != 0
        suggestion = nearest_divisor(N_steps, stride)
        @warn context * "stride ($stride) does not evenly divide N_steps ($N_steps). The \
          final output interval will be shorter. Consider using stride = $suggestion instead."
    end

    return stride
end

# ----------------------------- Storage Size Helper Functions ------------------------------

"""
    parse_storage_limit(limit::Integer)
    parse_storage_limit(limit::AbstractString)

  Parses a storage limit specified either as an integer (number of bytes) or as a string with 
  optional units (defaults to B). Supports both decimal and binary units up to "EiB", with a 
  decimal number of bytes. For example: `"10MB"`, `"5.2 GiB"`, `"1024"`.
"""
parse_storage_limit(limit::Integer) = limit

function parse_storage_limit(limit::AbstractString)
    # Units dictionary
    units = Dict("B" => 1, "KB" => 10^3, "KIB" => 1024, "MB" => 10^6, "MIB" => 1024^2,
                 "GB" => 10^9, "GIB" => 1024^3, "TB" => 10^12, "TIB" => 1024^4,
                 "PB" => 10^15,
                 "PIB" => 1024^5, "EB" => Int128(10^18), "EIB" => Int128(1024^6))

    # Match numeric value and optional unit
    # \s* matches 0 or more white spaces (\s) at the beginning ^ and end $
    # \d matches digits + -> 1 or more, * 0 or more
    # (?:\.\d) non-capturing group starting with a dot \., ? means matches zero or one
    # [a-zA-Z] matches all combinations of roman letters
    m = match(r"^\s*(\d*(?:\.\d*)?)\s*([a-zA-Z]*)\s*$", limit)
    if m === nothing
        throw(ArgumentError("Invalid storage limit format: $limit"))
    end

    # Get value and unit from Regex
    value = parse(Float64, m.captures[1])
    unit = uppercase(m.captures[2])
    unit = unit == "" ? "B" : unit # Converts no unit to mean byte

    if !haskey(units, unit)
        throw(ArgumentError("Unknown storage unit: $unit"))
    end

    multiplier = units[unit]
    return round(typeof(multiplier), value * multiplier)
end

"""
    format_bytes(nbytes::Integer)
  
  Formats the number of bytes into a "human readable" format. For example: `1032` ⇒ `"1.03 KB"`.
"""
function format_bytes(nbytes::Integer)
    units = ["B", "KB", "MB", "GB", "TB", "PB", "EB"]
    i = 1
    f = float(nbytes)
    while f >= 1000 && i < length(units)
        f /= 1000
        i += 1
    end
    return string(round(f; digits=2), " ", units[i])
end

"""
    next_divisor(N, target)

  Finds the smallest divisor of `N` that is ≥ `target`. If none exists (i.e. `target > N`), \
  `N` is returned.
"""
function next_divisor(N::Int, target::Int)
    recommended = N
    for d in 1:floor(Int, sqrt(N))
        if N % d == 0
            for c in (d, N ÷ d)
                if c >= target && c < recommended
                    recommended = c
                end
            end
        end
    end
    return recommended
end

"""
    nearest_divisor(N::Int, target::Int)
  Finds the nearest divisor of `N` to the target. If `N` is prime `1` is returned.
"""
function nearest_divisor(N::Int, target::Int)
    closest, smallest = 1, abs(target - 1)
    for d in 1:floor(Int, sqrt(N))
        if N % d == 0
            for c in (d, N ÷ d)
                dist = abs(c - target)
                if dist < smallest || (dist == smallest && c < closest)
                    closest, smallest = c, dist
                end
            end
        end
    end
    return closest
end

# ------------------------------- Diagnostic Initialization --------------------------------

"""
    initialize_diagnostics(prob, state, state_hat, t0)

  Build diagnostics from `prob.diagnostic_recipes` and sample.
"""
function initialize_diagnostics(prob::SpectralODEProblem, state, state_hat, t0)
    diagnostics = Diagnostic[]
    initial_samples = []

    for recipe in prob.diagnostic_recipes
        # Build diagnostic
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

"""
"""
function determine_strides(initial_samples, prob::SpectralODEProblem, total_storage_limit)
    strides = Int[]
    storage_requirements = Int[]
    total_storage_requirement = 0

    for (recipe, sample) in zip(prob.diagnostic_recipes, initial_samples)
        @unpack stride, storage_limit = recipe
        context = "($(recipe.method)) "
        # Determine the number of samples to be stored and the stride distance
        N_samples,
        stride = determine_sampling_strategy(sample, stride, storage_limit, prob;
                                             context=context)
        # Determine the needed storage
        storage_requirement = compute_storage_need(N_samples, stride, sample; context)
        # Accumulate
        total_storage_requirement += storage_requirement
        push!(strides, stride)
        push!(storage_requirements, storage_requirement)
    end

    # Nice to have
    println("Estimated filesize: ", format_bytes(total_storage_requirement))

    # Compare cumulative storage need to Output storage need
    if !isempty(total_storage_limit) &&
       parse_storage_limit(total_storage_limit) < total_storage_requirement
        # TODO add nice error message showing what each diagnostic requires, so the user can make up their mind
        error("The Output requires $(format_bytes(total_storage_requirement)), which is \
        more than the storage_limit: $total_storage_limit.")
    end
    return strides
end

""" 
  # In memory storage
  setup_local_storage
"""

# ---------------------------------- Handling Of Output ------------------------------------

"""
    handle_output!(output::O, step::Integer, u::T, prob::SOP, t::N)

  Handles output operations for a simulation step, including state storage, diagnostics 
  sampling, and mode removal. Ensures the spectral state is only transformed once per step. 
  In addition a check is performed to detect breakdowns, to throw an error. 
"""
function handle_output!(output::O, step::Integer, state::T, prob::SOP,
                        time::N) where {O<:Output,T<:AbstractArray,SOP<:
                                        SpectralODEProblem,N<:Number}

    # Keeps track such that state only transformed once
    output.transformed = false

    # Remove modes after each step using user defined function
    prob.remove_modes(state, prob.domain)

    # Handle diagnostics
    maybe_sample_diagnostics!(output, step, state, prob, time)

    # Handle flushing of file
    maybe_flush!(output)

    # Check if first value is NaN, if one value is NaN the whole Array will turn NaN after FFT
    assert_no_nan(state, time)
end

"""
    maybe_sample_diagnostics!(output, step::Integer, u, prob, t)

  Iterates through the list of diagnostics (`output.diagnostics`) and determines whether or
  not to sample the diagnostic.
"""
function maybe_sample_diagnostics!(output, step::Integer, state, prob, time)
    # Handle diagnostics
    for (diagnostic, stride) in zip(output.diagnostics, output.strides)
        if step % stride == 0
            # Calculate index
            idx = step ÷ stride + 1
            sample_diagnostic!(output, diagnostic, idx, state, prob, time)
        end
    end
end

"""
   TODO write actuall string: The spectral state `u` is transformed to the real
  state `U`, with the user defined `physical_transform` applied, before being stored.
"""
function sample_diagnostic!(output, diagnostic, idx::Integer, state, prob, time)
    # Check if diagnostic assumes physical field and transform if not yet done
    if !diagnostic.assumes_spectral_state && !output.transformed
        # Transform state (updates state_buffer)
        transform_state!(output, state, get_bwd(prob.domain))
    end

    # Apply diagnostic to the correct input (either physical or spectral state)
    input = diagnostic.assumes_spectral_state ? state : output.state_buffer
    sample = diagnostic(input, prob, time)

    # Store sample
    store_diagnostic!(output, diagnostic, idx, sample, time)
end

"""
    store_diagnostic!(output, step::Integer, sample, time)

  Stores sample to *HDF5* file and memory, depending on the state of the `output.store_hdf` 
  and `output.store_locally` respectively. The index is computed based on the step.
"""
function store_diagnostic!(output, diagnostic, idx::Integer, sample, time)
    if !isnothing(sample)
        if output.store_hdf
            write_to(output.simulation[diagnostic.name], idx, sample, time)
        end
        if output.store_locally
            #output.store_locally ? write_local_data(diagnostic, idx, sample, time) : nothing
        end
    end
end

"""
    transform_state!(output, u, p)

  Transforms spectral coefficients `u_hat` into real fields by applying `spectral_transform!` 
  to the `output.state_buffer` buffer. The user defined `physical_transform` is also applied to 
  the buffer, and the `output.transformed` flag is updated to not transform same field twice.
"""
function transform_state!(output, state_hat, plan)
    spectral_transform!(output.state_buffer, plan, state_hat)
    output.physical_transform(output.state_buffer)
    output.transformed = true
end

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

function maybe_flush!(output)
    # Time based flushing
    if now() - output.last_flush_time >= Minute(output.flush_interval)
        output.store_hdf ? flush(output.simulation.file) : nothing
        output.last_flush_time = now()
    end
end
# -------------------------------------- Checkpoint ----------------------------------------

"""
    save_checkpoint!(output::O, cache::C, step::Integer, t::N) where {O<:Output,
    C<:AbstractCache,N<:Number}

  Creates or opens a `checkpoint` and stors the `cache`at time `t` corresponding to step=`step`.
"""
function save_checkpoint!(output::O, cache::C, step::Integer,
                          time::N) where {O<:Output,
                                          C<:AbstractCache,N<:Number}
    if output.store_hdf
        # Create or open a h5group for the checkpoint
        checkpoint = create_or_open_group(output.simulation, "checkpoint")
        # Store checkpoint
        store_checkpoint!(checkpoint, cache, step, time)
    end
end

"""
    store_checkpoint!(checkpoint::G, cache::C, step::Integer, t::N) where {G<:HDF5.Group,
    C<:AbstractCache,N<:Number}

  Stores checkpoint by dumping the fields of the `cache` that are not of type 
  `AbstractTableau`. This overwrites the previously stored checkpoint.
""" # Perhaps one could look into HDF5 compound types in the future
function store_checkpoint!(checkpoint::G, cache::C, step::Integer,
                           time::N) where {G<:HDF5.Group,
                                           C<:AbstractCache,N<:Number}

    # Get all the attributes of cache
    keys = fieldnames(typeof(cache))

    # Dump cache
    for key in keys
        val = getfield(cache, key)
        # Do not want to bother with storing Tableau, easy to recover
        isa(val, AbstractTableau) ? continue : nothing

        # Convert to Array for easy handling, also downloads from GPU in case of CPUArray
        isa(val, AbstractArray) ? val = Array(val) : nothing

        # Backup the field
        rewrite_dataset(checkpoint, string(key), val)
    end

    # Backup step and time
    rewrite_dataset(checkpoint, "time", time)
    rewrite_dataset(checkpoint, "step", step) # TODO step might collide with cache step
end

"""
    restore_checkpoint(simulation::HDF5.Group, prob::SOP, scheme::SA) where {
    SOP<:SpectralODEProblem,SA<:AbstractODEAlgorithm}

  Restores Cache for `scheme` from checkpoint stored in `simulation`.
"""
function restore_checkpoint(simulation::HDF5.Group, prob::SOP,
                            scheme::SA) where {SOP<:SpectralODEProblem,
                                               SA<:AbstractODEAlgorithm}

    #validate_simulation_group TODO check that dt remains the same, and other parameters

    # Create cache container and get the fieldnames
    cache = get_cache(prob, scheme)
    keys = fieldnames(typeof(cache))

    # Restore cache
    for key in keys
        field = getfield(cache, key)
        # Skip restoring Tableau
        isa(field, AbstractTableau) ? continue : nothing

        # Adapt the data to the same type as the field
        data = adapt(typeof(field), read(simulation["checkpoint"], string(key)))
        setproperty!(cache, key, data)
    end

    return cache
end

# --------------------------------------- Utilities ----------------------------------------

# TODO detangle this mess
"""
    parameter_string(parameters::AbstractDict)
    parameter_string(parameters::AbstractDict)

  Creates a human readable string listing the parameters with values, for example: 
  `(κ=0.1,σ=0.02)` ⇒ `"κ=0.1, σ=0.02"`.
"""
function parameter_string(parameters::P) where {P<:AbstractDict}
    tmp = [string(key, "=", value) for (key, value) in sort(collect(parameters))]
    join(tmp, ", ")
end

function parameter_string(parameters::P) where {P<:NamedTuple}
    tmp = [string(key, "=", value) for (key, value) in sort(collect(pairs(parameters)))]
    join(tmp, ", ")
end

"""
    assert_no_nan(u::AbstractArray, t)
    assert_no_nan(u::AbstractGPUArray, t)

  Checks if the first entry in `u` is `NaN`, if so a breakdown occured and an error is thrown.
"""
assert_no_nan(state::AbstractArray, time) =
    if isnan(state[1])
        error("Breakdown occured at t=$time")
    end

function assert_no_nan(state::AbstractGPUArray, time)
    @allowscalar isnan(state[1]) ? error("Breakdown occured at t=$time") : nothing
end

"""
    create_or_open_group(parent::Union{HDF5.File,HDF5.Group}, path::AbstractString; properties...)

  Creates or opens an HDF5 group at the specified `path` within the given `parent` (which can 
  be an `HDF5.File` or `HDF5.Group`). If the group at `path` exists, it is opened using 
  `open_group`; otherwise, a new group is created using `create_group`. Additional keyword 
  arguments (`properties...`) are passed onto the group creation or opening functions.
"""
function create_or_open_group(parent::Union{HDF5.File,HDF5.Group}, path::AbstractString;
                              properties...)
    if haskey(parent, path)
        return open_group(parent, path, properties...)
    else
        return create_group(parent, path, properties...)
    end
end

"""
    rewrite_dataset(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, data; pv...)

  Deletes the dataset if it allready exists and creates and writes to a new dataset.
"""
function rewrite_dataset(parent::Union{HDF5.File,HDF5.Group}, name::AbstractString, data;
                         pv...)
    if haskey(parent, name)
        delete_object(parent[name])
    end
    write_dataset(parent, name, data; pv...)
end

# ---------------------------------------- Helpers -----------------------------------------

function Base.show(io::IO, m::MIME"text/plain", output::Output)
    print(io, "Output (store_hdf=", output.store_hdf,
          ", store_locally=", output.store_locally, ")")
    output.store_hdf ?
    print(io, ", simulation: ", HDF5.name(output.simulation), " (file: ",
          output.simulation.file.filename, "):") : print(":")

    if output.physical_transform !== identity
        print(io, "\nphysical_transform: ", nameof(output.physical_transform))
    end

    print(io, "\nDiagnostics:")
    for diagnostic in output.diagnostics
        print(io, "\n-")
        show(io, m, diagnostic)
    end
end

"""
    validate_resume_attributes(simulation::HDF5.Group, prob::SpectralODEProblem)

Checks if the existing HDF5 simulation group matches the current problem configuration
for critical parameters: Nx, Ny, dt, and real_transform.

!!! warning "User Responsibility"
    Advectra only checks a subset of critical numerical parameters (Nx, Ny, dt, real_transform). 
    It does **not** perform an exhaustive check of all physical input parameters.
    It is the responsibility of the user to ensure the simulations are physically compatible.
"""
function validate_resume_attributes(simulation::HDF5.Group, prob::SpectralODEProblem)
    @warn "Advectra only checks a subset of critical numerical parameters (Nx, Ny, dt, real_transform). " *
          "Advectra does not check if ALL input parameters are identical to the run " *
          "you are resuming from. It is the responsibility of the user to make sure " *
          "the simulations are compatible."

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

import Base.close
close(output::Output) = close(output.simulation.file)
