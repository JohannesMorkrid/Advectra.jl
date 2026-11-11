# ------------------------------------------------------------------------------------------
#                                        LocalWriter                                        
# ------------------------------------------------------------------------------------------
#
# NOTE: LocalWriter is scaffolding for future in-memory storage. The interface is complete
# but the implementation is intentionally minimal. See TODO comments throughout.

"""
    LocalWriter <: AbstractWriter

  Stores diagnostic samples in memory as Julia arrays. Useful for interactive sessions
  where the user wants to inspect results without writing to disk.

  Each diagnostic gets its own entry in `data` (keyed by `diagnostic.name`) holding a
  vector of samples, and a corresponding time vector in `t`.
"""
struct LocalWriter <: AbstractWriter
    data::Dict{String,Vector}
    t::Dict{String,Vector{Float64}}
    simulation_name::Ref{String}   # set by initialize!, for display purposes only
end

# By default empty dictionaries before `initialized!
LocalWriter() = LocalWriter(Dict{String,Vector}(), Dict{String,Vector{Float64}}(), Ref("advectra"))

# --------------------------------------- Interface ----------------------------------------

function initialize!(writer::LocalWriter, prob, diagnostics, initial_samples, strides, t0;
    simulation_name::String="advectra", kwargs...)
    writer.simulation_name[] = simulation_name
    N_steps = compute_number_of_steps(prob)
    for (diagnostic, sample) in zip(diagnostics, initial_samples)
        N_samples = N_steps #cld(N_steps, stride) + 1 # TODO integrate schedules here
        # TODO: allocate based on sample type (scalar vs array)
        writer.data[diagnostic.name] = Vector{typeof(sample)}(undef, N_samples)
        writer.t[diagnostic.name] = zeros(Float64, N_samples)
        # Store initial sample
        writer.data[diagnostic.name][1] = isa(sample, AbstractArray) ? copy(sample) : sample
        writer.t[diagnostic.name][1] = t0
    end
end

function write_to!(writer::LocalWriter, diagnostic, idx::Integer, sample, time)
    isinitialized(writer) ||
        error("LocalWriter has not been initialized — call initialize! before write_to!")
    writer.data[diagnostic.name][idx] = isa(sample, AbstractArray) ? copy(sample) : sample
    writer.t[diagnostic.name][idx] = time
end

# flush!, close, store_checkpoint! fall back to the no-op defaults in abstract_writer.jl

isinitialized(writer::LocalWriter) = !isempty(writer.data)
supports_checkpoint(::LocalWriter) = false

# --------------------------------------- Base.show ----------------------------------------

function Base.show(io::IO, ::MIME"text/plain", writer::LocalWriter)
    keys_str = isempty(writer.data) ? "(empty)" : join(keys(writer.data), ", ")
    print(io, "LocalWriter (simulation: ", writer.simulation_name[],
        ", diagnostics: ", keys_str, ")")
end