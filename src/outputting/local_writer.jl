# ------------------------------------------------------------------------------------------
#                                       Local Writer
# ------------------------------------------------------------------------------------------
#
# NOTE: LocalWriter is scaffolding for future in-memory storage. It is not yet wired into
# the Output constructor. See TODO comments throughout.

"""
    LocalWriter <: AbstractWriter

  Stores diagnostic samples in memory as Julia arrays. Useful for interactive sessions
  where the user wants to inspect results without writing to disk.

  Each diagnostic gets its own entry in `data` (keyed by `diagnostic.name`) holding a
  vector of samples, and a corresponding time vector in `t`.
"""
mutable struct LocalWriter <: AbstractWriter
    data::Dict{String,Vector}  # diagnostic.name => vector of samples
    t::Dict{String,Vector{Float64}}
end

LocalWriter() = LocalWriter(Dict{String,Vector}(), Dict{String,Vector{Float64}}())

# ------------------------------------ Interface ------------------------------------------

function initialize!(writer::LocalWriter, prob, diagnostics, initial_samples, strides, t0)
    N_steps = compute_number_of_steps(prob)
    for (diagnostic, sample, stride) in zip(diagnostics, initial_samples, strides)
        N_samples = cld(N_steps, stride) + 1
        # TODO allocate based on sample type (scalar vs array)
        writer.data[diagnostic.name] = Vector{typeof(sample)}(undef, N_samples)
        writer.t[diagnostic.name] = zeros(Float64, N_samples)

        # Store initial sample
        writer.data[diagnostic.name][1] = isa(sample, AbstractArray) ? copy(sample) : sample
        writer.t[diagnostic.name][1] = t0 #first(prob.tspan)
    end
end

function write!(writer::LocalWriter, diagnostic, idx::Integer, sample, time)
    writer.data[diagnostic.name][idx] = isa(sample, AbstractArray) ? copy(sample) : sample
    writer.t[diagnostic.name][idx] = time
end

# flush! and close fall back to the no-op defaults in abstract_writer.jl

# TODO implement checkpoints

# ----------------------------------- Base.show -------------------------------------------

function Base.show(io::IO, ::MIME"text/plain", writer::LocalWriter)
    keys_str = isempty(writer.data) ? "(empty)" : join(keys(writer.data), ", ")
    print(io, "LocalWriter (diagnostics: ", keys_str, ")")
end
