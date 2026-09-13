module AdvectraGPUArraysExt
# Inspired by: https://github.com/SciML/ComponentArrays.jl/blob/ff1c8c53a3bbb0dad84b902043ff4631ead74049/ext/ComponentArraysGPUArraysExt.jl

using Advectra: Advectra, State, _names, _unwrap, QuadraticTerm, fwd, bwd, pad!, unpad!
using GPUArrays: GPUArrays, AbstractGPUArray, generic_rmul!, @allowscalar

# ----------------------------------- General Interface ------------------------------------

# TODO move GPU code here...

# ------------------------------------- State Related --------------------------------------

const GPUState = State{Names,T,<:AbstractGPUArray} where {T,Names}
const GPUField = Field{T,N,D,<:AbstractGPUArray} where {T,N,D}

# (GPU-safe) mapreduce.
function Base.mapreduce(f, op, state::GPUState; kwargs...)
    mapreduce(f, op, get_data(state); kwargs...)
end
function Base.mapreduce(f, op, state::GPUState, args...; kwargs...)
    mapreduce(f, op, get_data(state), map(_unwrap, args)...; kwargs...)
end
function Base.mapreduce(f, op, state::GPUState,
                        args::Vararg{Union{Base.AbstractBroadcasted,AbstractArray}};
                        kwargs...,)
    mapreduce(f, op, get_data(state), map(_unwrap, args)...; kwargs...)
end

# Elementwise map
function Base.map(f, state::GPUState, args...)
    data = map(f, get_data(state), map(_unwrap, args)...)
    State(_names(state), data, getfield(state, :domain))
end
function Base.map(f, state::GPUState,
                  args::Vararg{Union{Base.AbstractBroadcasted,AbstractArray}})
    data = map(f, get_data(state), map(_unwrap, args)...)
    State(_names(state), data, getfield(state, :domain))
end

function Base.count(pred::Function, state::GPUState; dims=:, init=0)
    mapreduce(pred, Base.add_sum, get_data(state); init, dims)
end

# any/all
Base.any(state::GPUState{Bool}) = mapreduce(identity, |, get_data(state))
Base.all(state::GPUState{Bool}) = mapreduce(identity, &, get_data(state))
Base.any(f::Function, state::GPUState) = mapreduce(f, |, get_data(state))
Base.all(f::Function, state::GPUState) = mapreduce(f, &, get_data(state))

# -------------------------------- Linear Algebra Related ----------------------------------

import LinearAlgebra
LinearAlgebra.dot(x::GPUState, y::GPUState) = LinearAlgebra.dot(get_data(x), get_data(y))
LinearAlgebra.dot(x::GPUState, y::AbstractGPUArray) = LinearAlgebra.dot(get_data(x), y)
LinearAlgebra.dot(x::AbstractGPUArray, y::GPUState) = LinearAlgebra.dot(x, get_data(y))

LinearAlgebra.norm(state::GPUState, p::Real) = LinearAlgebra.norm(get_data(state), p)

# rmul! (`generic_rmul!` broadcasting internally (`x .*= b`))
LinearAlgebra.rmul!(state::GPUState, b::Number) = GPUArrays.generic_rmul!(state, b)

end