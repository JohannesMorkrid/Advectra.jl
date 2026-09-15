# Fundamental building block for Advectra code
struct Field{T,N,D<:AbstractDomain,A<:AbstractArray{T,N}} <: AbstractArray{T,N}
    data::A
    domain::D
end

get_data(field::Field) = getfield(field, :data)
get_domain(field::Field) = getfield(field, :domain)

# ------------------------------------ Array Interface -------------------------------------

Base.size(field::Field) = size(get_data(field))
Base.parent(field::Field) = get_data(field)
Base.IndexStyle(::Type{<:Field{T,N,D,A}}) where {T,N,D,A} = IndexStyle(A)

Base.@propagate_inbounds Base.getindex(field::Field, I...) = getindex(get_data(field), I...)
Base.@propagate_inbounds Base.setindex!(field::Field, v, I...) = setindex!(get_data(field), v, I...)
Base.fill!(field::Field, x) = (fill!(get_data(field), x); field)

Base.@propagate_inbounds Base.view(field::Field, I...) = view(get_data(field), I...)
# To be compatible with GPU Arrays
Base.print_array(io::IO, field::Field) = Base.print_array(io, get_data(field))

function Base.showarg(io::IO, field::Field, toplevel)
    print(io, "Field(", typeof(get_data(field)), ", ", typeof(get_domain(field)), ")")
end

# ------------------------------------- Similar/copy ---------------------------------------

function Base.similar(field::Field, ::Type{T}, dims::Dims) where {T}
    Field(similar(get_data(field), T, dims), get_domain(field))
end

Base.copyto!(dest::Field, src::Field) = (copyto!(get_data(dest), get_data(src)); dest)
Base.copyto!(dest::Field, src::AbstractArray) = (copyto!(get_data(dest), src); dest)
Base.copyto!(dest::AbstractArray, src::Field) = (copyto!(dest, get_data(src)))

function Base.deepcopy_internal(field::Field, stackdict::IdDict)
    haskey(stackdict, field) && return stackdict[field]
    data′ = Base.deepcopy_internal(get_data(field), stackdict)
    field′ = Field(data′, get_domain(field))
    stackdict[field] = field′
    field′
end

# -------------------------------- Broadcasting Machinery ----------------------------------

import Base.Broadcast: BroadcastStyle, Broadcasted

# Construction related
struct FieldStyle{S<:BroadcastStyle} <: BroadcastStyle end
FieldStyle(s::BroadcastStyle) = FieldStyle{typeof(s)}()
BroadcastStyle(::Type{<:Field{T,N,D,A}}) where {T,N,D,A} = FieldStyle(BroadcastStyle(A))

# Combination rules: two Fields combine via their inner styles
function BroadcastStyle(a::FieldStyle{S1}, ::FieldStyle{S2}) where {S1,S2}
    FieldStyle(Broadcast.result_style(S1(), S2()))
end
function BroadcastStyle(::FieldStyle{S}, other::BroadcastStyle) where {S}
    FieldStyle(Broadcast.result_style(S(), other))
end

# Taken from manual: https://docs.julialang.org/en/v1/manual/interfaces/#Broadcast-Styles
find_field(bc::Broadcasted) = find_field(bc.args)
find_field(args::Tuple) = find_field(find_field(args[1]), Base.tail(args))
find_field(x) = x
find_field(::Tuple{}) = nothing
find_field(f::Field, rest) = f
find_field(::Any, rest) = find_field(rest)

# Strip Field wrappers out of the expression tree to get inner style
@inline _unwrap(x) = x # Generic catch all
@inline _unwrap(x::Field) = parent(x)
@inline _unwrap(bc::Broadcasted{FieldStyle{S}}) where {S} = Broadcasted{S}(bc.f,
    map(_unwrap,
        bc.args),
    bc.axes)

# Walk a (possibly nested) Broadcasted tree and verify every Field operand's domain matches.
_check_domain(args::Tuple, domain) =
    (_check_domain(first(args), domain); _check_domain(Base.tail(args), domain))
_check_domain(::Tuple{}, domain) = nothing
_check_domain(bc::Broadcasted, domain) = _check_domain(bc.args, domain)
_check_domain(x::Field, domain) = get_domain(x) == domain ||
    throw(ArgumentError("Mismatched Field domains in broadcast"))
_check_domain(x, domain) = nothing

# Allocation of output Field
function Base.similar(bc::Broadcasted{FieldStyle{S}}, ::Type{ElType}) where {S,ElType}
    field = find_field(bc) # Get field to preserve Domain
    _check_domain(bc, get_domain(field))
    Field(similar(_unwrap(bc), ElType), get_domain(field))
end

# Usefull for GPU kernel-fusion
@inline function Base.copyto!(dest::Field, bc::Broadcasted{FieldStyle{S}}) where {S}
    _check_domain(bc, get_domain(dest))
    copyto!(get_data(dest), _unwrap(bc))
    dest
end

@inline function Base.copyto!(dest::AbstractArray, bc::Broadcasted{FieldStyle{S}}) where {S}
    copyto!(dest, _unwrap(bc))
    dest
end

Base.dataids(field::Field) = Base.dataids(get_data(field))
function Broadcast.broadcast_unalias(dest::Field, src)
    src isa Field && get_data(dest) === get_data(src) && return src
    return Base.unalias(dest, src)
end

# ------------------------------------ User Interface --------------------------------------

function Base.:(==)(f1::Field, f2::Field)
    get_domain(f1) == get_domain(f2) && get_data(f1) == get_data(f2)
end

function Base.isapprox(f1::Field, f2::Field; kwargs...)
    get_domain(f1) == get_domain(f2) && isapprox(get_data(f1), get_data(f2); kwargs...)
end

# ------------------------------------ BLAS Interface --------------------------------------

# Low-level array interface 
Base.pointer(field::Field) = pointer(get_data(field))
Base.unsafe_convert(::Type{P}, f::Field) where {P} = Base.unsafe_convert(P, get_data(f))
Base.strides(field::Field) = strides(get_data(field))
Base.stride(field::Field, k::Integer) = stride(get_data(field), k)

_rewrap(result::AbstractArray, template::Field) = Field(result, get_domain(template))
_rewrap(result, template) = result # scalars, etc. pass through untouched

LinearAlgebra.dot(x::Field, y::Field) = LinearAlgebra.dot(get_data(x), get_data(y))
LinearAlgebra.dot(x::Field, y::AbstractArray) = LinearAlgebra.dot(get_data(x), y)
LinearAlgebra.dot(x::AbstractArray, y::Field) = LinearAlgebra.dot(x, get_data(y))

LinearAlgebra.axpy!(alpha::Number, x::Field, y::Field) = (LinearAlgebra.axpy!(alpha, get_data(x), get_data(y)); y)
LinearAlgebra.axpy!(alpha::Number, x::Field, y::AbstractArray) = (LinearAlgebra.axpy!(alpha, get_data(x), y); y)
LinearAlgebra.axpy!(alpha::Number, x::AbstractArray, y::Field) = (LinearAlgebra.axpy!(alpha, x, get_data(y)); y)

LinearAlgebra.axpby!(alpha::Number, x::Field, beta::Number, y::Field) = (LinearAlgebra.axpby!(alpha, get_data(x), beta, get_data(y)); y)
LinearAlgebra.axpby!(alpha::Number, x::Field, beta::Number, y::AbstractArray) = (LinearAlgebra.axpby!(alpha, get_data(x), beta, y); y)
LinearAlgebra.axpby!(alpha::Number, x::AbstractArray, beta::Number, y::Field) = (LinearAlgebra.axpby!(alpha, x, beta, get_data(y)); y)

LinearAlgebra.rmul!(x::Field, b::Number) = (LinearAlgebra.rmul!(get_data(x), b); x)

# TODO add missing BLAS interface inspired by ComponentArrays.jl and FELTOR
#using LinearAlgebra
for f in (:mul!, :ldiv!, :rdiv!, :norm, :tr, :det, :inv, :cross, :qr, :lu, :cholesky)
    @eval function LinearAlgebra.$f(A::Field, args...; kwargs...)
        result = LinearAlgebra.$f(parent(A), map(_unwrap, args)...; kwargs...)
        _rewrap(result, A)
    end
end

# ---------------------------------- Adapt Compatibility -----------------------------------

function Adapt.adapt_structure(to, field::Field)
    Field(Adapt.adapt(to, get_data(field)), get_domain(field))
end

Adapt.parent_type(::Type{Field{T,N,D,A}}) where {T,N,D,A} = A

function Adapt.adapt_storage(::Type{Field{T,N,D,A}}, xs::AbstractArray) where {T,N,D,A}
    Adapt.adapt_storage(A, xs)
end

# --------------------------------- AdvectraGPUArraysExt -----------------------------------
const GPUField = Field{T,N,D,<:AbstractGPUArray} where {T,N,D}

function Base.copyto!(dest::Field, bc::Broadcasted{<:GPUArrays.AbstractGPUArrayStyle})
    _check_domain(bc, get_domain(dest))
    copyto!(get_data(dest), _unwrap(bc))
    dest
end
# # (GPU-safe) mapreduce.
function Base.mapreduce(f, op, field::GPUField, args...; kwargs...)
    mapreduce(f, op, get_data(field), map(_unwrap, args)...; kwargs...)
end
function Base.mapreduce(f, op, field::GPUField,
    args::Vararg{Union{Base.AbstractBroadcasted,AbstractArray}};
    kwargs...,)
    mapreduce(f, op, get_data(field), map(_unwrap, args)...; kwargs...)
end

# Elementwise map
function Base.map(f, field::GPUField, args...)
    data = map(f, get_data(field), map(_unwrap, args)...)
    Field(data, get_domain(field))
end
function Base.map(f, field::GPUField,
    args::Vararg{Union{Base.AbstractBroadcasted,AbstractArray}})
    data = map(f, get_data(field), map(_unwrap, args)...)
    Field(data, get_domain(field))
end

function Base.count(pred::Function, field::GPUField; dims=:, init=0)
    mapreduce(pred, Base.add_sum, get_data(field); init, dims)
end

# any/all
Base.any(field::GPUField{Bool}) = mapreduce(identity, |, get_data(field))
Base.all(field::GPUField{Bool}) = mapreduce(identity, &, get_data(field))
Base.any(f::Function, field::GPUField) = mapreduce(f, |, get_data(field))
Base.all(f::Function, field::GPUField) = mapreduce(f, &, get_data(field))