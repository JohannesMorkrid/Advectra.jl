"""
    State: a dense (size(domain)..., Nd) tensor of named Field's.
"""
struct State{Names,T,A<:AbstractArray{T,3}} <: AbstractArray{T,3}
    data::A
    domain::Domain
end

get_data(state::State) = getfield(state, :data)
get_domain(state::State) = getfield(state, :domain)
_names(::State{Names}) where {Names} = Names

@inline _getprop(state::State, ::Val{:data}) = get_data(state)
@inline _getprop(state::State, ::Val{:domain}) = get_domain(state)

@generated function _getprop(state::State{Names}, ::Val{name}) where {Names,name}
    idx = findfirst(==(name), Names)
    idx === nothing &&
        return :(throw(ArgumentError("State has no field :$($(QuoteNode(name))); valid fields: $Names")))
    return :(Field(view(get_data(state),:,:,($idx)), get_domain(state)))
    # TODO generalize the view
end

@inline Base.getproperty(state::State, name::Symbol) = _getprop(state, Val(name))

# Bare `S.n = val` is disallowed. Field is a view, so the only honest
# semantics would be an in-place copy, and doing that silently on `=` risks
# masking real bugs (e.g. typos). Force explicit `.=` instead.
function Base.setproperty!(::State, name::Symbol, val)
    throw(ArgumentError("State fields are views -- use `S.$name .= value` for in-place mutation, not `=`"))
end

# ------------------------------------ Array Interface -------------------------------------

Base.size(state::State) = size(get_data(state))
Base.axes(state::State) = axes(get_data(state))
Base.IndexStyle(::Type{<:State{Names,T,A}}) where {Names,T,A} = IndexStyle(A)

@inline Base.getindex(state::State, I...) = getindex(get_data(state), I...)
@inline Base.setindex!(state::State, v, I...) = setindex!(get_data(state), v, I...)
Base.fill!(state::State, x) = (fill!(get_data(state), x); state)

Base.@propagate_inbounds Base.view(state::State, I...) = view(get_data(state), I...)
Base.print_array(io::IO, state::State) = Base.print_array(io, state.data)

# ------------------------------------- Similar/copy ---------------------------------------

Base.similar(s::State) = State(_names(s), similar(get_data(s)), get_domain(s))
function Base.similar(state::State, ::Type{T}) where {T}
    State(_names(state), similar(get_data(state), T), get_domain(state))
end

function Base.similar(state::State{Names}, ::Type{T}, dims::Dims) where {Names,T}
    if length(dims) == 3 && dims[3] == length(Names)
        data = similar(get_data(state), T, dims)
        State{Names,T,typeof(data)}(data, get_domain(state))
    else
        # Requested shape no longer has Nd fields matching Names -- a State
        # would be a lie here, so a plain array is the honest result.
        similar(get_data(state), T, dims)
    end
end

function Base.copyto!(dest::State{Names}, src::State{Names}) where {Names}
    (copyto!(get_data(dest), get_data(src)); dest)
end
Base.copyto!(dest::State, src::AbstractArray) = (copyto!(get_data(dest), src); dest)
Base.copyto!(dest::AbstractArray, src::State) = copyto!(dest, get_data(src))

function Base.deepcopy_internal(state::State, stackdict::IdDict)
    haskey(stackdict, state) && return stackdict[state]
    data′ = Base.deepcopy_internal(get_data(state), stackdict)
    state′ = State(_names(state), data′, get_domain(state))
    stackdict[state] = state′
    state′
end

# -------------------------------- Broadcasting Machinery ----------------------------------

struct StateStyle{S<:Broadcast.BroadcastStyle} <: Broadcast.AbstractArrayStyle{3} end

StateStyle{S}(::Val) where {S} = StateStyle{S}()

function Base.BroadcastStyle(::Type{<:State{Names,T,A}}) where {Names,T,A}
    StateStyle{typeof(Broadcast.BroadcastStyle(A))}()
end

function Base.BroadcastStyle(::StateStyle{S1}, ::StateStyle{S2}) where {S1,S2}
    StateStyle{typeof(Broadcast.result_style(S1(), S2()))}()
end

@inline _unwrap(state::State) = get_data(state)
@inline function _unwrap(bc::Broadcast.Broadcasted{StateStyle{S}}) where {S}
    Broadcast.Broadcasted{S}(bc.f, map(_unwrap, bc.args))
end

# Walk a (possibly nested) Broadcasted tree and verify every State operand's
# Names matches `expected`. Runs once per broadcast
# call, not per element, so cost is O(number of operands).
function _check_consistency(args::Tuple, names, domain)
    (_check_consistency(first(args), names, domain);
     _check_consistency(Base.tail(args), names, domain))
end
_check_consistency(::Tuple{}, names, domain) = nothing
# essential
function _check_consistency(bc::Broadcast.Broadcasted, names, domain)
    _check_consistency(bc.args, names, domain)
end
function _check_consistency(x::State, names, domain)
    _names(x) === names ||
        throw(ArgumentError("Mismatched State field names in broadcast: expected $names, found $(_names(x))"))
    get_domain(x) == domain || throw(ArgumentError("Mismatched State domains in broadcast"))
    nothing
end
# Catch all
_check_consistency(x, names, domain) = nothing

# In-place: `S3 .= S1 .+ S2` -- no new State allocated.
function Base.copyto!(dest::State{Names},
                      bc::Broadcast.Broadcasted{StateStyle{S}}) where {Names,S}
    _check_consistency(bc, Names, get_domain(dest))
    copyto!(get_data(dest), _unwrap(bc))
    dest
end

_find_state(args::Tuple) = _find_state(first(args), Base.tail(args))
_find_state(x::State, rest) = x
function _find_state(bc::Broadcast.Broadcasted, rest)
    found = _find_state(bc.args)
    found === nothing ? _find_state(rest) : found
end
_find_state(x, rest) = _find_state(rest)
_find_state(::Tuple{}) = nothing

function Base.similar(bc::Broadcast.Broadcasted{StateStyle{S}}, ::Type{T}) where {S,T}
    st = _find_state(bc.args)
    _check_consistency(bc, _names(st), get_domain(st))
    data = similar(get_data(st), T)
    State{_names(st),T,typeof(data)}(data, get_domain(st))
end

Base.dataids(state::State) = Base.dataids(get_data(state))
function Broadcast.broadcast_unalias(dest::State, src)
    src isa State && get_data(dest) === get_data(src) && return src
    return Base.unalias(dest, src)
end

# ------------------------------------- Constructors ---------------------------------------

State(state::State) = state
State{T}(state::State) where {T} = T.(state)

"""
State(names::NTuple{Nd,Symbol}, domain::Domain, representation::Representation=Physical())

    Blank allocator, Representation-aware, rep defaults to Physical().
"""
function State(names::NTuple{Nd,Symbol}, domain::Domain,
               representation::Representation=Physical()) where {Nd}
    allunique(names) || throw(ArgumentError("Field names must be unique, got $names"))
    data = allocate(domain, Nd, representation)
    State{names,eltype(data),typeof(data)}(data, domain)
end

State(nt::NamedTuple, domain::Domain) = State(domain; nt...)

"""
State(args...)

    Attempt at making Vararg constructor such as State(:n, :p, :T, domain)
"""
function State(args...)
    if last(args) isa Representation
        rep, domain = args[end], args[end-1]
        domain isa Domain || throw(ArgumentError("expected (names..., domain, rep)"))
        names = args[1:(end-2)]
    else
        rep, domain = Physical(), args[end]
        domain isa Domain || throw(ArgumentError("expected (names..., domain[, rep])"))
        names = args[1:(end-1)]
    end
    all(name -> name isa Symbol, names) ||
        throw(ArgumentError("expected field names as Symbols, got $(names) -- did you mean State(nt::NamedTuple, domain)?"))
    State(Tuple(names), domain, rep)
end

"""
State(data::AbstractArray{T,3}, domain::Domain)

    Auto-named fields :f1, :f2, ...
"""
function State(data::AbstractArray{T,3}, domain::Domain) where {T}
    physsz, specsz = size(domain), spectral_size(domain)
    physet, specet = physical_eltype(domain), spectral_eltype(domain)
    leading = size(data)[1:2]

    is_physical = leading == physsz && T == physet
    is_spectral = leading == specsz && T == specet
    (is_physical || is_spectral) ||
        throw(DimensionMismatch("data (size=$(size(data)), eltype=$T) matches neither " *
                                "Physical ($physsz, $physet) nor Spectral ($specsz, $specet)"))

    representation = is_physical ? Physical() : Spectral()
    target = memory_type(domain, representation)
    data = target(data)

    names = ntuple(i -> Symbol(:f, i), size(data, 3))
    State{names,eltype(data),typeof(data)}(data, domain)
end

"""
State(names::NTuple{Nd,Symbol}, data::AbstractArray{T,3}, domain::Domain)

    Fast/programmatic constructor, no validation.
"""
function State(names::NTuple{Nd,Symbol}, data::AbstractArray{T,3},
               domain::Domain) where {Nd,T}
    State{names,T,typeof(data)}(data, domain)
end

# Different way of creating a Field in a State.
_materialize(f::Function, domain) = initial_condition(f, domain)
function _materialize((f, kw)::Tuple{<:Function,<:NamedTuple}, domain)
    initial_condition(f, domain; kw...)
end
_materialize(f::Field, domain) = get_data(f)
_materialize(a::AbstractArray, domain) = a

"""
    Mixed data / Field / function constructor, representation inferred.
"""
function State(domain::Domain; fields...)
    names = keys(fields)
    processed = map(v -> _materialize(v, domain), values(fields))
    shapes = map(size, processed)
    eltypes = map(eltype, processed)

    physsz, specsz = size(domain), spectral_size(domain)
    physet, specet = physical_eltype(domain), spectral_eltype(domain)

    is_physical = all(==(physsz), shapes) && all(==(physet), eltypes)
    is_spectral = all(==(specsz), shapes) && all(==(specet), eltypes)

    (is_physical || is_spectral) ||
        throw(DimensionMismatch("fields (shapes=$shapes, eltypes=$eltypes) match neither " *
                                "Physical ($physsz, $physet) nor Spectral ($specsz, $specet)"))

    representation = is_physical ? Physical() : Spectral()
    target = memory_type(domain, representation)
    processed = map(target, processed)
    data = stack(processed)
    State{names,eltype(data),typeof(data)}(data, domain)
end

# ------------------------------------ User Interface --------------------------------------

"""
unpack_state(S) -> (Field(view 1), Field(view 2), ..., Field(view Nd))

@generated so this expands to a flat tuple literal, no loop, no allocation
# beyond the (stack-allocatable) Field wrappers themselves.
"""
@generated function unpack_state(state::State{Names}) where {Names}
    n = length(Names)
    exprs = [:(Field(view(data,:,:,($i)), domain)) for i in 1:n]
    quote
        data = get_data(state)
        domain = get_domain(state)
        tuple($(exprs...))
    end
end
unpack_state(A) = eachslice(A; dims=ndims(A)) # Fallback

Base.propertynames(::State{Names}) where {Names} = (Names..., :data, :domain)
Base.keys(state::State) = _names(state)
Base.haskey(state::State, name::Symbol) = name in _names(state)

function map_field!(f!, dS::State{Names}, S::State{Names}) where {Names}
    ddata, sdata = get_data(dS), get_data(S)
    for i in 1:length(Names)
        f!(view(ddata,:,:,i), view(sdata,:,:,i))
    end
    dS
end

function map_field(f, S::State)
    dS = similar(S)
    map_field!(dS, S) do dest, src
        dest .= f(src)
    end
    dS
end

function spectral_transform!(dS::State{Names}, p::P,
                             S::State{Names}) where {Names,P<:FFTW.Plan}
    map_field!(dS, S) do dest, src
        LinearAlgebra.mul!(dest, p, src)
    end
end

function Base.:(==)(s1::State{N1}, s2::State{N2}) where {N1,N2}
    N1 == N2 && get_domain(s1) == get_domain(s2) && get_data(s1) == get_data(s2)
end

function Base.isapprox(s1::State{N1}, s2::State{N2}; kwargs...) where {N1,N2}
    N1 == N2 && get_domain(s1) == get_domain(s2) &&
        isapprox(get_data(s1), get_data(s2); kwargs...)
end

# ------------------------------------ BLAS Interface --------------------------------------

Base.parent(state::State) = get_data(state)

# Low-level array interface 
Base.pointer(state::State) = pointer(get_data(state))
Base.unsafe_convert(::Type{P}, s::State) where {P} = Base.unsafe_convert(P, get_data(s))
Base.strides(state::State) = strides(get_data(state))
Base.stride(state::State, k::Integer) = stride(get_data(state), k)

# Level-1 BLAS (axpy!/axpby!)

function LinearAlgebra.axpy!(alpha::Number, x::State{Names}, y::State{Names}) where {Names}
    LinearAlgebra.axpy!(alpha, get_data(x), get_data(y))
    y
end

function LinearAlgebra.axpby!(alpha::Number, x::State{Names}, beta::Number,
                              y::State{Names}) where {Names}
    LinearAlgebra.axpby!(alpha, get_data(x), beta, get_data(y))
    y
end

# ---------------------------------- Adapt Compatibility -----------------------------------

function Adapt.adapt_structure(to, state::State)
    State(_names(state), Adapt.adapt(to, get_data(state)), get_domain(state))
end

Adapt.parent_type(::Type{State{Names,T,A}}) where {Names,T,A} = A

function Adapt.adapt_storage(::Type{State{Names,T,A}}, xs::AbstractArray) where {Names,T,A}
    Adapt.adapt_storage(A, xs)
end

# --------------------------------- AdvectraGPUArraysExt -----------------------------------

const GPUState = State{Names,T,<:AbstractGPUArray} where {T,Names}

function Base.copyto!(dest::State, bc::Broadcasted{<:GPUArrays.AbstractGPUArrayStyle})
    (copyto!(get_data(dest), _unwrap(bc)); dest)
end

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
    State(_names(state), data, get_domain(state))
end
function Base.map(f, state::GPUState,
                  args::Vararg{Union{Base.AbstractBroadcasted,AbstractArray}})
    data = map(f, get_data(state), map(_unwrap, args)...)
    State(_names(state), data, get_domain(state))
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

LinearAlgebra.norm(state::GPUState) = LinearAlgebra.norm(get_data(state))
LinearAlgebra.norm(state::GPUState, p::Real) = LinearAlgebra.norm(get_data(state), p)

# rmul! (`generic_rmul!` broadcasting internally (`x .*= b`))
LinearAlgebra.rmul!(s::GPUState, b::Number) = (LinearAlgebra.rmul!(get_data(s), b); s)