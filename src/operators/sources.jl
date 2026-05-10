# ------------------------------------------------------------------------------------------
#                                     SpectralConstant                                      
# ------------------------------------------------------------------------------------------

struct SpectralConstant{T<:Number} <: SpectralOperator
    value::T
    SpectralConstant(; val) = new{typeof(val)}(val)
    SpectralConstant(domain::Domain; val) = new{typeof(val)}(val)
end

function build_operator(::Val{:spectral_constant}, domain::AbstractDomain; val=1, kwargs...)
    SpectralConstant(domain; val)
end

# ----------------------- SpectralConstant-SpectralConstant Methods ------------------------

function Base.:+(a::SpectralConstant, b::SpectralConstant)
    SpectralConstant(; val=a.value + b.value)
end
function Base.:-(a::SpectralConstant, b::SpectralConstant)
    SpectralConstant(; val=a.value - b.value)
end
function Base.:*(a::SpectralConstant, b::SpectralConstant)
    SpectralConstant(; val=a.value * b.value)
end
function Base.:/(a::SpectralConstant, b::SpectralConstant)
    SpectralConstant(; val=a.value / b.value)
end

# ---------------------------- SpectralConstant-Number Methods -----------------------------

Base.:*(sc::SpectralConstant, c::Number) = SpectralConstant(; val=sc.value * c)
Base.:*(c::Number, sc::SpectralConstant) = SpectralConstant(; val=c * sc.value)
Base.:/(sc::SpectralConstant, c::Number) = SpectralConstant(; val=sc.value * c)
Base.:/(c::Number, sc::SpectralConstant) = SpectralConstant(; val=c * sc.value)
Base.:-(sc::SpectralConstant) = SpectralConstant(; val=-sc.value)

# ---------------------------- SpectralConstant-Array Methods ------------------------------

function Base.:+(field::AbstractArray, sc::SpectralConstant)
    out = copy(field)
    @allowscalar out[1] += sc.value
    return out
end
Base.:+(sc::SpectralConstant, field::AbstractArray) = Base.:+(field, sc)

function Base.:-(field::AbstractArray, sc::SpectralConstant)
    out = copy(field)
    @allowscalar out[1] -= sc.value
    return out
end
Base.:-(sc::SpectralConstant, field::AbstractArray) = Base.:+(field, sc)

# ------------------------------------------------------------------------------------------
#                                          Sources                                          
# ------------------------------------------------------------------------------------------

struct Source <: SpectralOperator
    shape::AbstractArray
end

# ------------------------------------- Main Methods ---------------------------------------