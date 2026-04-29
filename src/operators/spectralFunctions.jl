# ------------------------------------------------------------------------------------------
#                                    Spectral Functions                                     
# ------------------------------------------------------------------------------------------

# ------------------------------------- Construction ---------------------------------------

struct SpectralFunction{F<:Function} <: SpectralOperator
    f::F
    quadratic_term::QuadraticTerm
    function SpectralFunction(f::Function, quadratic_term::QuadraticTerm)
        new{typeof(f)}(f, quadratic_term)
    end
end

# ------------------------------------- Constructors ---------------------------------------

function operator_dependencies(::Union{Val{:reciprocal},Val{:spectral_exp},
                                       Val{:spectral_expm1},
                                       Val{:spectral_log}}, ::Type)
    [OperatorRecipe(:quadratic_term)]
end

function build_operator(::Val{:reciprocal}, domain::AbstractDomain; quadratic_term,
                        kwargs...)
    SpectralFunction(u -> div(1, u), quadratic_term)
end

function build_operator(::Val{:spectral_exp}, domain::AbstractDomain; quadratic_term,
                        kwargs...)
    SpectralFunction(exp, quadratic_term)
end

function build_operator(::Val{:spectral_expm1}, domain::AbstractDomain; quadratic_term,
                        kwargs...)
    SpectralFunction(expm1, quadratic_term)
end

function build_operator(::Val{:spectral_log}, domain::AbstractDomain; quadratic_term,
                        kwargs...)
    SpectralFunction(log, quadratic_term)
end

# ------------------------------------- Main Methods ---------------------------------------

# In-place method
function (op::SpectralFunction)(du::AbstractArray, u::AbstractArray, args...; kwargs...)
    @unpack U, V, up, padded, transforms, dealiasing_coefficient = op.quadratic_term
    mul!(U, bwd(transforms), padded ? pad!(up, u, typeof(transforms)) : u)
    @. V = op.f(dealiasing_coefficient * U, args...; kwargs...)
    mul!(padded ? up : du, fwd(transforms), V)
    padded ? du .= unpad!(du, up, typeof(transforms)) ./ dealiasing_coefficient : up
end

# Out-of-place
function (op::SpectralFunction)(u::AbstractArray, args...; kwargs...)
    op(similar(u), u, args...; kwargs...)
end