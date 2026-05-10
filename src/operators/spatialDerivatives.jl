# ------------------------------------------------------------------------------------------
#                                    Spatial Derivatives                                    
# ------------------------------------------------------------------------------------------

# ------------------------------------ Fourier Domain --------------------------------------

# -------------------------------- Construction Interface ----------------------------------

function build_operator(::Val{:diff_x}, domain::Domain; order=1, kwargs...)
    ElwiseOperator(im .* transpose(domain.kx); order=order)
end

function build_operator(::Val{:diff_xx}, domain::Domain; order=1, kwargs...)
    build_operator(Val(:diff_x), domain; order=2 * order)
end

function build_operator(::Val{:diff_y}, domain::Domain; order=1, kwargs...)
    ElwiseOperator(im .* domain.ky; order=order)
end

function build_operator(::Val{:diff_yy}, domain::Domain; order=1, kwargs...)
    build_operator(Val(:diff_y), domain; order=2 * order)
end

# Helper function 
function get_laplacian(domain::Domain)
    ks = wave_vectors(domain)
    N = length(ks)
    reshaped = ntuple(i -> reshape(ks[i], ntuple(j -> j == i ? length(ks[i]) : 1, N)), N)
    return -mapreduce(k -> k .^ 2, (a, b) -> a .+ b, reshaped)
end

function build_operator(::Val{:laplacian}, domain::Domain; order=1, kwargs...)
    ElwiseOperator(get_laplacian(domain); order=order)
end

# ----------------------------------- Chebyshev Domain -------------------------------------

# Chebyshev
# function chebyshev_differentiation_matrix(x::AbstractVector)
#     Nx = length(x)
#     D = zeros(Nx, Nx)
#     for i in eachindex(x)
#         for j in eachindex(x)
#             if i == j
#                 D[j, i] = -x[j] / (2(1 - x[j]^2))
#             else
#                 ci = i != 1 && i != Nx ? 1 : 2
#                 cj = j != 1 && j != Nx ? 1 : 2
#                 D[i, j] = (ci / cj) * (-1)^(i + j) / (x[i] - x[j])
#             end
#         end
#     end
#     D[1, 1] = (2 * (Nx - 1)^2 + 1) / 6
#     D[Nx, Nx] = -(2 * (Nx - 1)^2 + 1) / 6

#     return D
# end

# ---------------------------------- Compound Operators ------------------------------------

# Quite similar to the poisson bracket method
struct GradDotGradOperator{T<:AbstractArray} <: NonLinearOperator
    diff_x::LinearOperator
    diff_y::LinearOperator
    quadratic_term::QuadraticTerm
    left::T
    right::T
    tmp::T
    function GradDotGradOperator(domain::AbstractDomain, diff_x::LinearOperator,
                                 diff_y::LinearOperator, quadratic_term::QuadraticTerm)
        tmp = zeros(spectral_size(domain)) |> domain.MemoryType{complex(domain.precision)}
        left = zero(tmp)
        right = zero(tmp)
        new{typeof(tmp)}(diff_x, diff_y, quadratic_term, left, right, tmp)
    end
end

function operator_dependencies(::Val{:grad_dot_grad}, ::Type{_}) where {_}
    [OperatorRecipe(:diff_x), OperatorRecipe(:diff_y), OperatorRecipe(:quadratic_term)]
end

function build_operator(::Val{:grad_dot_grad}, domain::Domain; diff_x, diff_y,
                        quadratic_term, kwargs...)
    GradDotGradOperator(domain, diff_x, diff_y, quadratic_term)
end

function (op::GradDotGradOperator)(out::AbstractArray, u::AbstractArray, v::AbstractArray)
    @unpack left, right, tmp, diff_x, diff_y, quadratic_term = op

    diff_x(left, u)
    diff_x(right, v)
    quadratic_term(out, left, right)
    diff_y(left, u)
    diff_y(right, v)
    quadratic_term(tmp, left, right)
    out .+= tmp
end

(op::GradDotGradOperator)(u::AbstractArray, v::AbstractArray) = op(similar(u), u, v)