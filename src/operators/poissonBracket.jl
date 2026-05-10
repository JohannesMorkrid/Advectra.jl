# ------------------------------------------------------------------------------------------
#                                      Poisson Bracket                                      
# ------------------------------------------------------------------------------------------

# ------------------------------------- Construction ---------------------------------------

struct PoissonBracket{T<:AbstractArray} <: NonLinearOperator
    diff_x::LinearOperator
    diff_y::LinearOperator
    quadratic_term::QuadraticTerm
    tmp::T
    qt_left::T
    qt_right::T

    function PoissonBracket(domain::AbstractDomain, diff_x::LinearOperator,
                            diff_y::LinearOperator, quadratic_term::QuadraticTerm)

        # Allocate
        tmp = zeros(spectral_size(domain)) |> domain.MemoryType{complex(domain.precision)}
        qt_left = zero(tmp)
        qt_right = zero(qt_left)

        new{typeof(tmp)}(diff_x, diff_y, quadratic_term, tmp, qt_left, qt_right)
    end
end

function operator_dependencies(::Val{:poisson_bracket}, ::Type{_}) where {_}
    [OperatorRecipe(:diff_x), OperatorRecipe(:diff_y), OperatorRecipe(:quadratic_term)]
end

function build_operator(::Val{:poisson_bracket}, domain::Domain; diff_x, diff_y,
                        quadratic_term, kwargs...)
    PoissonBracket(domain, diff_x, diff_y, quadratic_term)
end

# -------------------------------------- Main Method ---------------------------------------

# In-place
function (op::PoissonBracket)(out::AbstractArray, u::AbstractArray, v::AbstractArray)
    @unpack tmp, qt_left, qt_right, diff_x, diff_y, quadratic_term = op
    diff_x(out, u)
    diff_y(tmp, v)
    quadratic_term(qt_left, out, tmp)
    diff_x(out, v)
    diff_y(tmp, u)
    quadratic_term(qt_right, out, tmp)

    out .= qt_left .- qt_right
end

# Out-of-place
(op::PoissonBracket)(u::AbstractArray, v::AbstractArray) = op(similar(u), u, v)