# ------------------------------------------------------------------------------------------
#                                    Spectral Operators                                     
# ------------------------------------------------------------------------------------------

# ---------------------------------------- General -----------------------------------------

# Abstract class that may inherit from a more abstract one
abstract type SpectralOperator end

# TODO check if needed
# For broadcasting the operators for composite operators
Base.broadcastable(op::SpectralOperator) = Ref(op)

# ----------------------------------- Linear Operators -------------------------------------

# Abstract type that all Linear operators inherit from
abstract type LinearOperator{T} <: SpectralOperator end

# TODO add Algebraic rules of Linear operators
import Base: *, +, -, ^
*(a::Number, op::LinearOperator) = typeof(op)(a .* op.coeffs)

# Allows composite operators
# TODO add some sort of promotion rules
+(a::LinearOperator, b::LinearOperator) = typeof(op)(a.coeffs .+ b.coeffs)
-(a::LinearOperator, b::LinearOperator) = typeof(op)(a.coeffs .- b.coeffs)
^(op::LinearOperator, power::Number) = typeof(op)(op.coeffs .^ power)

# --------------------------------- Elementwise Operator -----------------------------------

struct ElwiseOperator{T<:AbstractArray} <: LinearOperator{T}
    coeffs::T
    order::Number

    ElwiseOperator(coeffs; order=1) = new{typeof(coeffs)}(coeffs .^ order, order)
end

# Out-of-place operator
(op::ElwiseOperator)(u::AbstractArray) = op.coeffs .* u

# In-place operator
(op::ElwiseOperator)(du::AbstractArray, u::AbstractArray) = du .= op.coeffs .* u

# To be able to use @. without applying LinearOperator to array element
import Base.Broadcast: broadcasted
broadcasted(op::ElwiseOperator, x) = broadcasted(*, op.coeffs, x)

# ------------------------------------ Matrix Operator -------------------------------------

struct MatrixOperator{T<:AbstractArray} <: LinearOperator{T}
    coeffs::T
    order::Number
end

# Out-of-place operator # TODO figure out what to do here
(op::MatrixOperator)(u::AbstractArray) = op.coeffs * u

# --------------------------------- Non-Linear Operators -----------------------------------

# Abstract type that all NonLinear operators inherit from
abstract type NonLinearOperator <: SpectralOperator end

include("quadraticTerm.jl")
include("spatialDerivatives.jl") # Technically both linear and non-linear
include("poissonBracket.jl")

# ---------------------------------------- Others ------------------------------------------

include("spectralFunctions.jl")
include("solvePhi.jl")
include("sources.jl")

# ------------------------------------ Operator Recipe -------------------------------------

struct OperatorRecipe{KwargsType<:NamedTuple}
    op::Symbol
    alias::Symbol
    kwargs::KwargsType

    function OperatorRecipe(op; alias=op, kwargs...)
        # TODO check that op is valid
        nt = NamedTuple(kwargs)
        new{typeof(nt)}(op, alias, nt)
    end
end

Base.hash(oprecipe::OperatorRecipe, h::UInt) = hash((oprecipe.op, oprecipe.kwargs), h)
function Base.:(==)(or1::OperatorRecipe, or2::OperatorRecipe)
    or1.op == or2.op && or1.kwargs == or2.kwargs
end

# function Base.get(oprecipe::OperatorRecipe, key, default)
#     hasproperty(oprecipe, key) ? getproperty(oprecipe, key) : default
# end

# ------------------------------------- Constructors ---------------------------------------

function build_operator(op::Symbol, domain::AbstractDomain; kwargs...)
    build_operator(Val(op), domain; kwargs...)
end

# TODO add @op for each operator

# Catch all method, can be overwritten with specilization
operator_dependencies(::Val{_}, ::Type{__}) where {_,__} = ()

function get_operator_recipes(operators::Symbol)
    if operators == :default
        return [OperatorRecipe(:diff_x; order=1),
                OperatorRecipe(:diff_y; order=1),
                OperatorRecipe(:laplacian; order=1),
                OperatorRecipe(:solve_phi),
                OperatorRecipe(:poisson_bracket)]
    elseif operators == :SOL
        return [OperatorRecipe(:diff_x),
                OperatorRecipe(:diff_y),
                OperatorRecipe(:laplacian),
                OperatorRecipe(:solve_phi),
                OperatorRecipe(:poisson_bracket),
                OperatorRecipe(:quadratic_term)]
    elseif operators == :all
        return [OperatorRecipe(:diff_x; order=1),
                OperatorRecipe(:diff_y; order=1),
                OperatorRecipe(:laplacian; order=1),
                OperatorRecipe(:solve_phi),
                OperatorRecipe(:poisson_bracket),
                OperatorRecipe(:quadratic_term),
                OperatorRecipe(:spectral_log),
                OperatorRecipe(:spectral_exp),
                OperatorRecipe(:spectral_expm1),
                OperatorRecipe(:reciprocal),
                OperatorRecipe(:grad_dot_grad)]
    elseif operators == :none
        OperatorRecipe[]
    else
        error()
    end
end

# --------------------------------- OperatorRecipe Macro -----------------------------------

# TODO implement
macro op()
end

## -------------------------------- TESTING BELOW ------------------------------------------

# @op ∂xx = diff_x(order=2) => OperatorRecipe(:diff_x, order=2, alias=∂xx)

# function prepare_operators(::Type{Domain}, operators, kx, ky, Nx, Ny; MemoryType=MemoryType,
#     precision=precision, real_transform=real_transform, dealiased=dealiased)

#     cache = Dict{Symbol,SpectralOperator}()

#     # Figure out aliases
#     #[:∂x, :∂y, :laplacian, :poisson_bracket] -> [:diff_x, :diff_y, :laplacian, :poisson_bracket]

#     # Link operator to alias
#     #[:∂x, :∂y, :laplacian, :poisson_bracket], [:diff_x, :diff_y, :laplacian, :poisson_bracket], [:diff_x, :diff_y, :laplacian, :quadratic_term, :poisson_bracket], [ElwiseOperator, ElwiseOperator, ElwiseOperator, QuadraticTerm, PoissonBracket]
#     #-> spectral_operators = (:∂x=ElwiseOperator, :∂y=ElwiseOperator, :laplacian=ElwiseOperator, :poisson_bracket=PoissonBracket)

# end