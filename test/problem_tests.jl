# ------------------------------------------------------------------------------------------
#                                  SpectralODEProblem Tests
# ------------------------------------------------------------------------------------------

using Test
using Advectra
import Advectra: isinplace

const domain = Domain(8; L=1)
const u0 = initial_condition(gaussian, domain)
Linear!(du, u, operators, p, t) = du .= 0
NonLinear!(du, u, operators, p, t) = du .= 0
Linear(u, operators, p, t) = zero(u)
NonLinear(u, operators, p, t) = zero(u)
problem(; kwargs...) = SpectralODEProblem(Linear!, NonLinear!, u0, domain, [0.0, 1.0];
                                          kwargs...)

@testset "Operator recipes" begin
    @test keys(problem().operators) ==
          (:diff_x, :diff_y, :laplacian, :solve_phi, :poisson_bracket)
    @test issubset((:quadratic_term, :spectral_exp, :reciprocal, :grad_dot_grad),
                   keys(problem(; operators=:all).operators))
    @test isempty(problem(; operators=:none).operators)
    @test_throws ErrorException problem(; operators=:unknown)

    # Additional operators are added under their alias
    biharmonic = OperatorRecipe(:laplacian; order=2, alias=:biharmonic)
    operators = problem(; operators=:none, additional_operators=[biharmonic]).operators
    @test keys(operators) == (:biharmonic,)
    @test operators.biharmonic.coeffs ≈ problem().operators.laplacian.coeffs .^ 2
end

@testset "Precision conversion" begin
    domain32 = Domain(8; L=1, precision=Float32)
    prob = SpectralODEProblem(Linear!, NonLinear!, initial_condition(gaussian, domain32),
                              domain32, [0, 1]; dt=0.1,
                              p=(ν=0.1, n=3, v=[1.0, 2.0], name="blob"))
    @test prob.p.ν isa Float32
    @test prob.p.n isa Int # Only floating point parameters are converted
    @test eltype(prob.p.v) == Float32
    @test prob.p.name == "blob"
    @test prob.dt isa Float32
    @test eltype(prob.tspan) == Float32
    @test eltype(prob.u0_hat) == ComplexF32
    @test_throws "exactly two" SpectralODEProblem(Linear!, NonLinear!, u0, domain,
                                                  [0.0, 0.5, 1.0])
end

@testset "In-place detection" begin
    @test isinplace(problem()) isa Val{true}
    prob = SpectralODEProblem(Linear, NonLinear, u0, domain, [0.0, 1.0])
    @test isinplace(prob) isa Val{false}
    @test_throws "Mismatch" SpectralODEProblem(Linear!, NonLinear, u0, domain, [0.0, 1.0])
    @test_throws "valid signature" isinplace((u, t) -> u)

    # Without a linear operator, it is assumed zero with the matching signature
    @test isinplace(SpectralODEProblem(NonLinear!, u0, domain, [0.0, 1.0])) isa Val{true}
    @test isinplace(SpectralODEProblem(NonLinear, u0, domain, [0.0, 1.0])) isa Val{false}
end

@testset "show" begin
    output = sprint(show, MIME"text/plain"(), problem(; p=(ν=0.1,), dt=0.05))
    @test occursin("SpectralODEProblem(", output)
    @test occursin("dt=0.05", output)
    @test occursin("p=(ν = 0.1,)", output)
    @test occursin("in-place: true", output)
    @test occursin("remove_modes: remove_nothing!", output)
    @test occursin("Domain(Nx:8", output)
end
