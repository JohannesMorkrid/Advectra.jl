# ------------------------------------------------------------------------------------------
#                                 Diagnostic Framework Tests
# ------------------------------------------------------------------------------------------

using Test
using Advectra
import Advectra: required_operators

@testset "Required operators" begin
    @testset "CFL" begin
        exb = [OperatorRecipe(:diff_x), OperatorRecipe(:diff_y), OperatorRecipe(:solve_phi)]
        @test required_operators(@diagnostics [cfl]) == exb
        @test required_operators(@diagnostics [cfl(; velocity=:ExB)]) == exb
        @test isempty(required_operators(@diagnostics [cfl(; velocity=:burger)]))

        # Operators needed by the diagnostic are built even when none are requested
        domain = Domain(16; L=10)
        u0 = initial_condition(isolated_blob, domain)
        Linear(du, u, operators, p, t) = du .= 0
        NonLinear(du, u, operators, p, t) = du .= 0
        prob = SpectralODEProblem(Linear, NonLinear, u0, domain, [0.0, 0.01]; dt=1e-3,
                                  operators=:none,
                                  diagnostics=@diagnostics [cfl(; silent=true)])
        @test issubset((:diff_x, :diff_y, :solve_phi), keys(prob.operators))
    end
end
