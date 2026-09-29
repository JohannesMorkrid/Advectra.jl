# ------------------------------------------------------------------------------------------
#                                 Diagnostic Framework Tests
# ------------------------------------------------------------------------------------------

using Test
using Advectra
import Advectra: required_operators, build_diagnostic, build_operator, sample_density,
                 sample_potential

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

@testset "@diagnostics macro" begin
    # Vector, block (with and without commas) and single diagnostic forms
    for recipes in (@diagnostics([progress, sample_density(; stride=2)]),
                    @diagnostics(begin
                                     progress
                                     sample_density(; stride=2)
                                 end),
                    @diagnostics(begin
                                     progress, sample_density(; stride=2)
                                 end))
        @test length(recipes) == 2
        @test recipes[1].method === progress && recipes[1].stride == -1
        @test recipes[2].method === sample_density && recipes[2].stride == 2
    end

    recipe = only(@diagnostics sample_density(; stride=4, storage_limit="1 MB", extra=1))
    @test (recipe.stride, recipe.storage_limit) == (4, "1 MB")
    @test recipe.kwargs == (; extra=1)

    # Keyword values and methods come from the caller's scope
    stride = 5
    my_diagnostic(state, prob, time) = sum(state)
    recipes = @diagnostics [sample_density(; stride), my_diagnostic(; stride=2stride)]
    @test recipes[1].stride == 5
    @test recipes[2].method === my_diagnostic && recipes[2].stride == 10

    # Aliases are not supported
    @test_throws "Aliases" macroexpand(@__MODULE__, :(@diagnostics [alias = progress]))

    # Operators required by several diagnostics are all collected
    recipes = @diagnostics [sample_potential, radial_flux, progress]
    @test required_operators(recipes) ==
          [OperatorRecipe(:solve_phi), OperatorRecipe(:solve_phi), OperatorRecipe(:diff_y)]
end

@testset "Sample diagnostics" begin
    domain = Domain(16, 8; Lx=10, Ly=6)
    state = randn(8, 16, 3)
    state_hat = spectral_transform(state, fwd_plan(domain))
    solve_phi = build_operator(Val(:solve_phi), domain)
    prob = (; domain=domain, operators=(; solve_phi))

    for (method, name, field) in [(:sample_density, "Density", 1),
                                  (:sample_vorticity, "Vorticity", 2),
                                  (:sample_temperature, "Temperature", 3)]
        diagnostic = build_diagnostic(Val(method))
        @test diagnostic.name == name
        @test !diagnostic.assumes_spectral_state
        @test diagnostic(state, prob, 0.0) == state[:, :, field]
    end

    diagnostic = build_diagnostic(Val(:sample_potential))
    @test diagnostic.name == "Potential"
    @test diagnostic.assumes_spectral_state
    ϕ_hat = solve_phi(state_hat[:, :, 1], state_hat[:, :, 2])
    @test diagnostic(state_hat, prob, 0.0) ≈ bwd_plan(domain) * ϕ_hat
    @test required_operators(@diagnostics [sample_potential]) == [OperatorRecipe(:solve_phi)]
end
