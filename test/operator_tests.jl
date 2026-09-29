using Test
using Advectra
using LinearAlgebra
import Advectra: SpectralConstant

# ------------------------------------------------------------------------------
# 1. Setup Self-Contained Domain Set
# ------------------------------------------------------------------------------
d1 = Domain(64, 64; Lx=2π, Ly=2π, real_transform=true, dealiased=true)
d2 = Domain(64; L=10.0, real_transform=false, dealiased=false)
d3 = Domain(128, 64; Lx=1.0, Ly=2.0, real_transform=true, dealiased=true)

Domain_set = [d1, d2, d3]

# ------------------------------------------------------------------------------
# 2. Operator Tests
# ------------------------------------------------------------------------------

@testset "Spectral Operator Tests" begin
    @testset "First Derivatives (Accuracy) - Domain: $(size(d))" for d in Domain_set
        T = Advectra.get_precision(d)
        fwd, bwd = Advectra.fwd_plan(d), Advectra.bwd_plan(d)

        m, n = 2, 3
        k0x = 2π / d.Lx
        k0y = 2π / d.Ly

        u_phys = @. sin(m * k0x * d.x') * cos(n * k0y * d.y)
        u_spec = fwd * (u_phys |> array_wrapper(d))

        # --- Test ∂x ---
        dx_op = build_operator(Val(:diff_x), d)
        du_spec = dx_op(u_spec)
        du_phys = Array(bwd * du_spec)

        expected_dx = @. (m * k0x) * cos(m * k0x * d.x') * cos(n * k0y * d.y)
        @test isapprox(du_phys, expected_dx; atol=1e-5)

        # --- Test ∂y ---
        dy_op = build_operator(Val(:diff_y), d)
        dv_spec = dy_op(u_spec)
        dv_phys = Array(bwd * dv_spec)

        expected_dy = @. -(n * k0y) * sin(m * k0x * d.x') * sin(n * k0y * d.y)
        @test isapprox(dv_phys, expected_dy; atol=1e-5)
    end

    @testset "Laplacian Consistency - Domain: $(size(d))" for d in Domain_set
        L_op = build_operator(Val(:laplacian), d)
        Dxx_op = build_operator(Val(:diff_xx), d)
        Dyy_op = build_operator(Val(:diff_yy), d)

        # Verify ∇² = ∂xx + ∂yy using the .coeffs field
        @test L_op.coeffs ≈ (Dxx_op.coeffs .+ Dyy_op.coeffs)

        @test all(real.(L_op.coeffs) .<= 0)
        @test all(isapprox.(imag.(L_op.coeffs), 0.0; atol=1e-12))
    end

    @testset "GradDotGrad Operator - Domain: $(size(d))" for d in Domain_set
        T = Advectra.get_precision(d)
        fwd, bwd = Advectra.fwd_plan(d), Advectra.bwd_plan(d)

        diff_x = build_operator(Val(:diff_x), d)
        diff_y = build_operator(Val(:diff_y), d)
        q_term = build_operator(Val(:quadratic_term), d)
        gdg = build_operator(Val(:grad_dot_grad), d;
                             diff_x=diff_x, diff_y=diff_y, quadratic_term=q_term)

        k0x = 2π / d.Lx
        k0y = 2π / d.Ly

        # Orthogonal gradients: ∇u = (∂u/∂x, 0) and ∇v = (0, ∂v/∂y) ⇒ ∇u⋅∇v = 0
        u_phys = @. cos(k0x * d.x') + 0 * d.y
        v_phys = @. sin(k0y * d.y) + 0 * d.x'
        res_phys = Array(bwd * gdg(fwd * (u_phys |> array_wrapper(d)),
                                   fwd * (v_phys |> array_wrapper(d))))
        @test all(isapprox.(res_phys, 0.0; atol=1e-10))

        # u = sin(k0x x), v = sin(k0x x) + cos(k0y y) ⇒ ∇u⋅∇v = k0x²cos²(k0x x)
        u_phys = @. sin(k0x * d.x') + 0 * d.y
        v_phys = @. sin(k0x * d.x') + cos(k0y * d.y)
        res_phys = Array(bwd * gdg(fwd * (u_phys |> array_wrapper(d)),
                                   fwd * (v_phys |> array_wrapper(d))))
        expected = @. k0x^2 * cos(k0x * d.x')^2 + 0 * d.y
        @test isapprox(res_phys, expected; atol=1e-10)
    end
end


# ------------------------------------------------------------------------------
# 3. Non-linear Operators and Operator Algebra
# ------------------------------------------------------------------------------

@testset "Non-linear operators - dealiased=$dealiased" for dealiased in (true, false)
    d = Domain(16, 12; Lx=2π, Ly=4π, dealiased)
    x, y = d.x', d.y
    F(u) = fwd_plan(d) * u
    B(u_hat) = bwd_plan(d) * u_hat
    q = build_operator(Val(:quadratic_term), d)

    @testset "Quadratic term" begin
        u = @. sin(x) + cos(y / 2)
        v = @. cos(x) + 0 * y
        @test B(q(F(u), F(v))) ≈ u .* v
    end

    @testset "Poisson bracket" begin
        pb = build_operator(Val(:poisson_bracket), d; diff_x=build_operator(Val(:diff_x), d),
                            diff_y=build_operator(Val(:diff_y), d), quadratic_term=q)
        # {ϕ, n} = ∂ϕ/∂x ∂n/∂y - ∂ϕ/∂y ∂n/∂x
        ϕ = @. sin(x) + 0 * y
        n = @. sin(y / 2) + 0 * x
        @test B(pb(F(ϕ), F(n))) ≈ @. cos(x) * cos(y / 2) / 2
        @test B(pb(F(n), F(ϕ))) ≈ @. -cos(x) * cos(y / 2) / 2
    end

    @testset "Spectral functions" begin
        # Smooth and positive, so the pseudo-spectral evaluation is accurate
        u = @. 1.5 + 0.3 * sin(x) * cos(y / 2)
        for (name, f) in [(:spectral_exp, exp), (:spectral_expm1, expm1),
                          (:spectral_log, log), (:reciprocal, inv)]
            op = build_operator(Val(name), d; quadratic_term=q)
            @test B(op(F(u))) ≈ f.(u) atol = 1e-5
        end
    end
end

@testset "Linear operator algebra" begin
    d = Domain(16; L=2π)
    diff_x = build_operator(Val(:diff_x), d)
    diff_y = build_operator(Val(:diff_y), d)
    @test (2 * diff_x).coeffs ≈ 2 .* diff_x.coeffs
    @test (diff_x + diff_y).coeffs ≈ diff_x.coeffs .+ diff_y.coeffs
    @test (diff_x - diff_y).coeffs ≈ diff_x.coeffs .- diff_y.coeffs
    @test (diff_x^2).coeffs ≈ build_operator(Val(:diff_xx), d).coeffs
    @test (diff_x^2 + diff_y^2).coeffs ≈ build_operator(Val(:laplacian), d).coeffs
end

@testset "Spectral constant" begin
    a, b = SpectralConstant(; val=6.0), SpectralConstant(; val=2.0)
    @test (a + b).value == 8
    @test (a - b).value == 4
    @test (a * b).value == 12
    @test (a / b).value == 3
    @test (a * 2).value == (2 * a).value == 12
    @test (a / 2).value == 3
    @test (12 / a).value == 2
    @test (-a).value == -6

    # Only the zeroth mode (first entry) is affected
    field = ones(ComplexF64, 3, 2)
    @test field + b == b + field == [3 1; 1 1; 1 1]
    @test field - b == [-1 1; 1 1; 1 1]
    @test b - field == [1 -1; -1 -1; -1 -1]

    d = Domain(16; L=2π)
    @test build_operator(Val(:spectral_constant), d; val=4).value == 4
end

@testset "Float32 precision" begin
    results = map((Float32, Float64)) do T
        d = Domain(16, 12; Lx=2π, Ly=4π, precision=T)
        x, y = d.x', d.y
        F(u) = fwd_plan(d) * u
        operators = build_operators(d; operators=:all)
        u_hat = F(@. 1.5 + 0.3 * sin(x) * cos(y / 2))
        v_hat = F(@. cos(x) + sin(y / 2))
        [operators.diff_x(u_hat), operators.laplacian(u_hat),
         operators.solve_phi(u_hat, v_hat), operators.quadratic_term(u_hat, v_hat),
         operators.poisson_bracket(u_hat, v_hat), operators.spectral_log(u_hat)]
    end
    for (result32, result64) in zip(results...)
        @test eltype(result32) == ComplexF32
        @test result32 ≈ result64 rtol = 1e-5
    end
end
