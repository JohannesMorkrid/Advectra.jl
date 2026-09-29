using Test
using Advectra
using LinearAlgebra

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

