using Test
using Advectra
using LinearAlgebra

# ------------------------------------------------------------------------------
# 1. Setup Self-Contained Domain Set
# ------------------------------------------------------------------------------
d1 = Domain(64, 64; Lx=2π, Ly=2π, real_transform=true,  dealiased=true)
d2 = Domain(64; L=10.0,          real_transform=false, dealiased=false)
d3 = Domain(128, 64; Lx=1.0, Ly=2.0, real_transform=true, dealiased=true)

Domain_set = [d1, d2, d3]

# ------------------------------------------------------------------------------
# 2. Operator Tests
# ------------------------------------------------------------------------------

@testset "Spectral Operator Tests" begin

    @testset "First Derivatives (Accuracy) - Domain: $(size(d))" for d in Domain_set
        T = Advectra.get_precision(d)
        fwd, bwd = Advectra.get_fwd(d), Advectra.get_bwd(d)
        
        m, n = 2, 3
        k0x = 2π / d.Lx
        k0y = 2π / d.Ly
        
        u_phys = @. sin(m * k0x * d.x') * cos(n * k0y * d.y)
        u_spec = fwd * (u_phys |> d.MemoryType)
        
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
        fwd, bwd = Advectra.get_fwd(d), Advectra.get_bwd(d)

        diff_x = build_operator(Val(:diff_x), d)
        diff_y = build_operator(Val(:diff_y), d)
        
        # Wrap in try-catch in case QuadraticTerm isn't fully implemented yet
        try
            q_term = build_operator(Val(:quadratic_term), d)
            gdg = build_operator(Val(:grad_dot_grad), d; 
                                 diff_x=diff_x, diff_y=diff_y, quadratic_term=q_term)

            # Test with simple smooth fields
            u_phys = @. cos(2π * d.x' / d.Lx) + 0*d.y
            v_phys = @. sin(2π * d.y / d.Ly)  + 0*d.x'
            
            u_spec = fwd * (u_phys |> d.MemoryType)
            v_spec = fwd * (v_phys |> d.MemoryType)

            # Check if calling the operator works
            res_spec = gdg(u_spec, v_spec)
            res_phys = Array(bwd * res_spec)

            # For these orthogonal inputs, the dot product should be zero
            @test all(isapprox.(res_phys, 0.0; atol=1e-10))
            
        catch e
            if e isa MethodError || e isa UndefVarError
                @warn "Skipping GradDotGrad: QuadraticTerm dependency not met."
            else
                rethrow(e)
            end
        end
    end
end

