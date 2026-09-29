# ------------------------------------------------------------------------------------------
#                                  Profile Diagnostic Tests
# ------------------------------------------------------------------------------------------

using Test
using Advectra
using Statistics
import Advectra: build_diagnostic, build_operator, required_operators

@testset "Profile diagnostics" begin
    # Nx ≠ Ny and Lx ≠ Ly to catch mixed up axes
    domain = Domain(32, 16; Lx=10, Ly=6)
    x, y = domain.x', domain.y

    # sin(2πy/Ly) averages to zero along y and gaussian_x has no y-dependence
    n = @. exp(-x^2 / 2) + sin(2π * y / domain.Ly)
    Ω = @. cos(2π * x / domain.Lx) * (1 + y^2)
    state = cat(n, Ω; dims=3)
    prob = (; domain=domain,
            operators=(; diff_x=build_operator(Val(:diff_x), domain),
                       diff_y=build_operator(Val(:diff_y), domain),
                       solve_phi=build_operator(Val(:solve_phi), domain)))

    @testset "Density and vorticity profiles" begin
        @test radial_density_profile(state, prob, 0.0) ≈ vec(exp.(-x .^ 2 / 2))
        @test poloidal_density_profile(state, prob, 0.0) ≈
              mean(exp.(-x .^ 2 / 2)) .+ sin.(2π * y / domain.Ly)
        @test length(radial_vorticity_profile(state, prob, 0.0)) == domain.Nx
        @test poloidal_vorticity_profile(state, prob, 0.0) ≈ zeros(domain.Ny) atol = 1e-12
    end

    @testset "Radial flux profile" begin
        # Single mode: Ω = cos(kx x)cos(ky y) ⇒ ϕ = -Ω/k² ⇒ v_x = -∂ϕ/∂y, so with
        # n = cos(kx x)sin(ky y): nv_x = -(ky/k²)cos²(kx x)sin²(ky y)
        kx, ky = 2π / domain.Lx, 2π / domain.Ly
        n_mode = @. cos(kx * x) * sin(ky * y)
        Ω_mode = @. cos(kx * x) * cos(ky * y)
        state_hat = spectral_transform(cat(n_mode, Ω_mode; dims=3), fwd_plan(domain))

        Γ = radial_flux_profile(state_hat, prob, 0.0)
        @test Γ ≈ vec(@. -ky / (kx^2 + ky^2) * cos(kx * x)^2 / 2)
        # Averaging the profile over x must give the average radial flux
        @test mean(Γ) ≈ real(build_diagnostic(Val(:radial_flux))(state_hat, prob, 0.0))
    end

    @testset "Construction through @diagnostics" begin
        recipes = @diagnostics [radial_density_profile, poloidal_density_profile,
                                radial_vorticity_profile, poloidal_vorticity_profile,
                                radial_flux_profile]
        @test length(required_operators(recipes)) == 2

        Linear(du, u, operators, p, t) = du .= 0
        NonLinear(du, u, operators, p, t) = du .= 0
        # TODO p is required, NullParameters can not be written to HDF5 (see TESTS_TODO A7)
        sim = SpectralODEProblem(Linear, NonLinear, state, domain, [0.0, 0.01];
                                 p=(ν=0.0,), dt=1e-3, diagnostics=recipes)
        mktempdir() do dir
            output = Output(sim; filename=joinpath(dir, "profiles.h5"))
            @test spectral_solve(sim, MSS3(), output; debug=true) isa Output
            for name in ["Radial density profile", "Poloidal density profile",
                         "Radial vorticity profile", "Poloidal vorticity profile",
                         "Radial flux profile"]
                @test haskey(output.simulation, name)
            end
            close(output.simulation.file)
        end
    end
end
