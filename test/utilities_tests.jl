# ------------------------------------------------------------------------------------------
#                                      Utilities Tests
# ------------------------------------------------------------------------------------------

using Test
using Advectra
using Statistics

# Nx ≠ Ny and Lx ≠ Ly to catch mixed up axes
const domain = Domain(16, 8; Lx=10, Ly=6)
const x, y = domain.x', domain.y

@nobroadcast constant_ic(domain; value=1.0) = fill(value, size(domain))

@testset "Initial conditions" begin
    @testset "Pointwise functions" begin
        ic(f; kwargs...) = initial_condition(f, domain; kwargs...)
        @test ic(gaussian; A=2, B=1, l=1.5) ≈ @. 1 + 2exp(-(x^2 + y^2) / (2 * 1.5^2))
        @test ic(gaussian; lx=1, ly=2, x0=1, y0=-1) ≈
              @. exp(-(x - 1)^2 / 2 - (y + 1)^2 / 8)
        @test ic(log_gaussian; A=2) ≈ @. log(1 + 2exp(-(x^2 + y^2) / 2))
        @test ic(gaussian_x; x0=1) ≈ @. exp(-(x - 1)^2 / 2) + 0 * y
        @test ic(gaussian_y; y0=1) ≈ @. exp(-(y - 1)^2 / 2) + 0 * x
        @test ic(sinusoidal; lx=10, ly=6) ≈ @. sin(2π * x / 10) * cos(2π * y / 6)
        @test ic(sinusoidal_x; L=10, N=2) ≈ @. sin(4π * x / 10) + 0 * y
        @test ic(sinusoidal_y; L=6) ≈ @. sin(2π * y / 6) + 0 * x
        @test ic(exponential_x; κ=0.5) ≈ @. exp(-0.5x) + 0 * y
        @test ic(quadratic_y) ≈ @. ifelse(abs(y) <= 1, 1 - y^2, 0.0) + 0 * x
        @test size(ic(white_noise)) == size(domain)
    end

    @testset "Blobs" begin
        blob = initial_condition(isolated_blob, domain; A=2)
        @test size(blob) == (size(domain)..., 2)
        @test blob[:, :, 1] ≈ initial_condition(gaussian, domain; A=2)
        @test iszero(blob[:, :, 2])
        blob = initial_condition(isolated_blob, domain; density=:log, ndims=3)
        @test size(blob) == (size(domain)..., 3)
        @test blob[:, :, 1] ≈ initial_condition(log_gaussian, domain)

        blob = initial_condition(isolated_temperature_blob, domain)
        @test size(blob) == (size(domain)..., 3)
        @test all(isone, blob[:, :, 1]) # Unity density
        @test iszero(blob[:, :, 2])
        @test blob[:, :, 3] ≈ initial_condition(gaussian, domain)
        blob = initial_condition(isolated_temperature_blob, domain; density=:log)
        @test iszero(blob[:, :, 1]) # log(1) = 0
        @test blob[:, :, 3] ≈ initial_condition(log_gaussian, domain)
    end

    @testset "Random modes" begin
        u = initial_condition(random_phase, domain; value=1e-3)
        @test size(u) == size(domain)
        # Zonal (ky = 0) and streamer (kx = 0) modes are removed
        @test maximum(abs, mean(u; dims=1)) < 1e-12
        @test maximum(abs, mean(u; dims=2)) < 1e-12
        @test size(initial_condition(random_phase, domain; ndims=3)) == (size(domain)..., 3)

        u = initial_condition(random_crossphased, domain)
        @test size(u) == (size(domain)..., 2)
        @test maximum(abs, mean(u[:, :, 1]; dims=1)) < 1e-12
    end

    @testset "@nobroadcast" begin
        # The domain is passed to the function instead of broadcasting over the grid
        @test initial_condition(constant_ic, domain; value=2.0) == fill(2.0, size(domain))
    end

    @testset "Precision" begin
        domain32 = Domain(16, 8; Lx=10, Ly=6, precision=Float32)
        for f in (gaussian, log_gaussian, gaussian_x, gaussian_y, sinusoidal, sinusoidal_x,
                  sinusoidal_y, exponential_x, quadratic_y, random_phase, random_crossphased,
                  isolated_blob, isolated_temperature_blob)
            @test eltype(initial_condition(f, domain32)) == Float32
        end
        @test eltype(initial_condition(white_noise, domain32)) == ComplexF32
    end
end

@testset "Mode removal - real_transform=$real_transform" for real_transform in (true, false)
    # Even sizes have Nyquist modes, odd sizes do not
    for (Nx, Ny) in ((16, 8), (15, 7))
        d = Domain(Nx, Ny; Lx=10, Ly=6, real_transform)
        u_hat = fwd_plan(d) * randn(real_transform ? Float64 : ComplexF64, Ny, Nx)
        remove(f!) = (v = copy(u_hat); f!(v, d); v)

        # Zonal modes have ky = 0 (first row), streamer modes kx = 0 (first column)
        @test iszero(remove(remove_zonal_modes!)[1, :])
        @test remove(remove_zonal_modes!)[2:end, :] == u_hat[2:end, :]
        @test iszero(remove(remove_streamer_modes!)[:, 1])
        @test remove(remove_streamer_modes!)[:, 2:end] == u_hat[:, 2:end]
        @test remove(remove_nothing!) == u_hat

        v = remove(remove_nyquist_modes!)
        @test iseven(Nx) ? iszero(v[:, Nx ÷ 2 + 1]) : v == u_hat
        @test iseven(Ny) ? iszero(v[Ny ÷ 2 + 1, :]) : v == u_hat
        @test remove_nyquist_modes! === remove_asymmetric_modes!

        # Multiple fields are handled the same way
        state_hat = cat(u_hat, u_hat; dims=3)
        remove_zonal_modes!(state_hat, d)
        @test iszero(state_hat[1, :, :])
    end
end

@testset "Other utilities" begin
    field = [1.0 2.0; 3.0 4.0]
    # Only the first (zeroth mode) entry is changed
    @test add_constant(field, 5) == [6.0 2.0; 3.0 4.0]
    @test field == [1.0 2.0; 3.0 4.0] # Out-of-place leaves the input unchanged
    out = zero(field)
    @test add_constant!(out, field, 5) === out
    @test out == [6.0 2.0; 3.0 4.0]
    @test add_constant!(field, -1) == [0.0 2.0; 3.0 4.0]

    @test logspace(-1, 2, 4) ≈ [0.1, 1, 10, 100]

    # Only meaningful while the SMTPClient extension is not loaded (see smtp_tests.jl)
    if isnothing(Base.get_extension(Advectra, :AdvectraSMTPClientExt))
        @test_throws "SMTPClient is not loaded" send_mail("Simulation done")
    end
end

@testset "spectral_sum" begin
    # Parseval: ∑|Â|² = N∑|A|², the real transform only stores half of the modes
    for (Nx, Ny) in ((16, 8), (15, 7), (16, 7))
        A = randn(Ny, Nx)
        rdomain = Domain(Nx, Ny; Lx=10, Ly=6)
        cdomain = Domain(Nx, Ny; Lx=10, Ly=6, real_transform=false)
        Ar_hat = fwd_plan(rdomain) * A
        Ac_hat = fwd_plan(cdomain) * complex.(A)
        N = length(rdomain)
        @test spectral_sum(abs2, Ar_hat, rdomain) ≈ N * sum(abs2, A)
        @test spectral_sum(abs2, Ac_hat, cdomain) ≈ N * sum(abs2, A)
        @test spectral_sum(abs2.(Ar_hat), rdomain) ≈ N * sum(abs2, A)
        @test sum(spectral_sum(abs2, Ar_hat, rdomain; dims=1)) ≈ N * sum(abs2, A)
        @test only(spectral_sum(abs2, Ar_hat, rdomain; dims=(1, 2))) ≈ N * sum(abs2, A)

        # A single row (ky = 0) or column (kx = 0) of modes
        @test spectral_sum(abs2, Ar_hat[1:1, :], rdomain) ≈ sum(abs2, Ac_hat[1:1, :])
        @test spectral_sum(abs2, Ar_hat[:, 1:1], rdomain) ≈ sum(abs2, Ac_hat[:, 1:1])
    end
end
