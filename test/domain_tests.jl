using Test
using Advectra
import Advectra: Domain, lengths, wave_vectors, differential_elements,
                 domain_kwargs, spectral_size, spectral_length, area, differential_area,
                 get_points

d1 = Domain(256, 256; Lx=1, Ly=1)
d2 = Domain(256; L=1)
d3 = Domain(256, 256; dealiased=true)
d4 = Domain(128, 256; Lx=2, Ly=1, real_transform=true)
d5 = Domain(128, 256; Lx=2, Ly=1, real_transform=false)
d6 = Domain(64, 128; Lx=1, Ly=1, dealiased=false)
d7 = Domain(128, 64; Lx=1, Ly=1, dealiased=false, real_transform=false)
Domain_set = [d1, d2, d3, d4, d5, d6, d7]

@testset "Domain outputs" begin
    # Checking default Domain construction
    d1 = Domain(256, 256; Lx=1, Ly=1)
    @test isa(d1, Domain)
    @test d1.Nx == 256
    @test d1.Ny == 256
    @test d1.Lx == 1.0
    @test d1.Ly == 1.0

    # SquareDomain should produce same Nx==Ny and matching lengths when given L
    d2 = Domain(256; L=1)
    @test d2.Nx == d2.Ny && d2.Lx == d2.Ly

    # Flags propagated
    d3 = Domain(256, 256; dealiased=true)
    @test getproperty(d3, :dealiased) === true

    d4 = Domain(128, 256; Lx=2, Ly=1, real_transform=true)
    @test d4.Nx == 128 && d4.Ny == 256
    @test getproperty(d4, :real_transform) === true

    @test lengths(d1) == (1.0, 1.0)
end

@testset "Wave vectors" for a_domain in Domain_set
    # indices from -Nx/2 .. Nx/2-1 as integers
    i = vcat(collect(0:(a_domain.Nx÷2-1)), collect(-a_domain.Nx÷2:-1))
    kx = (2 * π / a_domain.Lx) .* i
    @test wave_vectors(a_domain)[2] ≈ kx

    if a_domain.real_transform # for real transforms, kx only has non-negative values
        ky_real = (2 * π / a_domain.Ly) .* collect(0:(a_domain.Ny÷2))
        @test wave_vectors(a_domain)[1] ≈ ky_real
    else
        j = vcat(collect(0:(a_domain.Ny÷2-1)), collect(-a_domain.Ny÷2:-1))
        ky = (2 * π / a_domain.Ly) .* j
        # to specify tolerances for floating point comparisons
        @test isapprox(wave_vectors(a_domain)[1], ky; rtol=1e-12, atol=1e-12)
    end
end

@testset "Differential elements" begin
    d1 = Domain(256, 256; Lx=1, Ly=1)
    differential_elements(d1) == (d1.dx, d1.dy)
    @test diff(d1.x)[end] ≈ d1.dx && diff(d1.y)[end] ≈ d1.dy
end

@testset "Domain keyword arguments" for a_domain in Domain_set
    kwargs = domain_kwargs(a_domain)
    @test haskey(kwargs, :real_transform) && haskey(kwargs, :dealiased)
    @test kwargs[1] === a_domain.real_transform
    @test kwargs[2] === a_domain.dealiased
end

@testset "Spectral size and length" for a_domain in Domain_set
    spec_size = spectral_size(a_domain)
    spec_length = spectral_length(a_domain)

    if a_domain.real_transform
        expected_size = (a_domain.Ny ÷ 2 + 1, a_domain.Nx)
        expected_length = (a_domain.Ny ÷ 2 + 1) * a_domain.Nx
    else
        expected_size = (a_domain.Ny, a_domain.Nx)
        expected_length = a_domain.Ny * a_domain.Nx
    end
    @test spec_size == expected_size
    @test spec_length == expected_length
end

@testset "Area and Differential Area" for a_domain in Domain_set
    total_area = area(a_domain)
    diff_area = differential_area(a_domain)

    expected_area = a_domain.Lx * a_domain.Ly
    expected_diff_area = a_domain.dx * a_domain.dy

    @test isapprox(total_area, expected_area; rtol=1e-12, atol=1e-12)
    @test isapprox(diff_area, expected_diff_area; rtol=1e-12, atol=1e-12)
end

# Documentation should explain operators

# Should have operators._domain perhaps, incase user wants to use domain info in rhs