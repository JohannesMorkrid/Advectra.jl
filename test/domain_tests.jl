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

@testset "FFT Round Trip" for a_domain in Domain_set
    # Create a random physical field on the correct memory type
    T = Advectra.get_precision(a_domain)
    phys_in = rand(T, size(a_domain)...) |> a_domain.MemoryType
    
    fwd = Advectra.get_fwd(a_domain)
    bwd = Advectra.get_bwd(a_domain)
    
    # Physical -> Spectral -> Physical
    spec = fwd * phys_in
    phys_out = bwd * spec
    
    @test Array(phys_in) ≈ Array(phys_out)
end

@testset "Constructor Failures" begin
    # Test incorrect MemoryType (passing a parameterized type)
    @test_throws ArgumentError Domain(64; MemoryType=Array{Float64})
    
    # Test incorrect precision
    @test_throws ArgumentError Domain(64; precision=String)
    
    # Test non-positive dimensions
    @test_throws ArgumentError Domain(-64; L=1.0)
end

@testset "Precision Stability" begin
    for T in [Float32, Float64]
        for real_tr in [true, false]
            d = Domain(32; precision=T, real_transform=real_tr)
            
            @test eltype(d.x) === T
            @test eltype(d.kx) === T
            
            fwd = Advectra.get_fwd(d)
            
            # If it's a real transform, plan expects T (Float)
            # If it's a complex transform, plan expects Complex{T}
            expected_plan_eltype = real_tr ? T : Complex{T}
            @test eltype(fwd) === expected_plan_eltype
            
            # The result of a forward transform is ALWAYS complex
            # Note: for complex transforms, the input must be complex
            input_type = real_tr ? T : Complex{T}
            phys_tmp = zeros(input_type, size(d)...) |> d.MemoryType
            spec_tmp = fwd * phys_tmp
            @test eltype(spec_tmp) === Complex{T}
        end
    end
end

@testset "Domain Offsets" begin
    x0, y0 = 10.0, -5.0
    Lx, Ly = 2.0, 2.0
    Nx, Ny = 10, 10
    d = Domain(Nx, Ny; Lx=Lx, Ly=Ly, x0=x0, y0=y0)
    
    @test first(d.x) ≈ x0
    @test first(d.y) ≈ y0
    @test last(d.x) ≈ (x0 + Lx - d.dx)
end

@testset "IO and Show" begin
    d = Domain(32)
    @test_nowarn show(devnull, MIME"text/plain"(), d)
    
    # Test compact mode used in arrays
    @test_nowarn show(IOContext(devnull, :compact => true), d)
end
# Documentation should explain operators

# Should have operators._domain perhaps, incase user wants to use domain info in rhs
