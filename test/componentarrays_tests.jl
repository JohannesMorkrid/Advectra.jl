# ------------------------------------------------------------------------------------------
#                                ComponentArrays Extension Tests
# ------------------------------------------------------------------------------------------

using Test
using Advectra
using ComponentArrays

@testset "ComponentArrays extension" begin
    domain = Domain(16; L=10)
    n = initial_condition(gaussian, domain)
    Ω = initial_condition(sinusoidal, domain; lx=10, ly=10)
    parameters = (ν=0.1,)
    tspan = [0.0, 0.1]

    # Same diffusion problem, once as a ComponentArray and once as a plain 3D Array
    function Linear_ca(du, u, operators, p, t)
        du.n .= p.ν .* operators.laplacian(u.n)
        du.Ω .= p.ν .* operators.laplacian(u.Ω)
    end
    Linear(du, u, operators, p, t) = du .= p.ν .* operators.laplacian(u)
    NonLinear(du, u, operators, p, t) = du .= 0

    prob_ca = SpectralODEProblem(Linear_ca, NonLinear, ComponentArray(; n, Ω), domain,
                                 tspan; p=parameters, dt=1e-3)
    prob = SpectralODEProblem(Linear, NonLinear, cat(n, Ω; dims=3), domain, tspan;
                              p=parameters, dt=1e-3)

    @testset "Spectral coefficients keep their structure" begin
        @test prob_ca.u0 isa ComponentArray
        @test prob_ca.u0_hat isa ComponentArray
        @test size(prob_ca.u0_hat.n) == spectral_size(domain)
        @test prob_ca.u0_hat.n ≈ prob.u0_hat[:, :, 1]
        @test prob_ca.u0_hat.Ω ≈ prob.u0_hat[:, :, 2]
    end

    @testset "Solution matches plain Array" begin
        mktempdir() do dir
            output_ca = Output(prob_ca; filename=joinpath(dir, "ca.h5"))
            output = Output(prob; filename=joinpath(dir, "array.h5"))
            spectral_solve(prob_ca, MSS3(), output_ca; debug=true)
            spectral_solve(prob, MSS3(), output; debug=true)

            # Fields are stored contiguously (n then Ω) in both cases
            u_ca = read(output_ca.simulation, "checkpoint/u")
            u = read(output.simulation, "checkpoint/u")
            @test u_ca ≈ vec(u)
            @test !(u ≈ prob.u0_hat) # Sanity check that the state evolved

            close(output_ca.simulation.file)
            close(output.simulation.file)
        end
    end
end
