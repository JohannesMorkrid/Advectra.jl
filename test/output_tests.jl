# ------------------------------------------------------------------------------------------
#                                        Output Tests
# ------------------------------------------------------------------------------------------

using Test
using Advectra
using HDF5

@testset "Output" begin
    domain = Domain(16; L=10)
    u0 = initial_condition(isolated_blob, domain)
    Linear(du, u, operators, p, t) = du .= 0
    NonLinear(du, u, operators, p, t) = du .= 0

    @testset "Problem without parameters" begin
        prob = SpectralODEProblem(Linear, NonLinear, u0, domain, [0.0, 0.01]; dt=1e-3)
        mktempdir() do dir
            output = Output(prob; filename=joinpath(dir, "no parameters.h5"))
            @test spectral_solve(prob, MSS3(), output; debug=true) isa Output
            @test read_attribute(output.simulation, "Nx") == domain.Nx
            close(output.simulation.file)

            output = Output(prob; filename=joinpath(dir, "named.h5"),
                            simulation_name=:parameters)
            @test haskey(output.simulation.file, "no parameters")
            close(output.simulation.file)
        end
    end

    @testset "Without HDF5 storage" begin
        prob = SpectralODEProblem(Linear, NonLinear, u0, domain, [0.0, 0.01]; p=(ν=0.1,),
                                  dt=1e-3)
        output = Output(prob; store_hdf=false)
        @test isnothing(output.simulation)
        @test spectral_solve(prob, MSS3(), output; debug=true) isa Output
    end
end

@testset "HDF5 layout" begin
    domain = Domain(8; L=1)
    u0 = initial_condition(isolated_blob, domain)
    Linear(du, u, operators, p, t) = du .= 0
    NonLinear(du, u, operators, p, t) = du .= 0
    # 8 steps, dt is exact in floating point
    problem(diagnostics) = SpectralODEProblem(Linear, NonLinear, u0, domain, [0.0, 0.5];
                                              p=(ν=0.1,), dt=0.0625,
                                              diagnostics=diagnostics)
    # A 8×8 Float64 sample is 512 bytes, so 1.5 KB only fits the first and last sample
    prob = problem(@diagnostics [sample_density(; stride=4),
                                 sample_vorticity(; storage_limit="1.5 KB")])
    double_density!(U) = (selectdim(U, 3, 1) .*= 2; U)

    mktempdir() do dir
        output = Output(prob; filename=joinpath(dir, "layout.h5"),
                        physical_transform=double_density!)
        spectral_solve(prob, MSS3(), output; debug=true)
        simulation = output.simulation

        # Named after the file by default, with the problem as attributes
        @test HDF5.name(simulation) == "/layout"
        for (key, value) in [("dt", 0.0625), ("Nx", 8), ("Ny", 8), ("Lx", 1.0),
                             ("real_transform", true), ("ν", 0.1)]
            @test read_attribute(simulation, key) == value
        end

        @test read_attribute(simulation["Density"], "metadata") == "Sampled density field"
        @test simulation["Density/t"][:] == [0.0, 0.25, 0.5]
        @test size(simulation["Density/data"]) == (8, 8, 3)
        @test simulation["Density/data"][:, :, 1] ≈ 2 .* u0[:, :, 1]
        @test simulation["Vorticity/t"][:] == [0.0, 0.5]
        @test haskey(simulation, "checkpoint")
        close(simulation.file)

        for (simulation_name, group) in [("custom", "custom"), (:parameters, "ν=0.1")]
            output = Output(prob; filename=joinpath(dir, "names.h5"), simulation_name)
            @test HDF5.name(output.simulation) == "/" * group
            close(output.simulation.file)
        end

        # The total storage limit covers all samples: 3 × 512 bytes
        prob = problem(@diagnostics [sample_density(; stride=4)])
        @test_throws ErrorException Output(prob; filename=joinpath(dir, "limit.h5"),
                                           storage_limit="1.5 KB")
        output = Output(prob; filename=joinpath(dir, "limit.h5"), storage_limit="1.6 KB")
        close(output.simulation.file)
    end
end

@testset "Checkpoint and resume" begin
    domain = Domain(8; L=10)
    u0 = initial_condition(isolated_blob, domain)
    Linear(du, u, operators, p, t) = du .= p.ν .* operators.laplacian(u)
    function NonLinear(du, u, operators, p, t)
        n, Ω = eachslice(u; dims=3)
        dn, dΩ = eachslice(du; dims=3)
        ϕ = operators.solve_phi(n, Ω)
        dn .= .-operators.poisson_bracket(ϕ, n)
        dΩ .= .-operators.poisson_bracket(ϕ, Ω)
    end
    problem(T) = SpectralODEProblem(Linear, NonLinear, u0, domain, [0.0, T]; p=(ν=0.1,),
                                    dt=0.0625, diagnostics=@diagnostics([sample_density]))

    mktempdir() do dir
        full = Output(problem(1.0); filename=joinpath(dir, "full.h5"))
        spectral_solve(problem(1.0), MSS3(), full; debug=true)

        # Solve to T/2, then resume from the checkpoint to T
        half = Output(problem(0.5); filename=joinpath(dir, "resumed.h5"))
        spectral_solve(problem(0.5), MSS3(), half; debug=true)
        @test read(half.simulation, "checkpoint/step") == 8
        close(half.simulation.file)

        resumed = @test_logs (:warn,) (:info,) Output(problem(1.0);
                                                     filename=joinpath(dir, "resumed.h5"),
                                                     resume=true)
        spectral_solve(problem(1.0), MSS3(), resumed; debug=true)
        @test read(resumed.simulation, "checkpoint/step") == 16
        @test read(resumed.simulation, "checkpoint/u") ≈ read(full.simulation, "checkpoint/u")
        @test resumed.simulation["Density/t"][:] == full.simulation["Density/t"][:]
        @test resumed.simulation["Density/data"][:, :, :] ≈
              full.simulation["Density/data"][:, :, :]
        close(resumed.simulation.file)
        close(full.simulation.file)

        # Resuming with a different domain is rejected
        domain2 = Domain(16; L=10)
        prob2 = SpectralODEProblem(Linear, NonLinear,
                                   initial_condition(isolated_blob, domain2), domain2,
                                   [0.0, 1.0]; p=(ν=0.1,), dt=0.0625)
        # Silences the expected warning about only checking a subset of the parameters
        quiet(f) = Base.CoreLogging.with_logger(f, Base.CoreLogging.NullLogger())
        @test_throws ErrorException quiet() do
            Output(prob2; filename=joinpath(dir, "resumed.h5"), resume=true)
        end
    end
end

# TODO test the storage helpers (see TESTS_TODO D8):
# * parse_storage_limit, e.g. "1 MB", "", negative and invalid input
# * determine_sampling_strategy, validate_stride and recommend_stride, including the cases
#   where N_steps/stride is not an integer, there is not enough storage for two samples,
#   the stride is negative and the stride is 0
# * a simulation with the same name can have a new domain size
