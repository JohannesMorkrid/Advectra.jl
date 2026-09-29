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

# TODO test the storage helpers (see TESTS_TODO D8):
# * parse_storage_limit, e.g. "1 MB", "", negative and invalid input
# * determine_sampling_strategy, validate_stride and recommend_stride, including the cases
#   where N_steps/stride is not an integer, there is not enough storage for two samples,
#   the stride is negative and the stride is 0
# * a simulation with the same name can have a new domain size
