# ------------------------------------------------------------------------------------------
#                                         GPU Tests
# ------------------------------------------------------------------------------------------

# Not part of the default test run, enable with ADVECTRA_TEST_GPU=true. CUDA.jl and Plots.jl
# are not test dependencies and must be available in the active environment, e.g.
#   ADVECTRA_TEST_GPU=true julia --project -e 'using Pkg; Pkg.test()'
# after adding CUDA and Plots to test/Project.toml locally.

using Test
using CUDA

const GPU_TEST_FILES = ["cfl_tests.jl",
                        "COM_tests.jl",
                        "energy_integrals_tests.jl",
                        "fluxes_tests.jl",
                        "probe_tests.jl",
                        "spectral_tests.jl",
                        "vorticity_tests.jl"]

@testset "GPU" begin
    if CUDA.functional()
        # These files have no assertions yet, they only check that the code runs
        @testset "$file" for file in GPU_TEST_FILES
            Base.include(Module(), joinpath(@__DIR__, file))
        end
    else
        @warn "ADVECTRA_TEST_GPU=true, but CUDA is not functional. Skipping GPU tests."
    end
end
