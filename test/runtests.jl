# Render plots off-screen, must be set before Plots is loaded
ENV["GKSwstype"] = "100"

using Test

const TEST_FILES = ["domain_tests.jl",
                    "operator_tests.jl",
                    "problem_tests.jl",
                    "schemes_tests.jl",
                    "spectrum_tests.jl",
                    "profile_tests.jl",
                    "diagnostics_tests.jl",
                    "display_tests.jl",
                    "progress_tests.jl",
                    "output_tests.jl",
                    "integration_tests.jl",
                    "componentarrays_tests.jl"]

# Show the testsets inside each file as well, e.g. ADVECTRA_TEST_VERBOSE=true
const VERBOSE = get(ENV, "ADVECTRA_TEST_VERBOSE", "false") == "true"

# Each file is included in its own module, so files can not depend on each other's imports
@testset "Advectra" verbose=true begin
    @testset "$file" verbose=VERBOSE for file in TEST_FILES
        Base.include(Module(), joinpath(@__DIR__, file))
    end

    # GPU tests require CUDA.jl and a functional GPU, see gpu/runtests.jl
    if get(ENV, "ADVECTRA_TEST_GPU", "false") == "true"
        include(joinpath(@__DIR__, "gpu", "runtests.jl"))
    end
end
