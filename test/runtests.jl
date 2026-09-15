using Advectra
using Test

@testset "Advectra" begin
    # Write your tests here.
    include("domain_tests.jl")
    include("display_tests.jl")
    include("integration_tests.jl")
    include("progressbar_test.jl")
    include("operator_tests.jl")
    include("scheme_tests.jl")
end
