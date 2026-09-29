# ------------------------------------------------------------------------------------------
#                                  Spectrum Diagnostic Tests
# ------------------------------------------------------------------------------------------

using Test
using Advectra
import Advectra: build_operator

@testset "Spectrum default arguments" begin
    domain = Domain(32; L=10)
    state_hat = spectral_transform(initial_condition(isolated_blob, domain), fwd_plan(domain))
    prob = (; domain=domain,
            operators=(; diff_x=build_operator(Val(:diff_x), domain),
                       diff_y=build_operator(Val(:diff_y), domain),
                       solve_phi=build_operator(Val(:solve_phi), domain)))

    # Calling without a spectrum must equal calling with the documented default
    @testset "$method" for (method, default) in [(potential_energy_spectrum, :radial),
                                                 (kinetic_energy_spectrum, :radial),
                                                 (flux_spectrum, :poloidal),
                                                 (enstrophy_spectrum, :radial),
                                                 (electrostatic_potential_spectrum,
                                                  :radial)]
        expected = method(state_hat, prob, 0.0, Val(default))
        @test method(state_hat, prob, 0.0) == expected
        @test method(state_hat, prob, 0.0; spectrum=default) == expected
    end
end
