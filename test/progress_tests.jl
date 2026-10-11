# ------------------------------------------------------------------------------------------
#                                 Progress Diagnostic Tests
# ------------------------------------------------------------------------------------------

using Test
using Advectra
import Advectra: build_diagnostic

@testset "Progress diagnostic" begin
    domain = Domain(16)
    prob = (; domain=domain)
    tspan = (0.0, 1.0)
    dt = 1e-1

    state_hat = spectral_transform(initial_condition(isolated_blob, domain), fwd_plan(domain))

    progress = build_diagnostic(Val(:progress); tspan=tspan, dt=dt)
    @test progress.name == "Progress"
    @test !progress.stores_data
    @test progress.assumes_spectral_state

    # The progress bar counts the number of steps taken since tspan[1]
    progress_bar = only(progress.args)
    @test progress_bar.n == 10
    for step in 0:10
        progress(state_hat, prob, first(tspan) + step * dt)
        @test progress_bar.counter == step
    end
end
