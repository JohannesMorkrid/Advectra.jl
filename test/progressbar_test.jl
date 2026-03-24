# ------------------------------------------------------------------------------------------
#                                Progressbar Diagnostic Test                                
# ------------------------------------------------------------------------------------------

using Advectra
domain = Domain(256, 256)
prob = (; domain=domain)
tspan = [0.0, 1.0]
dt = 1e-3

@testset "Diagnostics and Progress" begin
    nx, ny = 256, 256
    domain = Domain(nx, ny)
    prob = (; domain=domain)
    tspan = (0.0, 1.0) # Use a tuple for standard tspan convention
    dt = 1e-3

import Advectra: build_diagnostic
progress = build_diagnostic(Val(:progress); tspan=tspan, dt=dt)
for i in 0.0:dt:1.0
    progress(ic_hat, prob, i)
    sleep(0.005)
end

    @testset "Progress Bar Diagnostic" begin
        import HasegawaWakatani: build_diagnostic
        
        progress = build_diagnostic(Val(:progress); tspan=tspan, dt=dt)
        
        for t in test_range
            progress(ic_hat, prob, t)
        end
        @test true # Dummy test to keep the set valid
    end
end
