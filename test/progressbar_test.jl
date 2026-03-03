# ------------------------------------------------------------------------------------------
#                                Progressbar Diagnostic Test                                
# ------------------------------------------------------------------------------------------

using Test
using HasegawaWakatani

@testset "Diagnostics and Progress" begin
    nx, ny = 256, 256
    domain = Domain(nx, ny)
    prob = (; domain=domain)
    tspan = (0.0, 1.0) # Use a tuple for standard tspan convention
    dt = 1e-3

    ic = initial_condition(isolated_blob, domain)
    fwd_plan = get_fwd(domain)
    ic_hat = spectral_transform(ic, fwd_plan)

    test_range = 0.0:dt:0.05 

    @testset "Progress Bar Diagnostic" begin
        import HasegawaWakatani: build_diagnostic
        
        progress = build_diagnostic(Val(:progress); tspan=tspan, dt=dt)
        
        for t in test_range
            progress(ic_hat, prob, t)
        end
        @test true # Dummy test to keep the set valid
    end
end
