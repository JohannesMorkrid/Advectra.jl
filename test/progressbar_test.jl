# ------------------------------------------------------------------------------------------
#                                Progressbar Diagnostic Test                                
# ------------------------------------------------------------------------------------------
using Advectra

@testset "Diagnostics and Progress" begin
    Nx, Ny = 256, 256
    domain = Domain(Nx, Ny)
    prob = (; domain=domain)
    tspan = (0.0, 1.0) # Use a tuple for standard tspan convention
    dt = 1e-1

    ic = initial_condition(isolated_blob, domain)
    ic_hat = spectral_transform(ic, get_fwd(domain))

    progress = build_diagnostic(Val(:progress); tspan=tspan, dt=dt)
    for i in 0.0:dt:1.0
        progress(ic_hat, prob, i*dt)
        sleep(0.005)
    end

    @test true # Dummy test to keep the set valid
end
