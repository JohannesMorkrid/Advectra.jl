# ------------------------------------------------------------------------------------------
#                                Progressbar Diagnostic Test                                
# ------------------------------------------------------------------------------------------

using Advectra
domain = Domain(256, 256)
prob = (; domain=domain)
tspan = [0.0, 1.0]
dt = 1e-3

ic = initial_condition(isolated_blob, domain)
ic_hat = spectral_transform(ic, get_fwd(domain))

import Advectra: build_diagnostic
progress = build_diagnostic(Val(:progress); tspan=tspan, dt=dt)
for i in 0.0:dt:1.0
    progress(ic_hat, prob, i)
    sleep(0.005)
end

"""
Test that no error is thrown.
"""