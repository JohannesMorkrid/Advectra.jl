# ------------------------------------------------------------------------------------------
#                                   CFL Diagnostic Tests                                    
# ------------------------------------------------------------------------------------------

using Advectra
using CUDA
import Advectra: build_diagnostic, build_operator

# Minimal construction
domain = Domain(256, 256; MemoryType=CuArray)
ic = initial_condition(random_crossphased, domain) |> Advectra.memory_type(domain)
dt = 0.0001

# Emulates SpectralODEProblem
prob = (; domain=domain,
        operators=(; diff_x=build_operator(Val(:diff_x), domain),
                   diff_y=build_operator(Val(:diff_y), domain),
                   solve_phi=build_operator(Val(:solve_phi), domain)),
        dt=dt)

ic_hat = spectral_transform(ic, get_fwd(domain))

velocities = Advectra.compute_velocity(ic_hat, prob, 0.0)
velocities = Advectra.compute_velocity(ic, prob, 0.0; velocity=:burger)

CFLs = Advectra.compute_cfl(velocities, prob, 0.0, Val(:x))
CFLs = Advectra.compute_cfl(velocities, prob, 0.0, Val(:y))
CFLs = Advectra.compute_cfl(velocities, prob, 0.0, Val(:both))
CFLs = Advectra.compute_cfl(velocities, prob, 0.0, Val(:magnitude))
CFLs = Advectra.compute_cfl(velocities, prob, 0.0, Val(:something))

Advectra.cfl(ic_hat, prob, 0.00023, Val(:ExB), Val(:x));
Advectra.cfl(ic_hat, prob, 0.00023, Val(:ExB), Val(:y));
Advectra.cfl(ic_hat, prob, 0.00023, Val(:ExB), Val(:both); silent=false);
Advectra.cfl(ic_hat, prob, 0.00023, Val(:ExB), Val(:magnitude); silent=true)

cfl = build_diagnostic(Val(:cfl); velocity=:burger, component=:magnitude, silent=false)
cfl(ic, prob, 0.0)
cfl.metadata
cfl = build_diagnostic(Val(:cfl); silent=true)
cfl(ic_hat, prob, 0.0)
cfl.metadata

"""
* Test that each component mode works and that can use both ExB and burger velocity.
* Test that computes correctly for some test case 
* Test that `silent` flag works
"""