# ------------------------------------------------------------------------------------------
#                                        Flux Tests                                         
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

radial_flux = build_diagnostic(Val(:radial_flux))

ic_hat = cat(get_fwd(domain) * ic[:, :, 1], get_fwd(domain) * ic[:, :, 2]; dims=3)

radial_flux(ic_hat, prob, 0.0)

poloidal_flux = build_diagnostic(Val(:poloidal_flux))
poloidal_flux(ic_hat, prob, 0.0)

"""
* Test that flux works for both GPU and CPU arrays
* Test that flux is indeed correctly computed
"""