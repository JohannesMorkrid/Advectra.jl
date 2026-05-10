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
prob = (; domain=domain, dt=dt)

COM_diagnostic = build_diagnostic(Val(:radial_COM))
COM_diagnostic(ic, prob, 0.0)
ic = initial_condition(random_crossphased, domain) |> Advectra.memory_type(domain)
COM_diagnostic(ic, prob, 1.0)

"""
* Test that COM diagnostic work.
* Test that velocity is computed when giving two different arrays at different times.
"""