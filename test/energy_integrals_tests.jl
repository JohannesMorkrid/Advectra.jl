# ------------------------------------------------------------------------------------------
#                                  Energy Diagnostic Tests                                  
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
                   solve_phi=build_operator(Val(:solve_phi), domain),
                   laplacian=build_operator(Val(:laplacian), domain),
                   quadratic_term=build_operator(Val(:quadratic_term), domain),
                   hyper_laplacian=build_operator(Val(:laplacian), domain; order=3)),
        p=(c=0.01,),
        dt=dt)

ic_hat = spectral_transform(ic, get_fwd(domain))

# Parseval energy integrals:

kinetic_energy = build_diagnostic(Val(:kinetic_energy_integral))
kinetic_energy(ic_hat, prob, 0.0)
potential_energy = build_diagnostic(Val(:potential_energy_integral))
potential_energy(ic_hat, prob, 0.0)
total_energy = build_diagnostic(Val(:total_energy_integral))
total_energy(ic_hat, prob, 0.0)
enstropy_energy = build_diagnostic(Val(:enstropy_energy_integral))
enstropy_energy(ic_hat, prob, 0.0)

# Dissipative energy integrals:

resistive_dissipation = build_diagnostic(Val(:resistive_dissipation_integral);
                                         adiabaticity_symbol=:c)
resistive_dissipation(ic_hat, prob, 0.0)
potential_dissipation = build_diagnostic(Val(:potential_dissipation_integral);
                                         diffusivity_symbol=:c)
potential_dissipation(ic_hat, prob, 0.0)
kinetic_dissipation = build_diagnostic(Val(:kinetic_dissipation_integral);
                                       viscosity_symbol=:c)
kinetic_dissipation(ic_hat, prob, 0.0)
viscous_dissipation = build_diagnostic(Val(:viscous_dissipation_integral);
                                       viscosity_symbol=:c, diffusivity_symbol=:c)
viscous_dissipation(ic_hat, prob, 0.0)
enstropy_dissipation = build_diagnostic(Val(:enstropy_dissipation_integral);
                                        viscosity_symbol=:c, diffusivity_symbol=:c)
enstropy_dissipation(ic_hat, prob, 0.0)

# Energy evolution integrals:
energy_evolution = build_diagnostic(Val(:energy_evolution_integral); adiabaticity_symbol=:c,
                                    diffusivity_symbol=:c, viscosity_symbol=:c)
energy_evolution(ic_hat, prob, 0.0)
enstropy_evolution = build_diagnostic(Val(:enstropy_evolution_integral);
                                      diffusivity_symbol=:c, viscosity_symbol=:c)
enstropy_evolution(ic_hat, prob, 0.0)

import Advectra: parsevals_theorem, integral_of_quadratic_term
# Test robustness, should give 1
domain = Domain(256, 256; real_transform=false, Lx=10, Ly=10)
parsevals_theorem(ones(256, 256), domain)
domain = Domain(256, 256; real_transform=true, Lx=50, Ly=500)
parsevals_theorem(ones(129, 256), domain; compute_density=true)

domain = Domain(8, 4; MemoryType=CuArray, real_transform=true, Lx=50, Ly=500,
                dealiased=true)
prob = (; domain=domain,
        operators=(; diff_x=build_operator(Val(:diff_x), domain),
                   diff_y=build_operator(Val(:diff_y), domain),
                   solve_phi=build_operator(Val(:solve_phi), domain),
                   laplacian=build_operator(Val(:laplacian), domain),
                   quadratic_term=build_operator(Val(:quadratic_term), domain)),
        p=(c=0.01,),
        dt=dt)
A = ones(4, 8) |> Advectra.memory_type(domain)
A_hat = get_fwd(domain) * A
integral_of_quadratic_term(A_hat, A_hat, domain, prob.operators.quadratic_term)
prob.operators.quadratic_term.dealiasing_coefficient
#2.25

parsevals_theorem(ones(129, 256), domain; compute_density=false)

"""
* Test that all energy integrals work:
    - Kinetic energy
    - Potential energy
    - Total energy
    - Enstropy energy
* Test that all dissipative energy integrals work:
    - Resistive dissipation
    - Potential dissipation
    - Kinetic dissipation
    - Viscous dissipation
    - Enstropy dissipation
* Test that all energy evolution integrals work:
    - Energy evolution 
    - Enstrophy evolution
* Test that can construct using different coefficient symbols
* Test that parsevals_theorem and integral_of_quadratic_term are robust to changes in size 
and lengths and if real_transform or not.
"""