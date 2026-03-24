# ------------------------------------------------------------------------------------------
#                                 Spectral Diagnostic Tests                                 
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

modes = build_diagnostic(Val(:get_modes); axis=:kx)
modes(ic_hat, prob, 0.0)
modes = build_diagnostic(Val(:get_modes); axis=:ky)
modes(ic_hat, prob, 0.0)
modes = build_diagnostic(Val(:get_modes); axis=:both)
modes(ic_hat, prob, 0.0)
modes = build_diagnostic(Val(:get_modes); axis=:diag)
modes(ic_hat, prob, 0.0)

# Log(|modes|)
modes = build_diagnostic(Val(:get_log_modes); axis=:ky)
modes(ic_hat, prob, 0.0)
modes = build_diagnostic(Val(:get_modes); axis=:all)
modes(ic_hat, prob, 0.0)

modes = build_diagnostic(Val(:get_modes); axis=:diag)
modes(ic_hat[:, 1, 1], prob, 0.0)

psd = build_diagnostic(Val(:potential_energy_spectrum); spectrum=:radial)
k, S = psd(ic_hat, prob, 0.0)
psd = build_diagnostic(Val(:potential_energy_spectrum); spectrum=:poloidal)
k, S = psd(ic_hat, prob, 0.0)
psd = build_diagnostic(Val(:potential_energy_spectrum); spectrum=:wavenumber)
k, S = psd(ic_hat, prob, 0.0)
psd = build_diagnostic(Val(:kinetic_energy_spectrum); spectrum=:radial)
k, S = psd(ic_hat, prob, 0.0)
psd = build_diagnostic(Val(:kinetic_energy_spectrum); spectrum=:poloidal)
k, S = psd(ic_hat, prob, 0.0)
psd = build_diagnostic(Val(:kinetic_energy_spectrum); spectrum=:wavenumber)
k, S = psd(ic_hat, prob, 0.0)

psd = build_diagnostic(Val(:kinetic_energy_spectrum); spectrum=:wavenumbers)

"""
* Test that get_modes works for axis=:kx, :ky, :both and :diag for state with one field and
for state with multiple fields.
* Test that log_get_modes is computed correctly and perhaps also works for all cases 
* Test that axis=:<something> throws an error
* Test that axis=:diag on 1D domain throws error

* Test that kinetic and potential energy spectrum are computed with spectrum=:radial,
:poloidal and :wavenumber
* Test that other spectrum symbols throws an error
"""