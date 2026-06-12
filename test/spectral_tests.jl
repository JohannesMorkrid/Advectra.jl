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

using Advectra
rdomain = Domain(257; L=48, real_transform=true)
domain = Domain(257; L=48, real_transform=false)

A = randn(257, 257)
Ar_hat = get_fwd(rdomain)*A
A_hat = get_fwd(domain)*A
f = abs2 #identity
dims = 2

using Statistics
A_zonal = mean(A; dims=1)
A_streamer = mean(A; dims=2)

A_zonal_hat = selectdim(A_hat, 1, 1:1)
A_zonal_hat == selectdim(Ar_hat, 1, 1:1)

A_streamer_hat = selectdim(A_hat, 2, 1:1)
A_streamer_hat == selectdim(Ar_hat, 2, 1:1)

# Using normal transform

# Build operators
diff_x = build_operator(:diff_x, domain)
diff_y = build_operator(:diff_y, domain)

A_zonal = mean(A; dims=1) .+ 0.0*domain.y
A_streamer = mean(A; dims=2) .+ 0.0*domain.x'

A_zonal_hat = get_fwd(domain)*A_zonal
A_streamer_hat = get_fwd(domain)*A_streamer

#A_zonal_hat = selectdim(A_hat, 1, 1:1)
#A_streamer_hat = selectdim(A_hat, 2, 1:1)

vx_hat = -diff_y(A_streamer_hat)
vy_hat = diff_x(A_zonal_hat)

vx = get_bwd(domain)*vx_hat
vy = get_bwd(domain)*vy_hat

mean(abs2, vx)
mean(abs2, vy)

A_hat = get_fwd(domain)*A
A_zonal_hat = selectdim(A_hat, 1, 1:1)
A_streamer_hat = selectdim(A_hat, 2, 1:1)

vx_hat = -diff_y(A_streamer_hat)
vy_hat = diff_x(A_zonal_hat)

spectral_sum(abs2, vx_hat, domain)/length(domain)^2 ≈ mean(abs2, vx)
spectral_sum(abs2, vy_hat, domain)/length(domain)^2 ≈ mean(abs2, vy)

# Using real transform

rdiff_x = build_operator(:diff_x, rdomain)
rdiff_y = build_operator(:diff_y, rdomain)

#Ar_zonal = mean(A; dims=1) .+ 0.0*rdomain.y
#Ar_streamer = mean(A; dims=2) .+ 0.0*rdomain.x'

#Ar_zonal_hat = get_fwd(rdomain)*Ar_zonal
#Ar_streamer_hat = get_fwd(rdomain)*Ar_streamer

Ar_hat = get_fwd(rdomain)*A
Ar_zonal_hat = selectdim(Ar_hat, 1, 1:1)
Ar_streamer_hat = selectdim(Ar_hat, 2, 1:1)

rvx_hat = -rdiff_y(Ar_streamer_hat)
rvy_hat = rdiff_x(Ar_zonal_hat)

spectral_sum(abs2, rvx_hat, rdomain)/length(domain)^2 ≈ mean(abs2, vx)
spectral_sum(abs2, rvy_hat, rdomain)/length(domain)^2 ≈ mean(abs2, vy)

# Tests
spectral_sum(f, A, domain)
spectral_sum(spectral_sum(f, A, domain; dims=1), domain)
spectral_sum(spectral_sum(f, A, domain; dims=2), domain)
spectral_sum(spectral_sum(f, A, domain; dims=2), domain; dims=1)[1]
spectral_sum(spectral_sum(f, A, domain; dims=1), domain; dims=2)[1]
spectral_sum(spectral_sum(f, A, domain; dims=1), domain; dims=1)
spectral_sum(spectral_sum(f, A, domain; dims=2), domain; dims=2)

spectral_sum(f, selectdim(A, 1, 1), domain)
spectral_sum(f, selectdim(A, 1, 1:1), domain)
spectral_sum(f, selectdim(A, 2, 1:1), domain)
spectral_sum(f, selectdim(A, 2, 1), domain)