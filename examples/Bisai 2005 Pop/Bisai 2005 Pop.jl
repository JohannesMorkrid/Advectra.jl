## Run all (alt+enter)
using Advectra
using CUDA

## Run scheme test for Burgers equation
domain = Domain(128, 128; Lx=160, Ly=160, MemoryType=CuArray, precision=Float32)
ic = initial_condition(random_crossphased, domain; value=1e-3)
ic[:, :, 1] .+= 0.5

# Linear operator
function Linear(du, u, operators, p, t)
    @unpack laplacian = operators
    n, Ω = eachslice(u; dims=3)
    dn, dΩ = eachslice(du; dims=3)
    @unpack ν, μ = p

    dn .= ν * laplacian(n)
    dΩ .= μ * laplacian(Ω)
end

# Current implementation of source
source(x, y, S_0, λ_s) = @. S_0 * exp(-((x + 80 - λ_s) / λ_s)^2) + 0 * y

S = get_fwd(domain) * CuArray(source(domain.x', domain.y, 5e-4, 5))
using Plots
heatmap(Array(get_bwd(domain) * S); aspect_ratio=:equal)

# Non-linear operator, fully non-linear
function NonLinear(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, diff_y, quadratic_term = operators
    @unpack spectral_exp, spectral_expm1, spectral_log = operators
    @unpack g, σ_0 = p
    n, Ω = eachslice(u; dims=3)
    dn, dΩ = eachslice(du; dims=3)
    ϕ = solve_phi(Ω)

    dn .= poisson_bracket(n, ϕ) - g * diff_y(n) + g * quadratic_term(n, diff_y(ϕ)) -
          σ_0 * quadratic_term(n, spectral_exp(-ϕ)) .+ S#.+ source()
    dΩ .= poisson_bracket(Ω, ϕ) - g * diff_y(spectral_log(n)) - σ_0 * spectral_expm1(-ϕ)
end

# Parameters
parameters = (μ=1e-2, ν=1e-2, g=8e-4, σ_0=2e-4, S_0=5e-4, λ_s=5.0)

tspan = [0.0, 50_000_000.0]

# Diagnostics
diagnostics = @diagnostics [
    progress(; stride=1000),
    #probe_all(; positions=[(x, 0) for x in LinRange(-80, 64, 10)], stride=10),
    plot_density(; stride=50)
    #radial_flux(; stride=50),
    #kinetic_energy_integral(; stride=50),
    #potential_energy_integral(; stride=50),
    #enstropy_energy_integral(; stride=50),
    #get_log_modes(; stride=50, axis=:diag),
    #cfl(; stride=50)
    #potential_energy_spectrum(; spectrum=:radial, stride=50),
    #potential_energy_spectrum(; spectrum=:poloidal, stride=50),
    #kinetic_energy_spectrum(; spectrum=:radial, stride=50),
    #kinetic_energy_spectrum(; spectrum=:poloidal, stride=50)
]

prob = SpectralODEProblem(Linear, NonLinear, ic, domain, tspan; p=parameters, dt=1,
                          operators=:all, diagnostics=diagnostics)

# Output
output_file_name = joinpath(@__DIR__, "output", "Bisai.h5")
output = Output(prob; filename=output_file_name, simulation_name=:parameters,
                store_locally=false, resume=true)

# Solve and plot
sol = spectral_solve(prob, MSS3(), output;)

# bisai_php_19_052509 gives Lₓ = L_y = 160ρₛ, g = 1e-3, 128-modes
# bisai_php_12_072520 gives dt = \omega_c^{-1}, g = 8e-4
