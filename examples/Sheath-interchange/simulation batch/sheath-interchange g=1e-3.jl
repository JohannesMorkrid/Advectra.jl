## Run all (alt+enter)
using Advectra
using CUDA

domain = Domain(128, 128; Lx=100, Ly=100, MemoryType=CuArray)
ic = initial_condition(random_crossphased, domain; value=1e-3)

# Linear operator
function Linear!(du, u, operators, p, t)
    @unpack laplacian = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ν, μ = p
    dη .= ν .* laplacian(η)
    dΩ .= μ .* laplacian(Ω)
end

# Non-linear operator, linearized
function NonLinear!(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack g, σ, κ = p
    ϕ = solve_phi(Ω)

    dη .= poisson_bracket(η, ϕ) - (κ - g) * diff_y(ϕ) - g * diff_y(η) + σ * ϕ
    dΩ .= poisson_bracket(Ω, ϕ) - g * diff_y(η) + σ * ϕ
end

# Parameters
parameters = (ν=1e-2, μ=1e-2, g=1e-3, σ=1e-3, κ=sqrt(1e-1))

# Time interval
tspan = [0, 5_000_000]

# Diagnostics
diagnostics = @diagnostics [
    progress(; stride=1000),
    probe_all(; positions=[(x, 0) for x in LinRange(-50, 40, 10)], stride=100),
    plot_density(; stride=5000),
    radial_flux(; stride=50),
    kinetic_energy_integral(; stride=50),
    potential_energy_integral(; stride=50),
    enstropy_energy_integral(; stride=50),
    get_log_modes(; stride=50, axis=:diag),
    cfl(; stride=500)
    #potential_energy_spectrum(; spectrum=:radial, stride=500),
    #potential_energy_spectrum(; spectrum=:poloidal, stride=500),
    #kinetic_energy_spectrum(; spectrum=:radial, stride=500),
    #kinetic_energy_spectrum(; spectrum=:poloidal, stride=500)
]

# Collection of specifications defining the problem to be solved
prob = SpectralODEProblem(Linear!, NonLinear!, ic, domain, tspan; p=parameters, dt=1e-1,
                          operators=:all, diagnostics=diagnostics)

# Output
output_file_name = joinpath(@__DIR__, "../output", "sheath-interchange long time series.h5")
output = Output(prob; filename=output_file_name, simulation_name=:parameters,
                storage_limit="1 GB", store_locally=false, resume=true)

## Solve and plot
sol = spectral_solve(prob, MSS3(), output;)

using SMTPClient
send_mail("g=1e-3 finnished, go analyse the data and see if it has different PDF!")
close(output)
