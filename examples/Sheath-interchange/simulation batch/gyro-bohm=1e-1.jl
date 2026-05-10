## Run all (alt+enter)
using Advectra
using CUDA

domain = Domain(256, 256; Lx=48, Ly=48, MemoryType=CuArray, precision=Float64)
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
    @unpack g, σ = p
    ϕ = solve_phi(Ω)

    dη .= poisson_bracket(η, ϕ) - (1 - g) * diff_y(ϕ) - g * diff_y(η) + σ * ϕ
    dΩ .= poisson_bracket(Ω, ϕ) - g * diff_y(η) + σ * ϕ
end

# Parameters
parameters = (ν=1e-2, μ=1e-2, g=1e-1, σ=1e-1)

# Time interval
tspan = [0, 500_000]

# Diagnostics
diagnostics = @diagnostics [
    progress(; stride=9711),
    probe_all(; positions=[(x, 0) for x in range(-24, 19.6, 10)], stride=10),
    plot_density(; stride=500),
    radial_flux(; stride=100),
    kinetic_energy_integral(; stride=100),
    potential_energy_integral(; stride=100),
    enstropy_energy_integral(; stride=100),
    get_log_modes(; stride=500, axis=:diag),
    cfl(; stride=500, silent=true)
    #potential_energy_spectrum(; spectrum=:radial, stride=50),
    #potential_energy_spectrum(; spectrum=:poloidal, stride=50),
    #kinetic_energy_spectrum(; spectrum=:radial, stride=50),
    #kinetic_energy_spectrum(; spectrum=:poloidal, stride=50)
]

# Collection of specifications defining the problem to be solved
prob = SpectralODEProblem(Linear!, NonLinear!, ic, domain, tspan; p=parameters, dt=2e-3,
                          diagnostics=diagnostics)

# Output
output_file_name = joinpath(@__DIR__, "../output", "gyro-bohm=1e-1.h5")
output = Output(prob; filename=output_file_name, simulation_name=:parameters,
                storage_limit="1.2 GB", store_locally=false, resume=true)

## Solve and plot
sol = spectral_solve(prob, MSS3(), output;)

using SMTPClient
send_mail("sigma=1e-1 finnished, go analyse the data!")
close(output)
