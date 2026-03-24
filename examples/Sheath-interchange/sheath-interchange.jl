## Run all (alt+enter)
using Advectra
using CUDA

domain = Domain(256, 256; Lx=48, Ly=48, MemoryType=CuArray, precision=Float64)
ic = initial_condition(random_crossphased, domain; value=1e-3)

# Linear operator
function Linear(du, u, operators, p, t)
    @unpack laplacian = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ν, μ = p
    dη .= ν * laplacian(η)
    dΩ .= μ * laplacian(Ω)
end

# Non-linear operator, linearized potential and diffusion
function NonLinear(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, diff_y, quadratic_term = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack g, σ = p
    ϕ = solve_phi(Ω)

    dη .= poisson_bracket(η, ϕ) - (1 - g) * diff_y(ϕ) - g * diff_y(η) + σ * ϕ
    dΩ .= poisson_bracket(Ω, ϕ) - g * diff_y(η) + σ * ϕ
end

# Parameters
parameters = (κ=sqrt(1e1), g=1e-1, σ=1e-1, ν=1e-2, μ=1e-2)

# Time interval
tspan = [0.0, 5000.0]

# Diagnostics
diagnostics = @diagnostics [
    progress(; stride=1000),
    probe_all(; positions=[(x, 0) for x in LinRange(-10, 10, 10)], stride=100),
    plot_density(; stride=1000),
    radial_flux(; stride=50),
    kinetic_energy_integral(; stride=50),
    potential_energy_integral(; stride=50),
    enstropy_energy_integral(; stride=50),
    get_log_modes(; stride=50, axis=:diag),
    cfl(; stride=5000),
    sample_density(; storage_limit="1 GB"),
    sample_vorticity(; storage_limit="1 GB")
    #sample_potential(; storage_limit="1 GB")
    #potential_energy_spectrum(; spectrum=:radial, stride=50),
    #potential_energy_spectrum(; spectrum=:poloidal, stride=50),
    #kinetic_energy_spectrum(; spectrum=:radial, stride=50)
    #kinetic_energy_spectrum(; spectrum=:poloidal, stride=50)
]

# Collection of specifications defining the problem to be solved
prob = SpectralODEProblem(Linear, NonLinear, ic, domain, tspan; p=parameters, dt=1e-3,
                          operators=:all, diagnostics=diagnostics)

# Inverse transform
inverse_transformation!(u) = @. u[:, :, 1] = exp(1e-2 * u[:, :, 1]) - 1

# Output
output_file_name = joinpath(@__DIR__, "output", "sheath-interchange_example.h5")
output = Output(prob; filename=output_file_name, simulation_name=:parameters,
                #physical_transform=inverse_transformation!,
                storage_limit="1 GB",
                store_locally=false, resume=true)

## Solve and plot
sol = spectral_solve(prob, MSS3(), output;)

using SMTPClient
send_mail("Simulation finnished!"; attachment="benkadda.gif")
close(output)
