## Run all (alt+enter)
using Advectra
using CUDA

domain = Domain(256, 256; Lx=50, Ly=50, MemoryType=CuArray, precision=Float32)
# Check documentation to see other initial conditions
ic = initial_condition(isolated_blob, domain)

# Linear operator
function Linear(du, u, operators, p, t)
    @unpack κ, ν = p
    θ, Ω = eachslice(u; dims=3)
    dθ, dΩ = eachslice(du; dims=3)
    @unpack laplacian = operators
    dθ .= κ .* laplacian(θ)
    dΩ .= ν .* laplacian(Ω)
end

# Non-linear operator
function NonLinear(du, u, operators, p, t)
    θ, Ω = eachslice(u; dims=3)
    dθ, dΩ = eachslice(du; dims=3)
    @unpack diff_y, poisson_bracket, solve_phi = operators
    ϕ = solve_phi(Ω)
    poisson_bracket(dΩ, Ω, ϕ)
    diff_y(dθ, θ)
    dΩ .-= dθ
    poisson_bracket(dθ, θ, ϕ)
end

# Parameters
parameters = (ν=1e-2, κ=1e-2)

# Time interval
tspan = [0.0, 20.0]

# Array of diagnostics want
diagnostics = @diagnostics [
    probe_density(; positions=[(5, 0), (8.5, 0), (11.25, 0), (14.375, 0)], stride=10),
    radial_COM(; stride=1),
    progress(; stride=-1),
    # cfl(; stride=250, silent=true, storage_limit="2KB"),
    sample_density(stride=10),
    plot_vorticity(; stride=1000),
    plot_potential(; stride=1000),
    plot_density(; stride=1000)
]

# Collection of specifications defining the problem to be solved
prob = SpectralODEProblem(Linear, NonLinear, ic, domain, tspan; p=parameters, dt=2.5e-3,
                          boussinesq=true, diagnostics=diagnostics)

# The output
output = Output(prob; filename="Garcia 2005 PoP.h5", simulation_name=:parameters,
                storage_limit="5 GB", store_locally=false, resume=false)

# Solve and plot
sol = spectral_solve(prob, MSS3(), output;)

using Plots
using LaTeXStrings
plot(sol.simulation["Density probe/t"][:], sol.simulation["Density probe/data"][:, :]';
     xlabel=L"$t$",
     ylabel=L"$\theta$", label=["(5.0, 0.0)" "(8.5, 0.0)" "(11.25, 0.0)" "(14.375, 0.0)"])
