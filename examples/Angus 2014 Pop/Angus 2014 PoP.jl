## Run all (alt+enter)
using Advectra
using CUDA

domain = Domain(256, 256; Lx=50, Ly=50, MemoryType=CuArray, precision=Float64)

# Blob amplitude
A = 2

# Check documentation to see other initial conditions
ic = initial_condition(isolated_blob, domain; A=A, B=1)

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
    @unpack diff_y, poisson_bracket, solve_phi, grad_dot_grad = operators
    ϕ = solve_phi(θ, Ω)
    # Vorticity differential
    grad_dot_grad(dΩ, ϕ, ϕ)
    poisson_bracket(dθ, θ, dΩ)
    dθ .*= 0.5
    poisson_bracket(dΩ, Ω, ϕ)
    dΩ .+= dθ
    diff_y(dθ, θ)
    dΩ .-= dθ
    # Density differential
    poisson_bracket(dθ, θ, ϕ)
end

# Parameters
parameters = (ν=1e-2, κ=1e-2, A=A)

# Time interval
tspan = [0.0, 20.0]

# Array of diagnostics want
diagnostics = @diagnostics [
    probe_density(; positions=[(5, 0), (8.5, 0), (11.25, 0), (14.375, 0)], stride=10),
    radial_COM(; stride=1),
    progress(; stride=-1),
    cfl(; stride=250, silent=true, storage_limit="2KB"),
    plot_vorticity(; stride=1000),
    plot_potential(; stride=1000),
    plot_density(; stride=1000)
]

# Collection of specifications defining the problem to be solved
prob = SpectralODEProblem(Linear, NonLinear, ic, domain, tspan; p=parameters, dt=2.5e-3,
                          boussinesq=false, diagnostics=diagnostics, density=:linear,
                          operators=:all)

inverse_transformation!(u) = @. u[:, :, 1] = u[:, :, 1] - 1

# The output
output = Output(prob; filename="Angus 2014 PoP.h5", simulation_name=:parameters,
                physical_transform=inverse_transformation!, storage_limit="0.5 GB",
                store_locally=false, resume=false)

# Solve and plot
sol = spectral_solve(prob, MSS3(), output; debug=true)