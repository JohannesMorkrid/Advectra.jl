## Run all (alt+enter)
using Advectra
using CUDA

domain = Domain(1024, 1024; Lx=50, Ly=50, MemoryType=CuArray, precision=Float32)
n0 = initial_condition(log_gaussian, domain; A=1, B=1, l=1)
ic = cat(n0, zero(n0); dims=3)

# Linear operator
function Linear!(du, u, operators, p, t)
    @unpack laplacian = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack κ, ν = p

    dη .= κ * laplacian(η)
    dΩ .= ν * laplacian(Ω)
end

# Non-linear operator
function NonLinear!(du, u, operators, p, t)
    @unpack laplacian, solve_phi, poisson_bracket, quadratic_term, diff_x,
    diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack κ, ν = p

    ϕ = solve_phi(Ω)
    dη .= poisson_bracket(η, ϕ) + κ * quadratic_term(diff_x(η), diff_x(η)) +
          κ * quadratic_term(diff_y(η), diff_y(η))
    dΩ .= poisson_bracket(Ω, ϕ) - diff_y(η)
end

# Parameters
parameters = (ν=1e-3, κ=1e-3)

# Time interval
tspan = [0, 16]

# Array of diagnostics want
diagnostics = @diagnostics [
    radial_COM(; stride=1),
    progress(; stride=100),
    plot_density(; stride=1000)
]

# The problem
prob = SpectralODEProblem(Linear!, NonLinear!, ic, domain, tspan; p=parameters, dt=1e-3,
                          diagnostics=diagnostics, operators=:all)

# Inverse transform
inverse_transformation!(u) = @. u[:, :, 1] = exp(u[:, :, 1]) - 1

# The output
output_file_name = joinpath(@__DIR__, "output", "Kube 2011 PoP.h5")
output = Output(prob; filename=output_file_name, simulation_name=:parameters,
                physical_transform=(inverse_transformation!), storage_limit="0.5 GB",
                store_locally=false, resume=false)

# Solve
sol = spectral_solve(prob, MSS3(), output;)

using Plots
plot(output.simulation["Radial COM/t"][1:(end-1)],
     output.simulation["Radial COM/data"][2, 1:(end-1)])
