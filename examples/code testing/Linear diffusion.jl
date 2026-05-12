## Run all (alt+enter)
using Advectra

## Run scheme test
domain = Domain(32, 32; Lx=50, Ly=50)
u0 = initial_condition(gaussian, domain)

# Diffusion 
Linear(du, u, operators, p, t) = du .= p.ν * operators.laplacian(u)
NonLinear(du, u, operators, p, t) = du .= zero(u)

# Parameters
parameters = (ν=0.5,)

tspan = [0.0, 2.0]

diagnostics = @diagnostics [progress(; stride=10), sample_density(; stride=10)]

prob = SpectralODEProblem(Linear, NonLinear, u0, domain, tspan; p=parameters, dt=1e-3,
                          diagnostics=diagnostics, operators=:all)

output_file_name = joinpath(@__DIR__, "output", "linear diffusion.h5")
output = Output(prob; filename=output_file_name, simulation_name=:parameters)

## Solve and plot
sol = spectral_solve(prob, MSS3(), output)

# println(sol.simulation["Density/data"][end, :, :, 1])
println(sol.simulation["Density/data"][:, :, end])
