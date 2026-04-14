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

# ## TODO update the plots below
# #norm(sol.u[end]-HeatEquationAnalyticalSolution2(u0, domain, parameters, last(sol.t)))
#
# ## Time convergence test
# timesteps = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5]
# _,
# convergence1 = test_timestep_convergence(prob, HeatEquationAnalyticalSolution2,
#                                          timesteps, MSS1())
# _,
# convergence2 = test_timestep_convergence(prob, HeatEquationAnalyticalSolution2,
#                                          timesteps, MSS2())
# _,
# convergence3 = test_timestep_convergence(prob, HeatEquationAnalyticalSolution2,
#                                          timesteps, MSS3())
# plot(timesteps, convergence1; xaxis=:log, yaxis=:log, label="MSS1")
# plot!(timesteps, convergence2; xaxis=:log, yaxis=:log, label="MSS2", color="dark green")
# plot!(timesteps, convergence3; xaxis=:log, yaxis=:log, label="MSS3", color="orange")
# plot!(timesteps, 0.5 * timesteps .^ 3; linestyle=:dash, label=L"\frac{1}{2}dt^3")
# plot!(timesteps, 0.5 * timesteps .^ 2; linestyle=:dash, label=L"\frac{1}{2}dt^2",
#       xlabel="dt",
#       ylabel=L"||U-u_a||", title="Timestep convergence, Lin-Diffusion (N =$(domain.Nx))",
#       xticks=timesteps)
# savefig("figures/Timestep convergence, Lin-Diffusion (N =$(domain.Nx)).pdf")
#
# using JLD
# jldopen("output/linear diffusion test.jld", "w") do file
#     g = create_group(file, "data")
#     g["convergence1"] = convergence1
#     g["convergence2"] = convergence2
#     g["convergence3"] = convergence3
#     g["timesteps"] = timesteps
#     #g["colors"] = "#".*hex.(getindex.(p.series_list[1:end], :seriescolor))
# end
#
# ## Resolution convergence test
# resolutions = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192]
# _,
# convergence1 = test_resolution_convergence(prob, gaussian,
#                                            HeatEquationAnalyticalSolution2, resolutions,
#                                            MSS1())
# _,
# convergence2 = test_resolution_convergence(prob, gaussian,
#                                            HeatEquationAnalyticalSolution2, resolutions,
#                                            MSS2())
# _,
# convergence3 = test_resolution_convergence(prob, gaussian,
#                                            HeatEquationAnalyticalSolution2, resolutions,
#                                            MSS3())
#
# plot(resolutions, convergence1; xaxis=:log2, yaxis=:log, label="MSS1")
# plot!(resolutions, convergence2; xaxis=:log2, yaxis=:log, label="MSS2", color="dark green")
# plot!(resolutions, convergence3; xaxis=:log2, yaxis=:log, label="MSS3", color="orange")
# plot!(resolutions[1:(end-4)], 0.5 * exp.(-0.5 * resolutions)[1:(end-4)];
#       label=L"\frac{1}{2}\exp\left(-\frac{N}{2}\right)", linestyle=:dash,
#       xaxis=:log2, yaxis=:log, xticks=resolutions, xlabel=L"N_x \wedge N_y",
#       ylabel=L"||U-u_a||/N_xN_y",
#       title="Resolution convergence, Lin-Diffusion (dt=$(prob.dt))")
# savefig("figures/Resolution convergence, Lin-Diffusion (dt=$(prob.dt)).pdf")
#
# jldopen("output/linear diffusion resolution test.jld", "w") do file
#     g = create_group(file, "data")
#     g["convergence1"] = convergence1
#     g["convergence2"] = convergence2
#     g["convergence3"] = convergence3
#     g["resolutions"] = resolutions
#     #g["colors"] = "#".*hex.(getindex.(p.series_list[1:end], :seriescolor))
# end
