## Run all (alt+enter)
using Advectra

## Run test for Burgers equation
domain = Domain(1, 1024; Lx=1, Ly=20)
u0 = initial_condition(gaussianWallY, domain)

# Diffusion 
function Linear(u, operators, p, t)
    @unpack laplacian = operators
    @unpack ν = p
    ν * laplacian(u)
end

# Burgers equation 
function NonLinear(u, operators, p, t)
    @unpack quadratic_term, diff_y = operators
    return -quadratic_term(u, diff_y(u))
end

# Parameters
parameters = (ν=0.0,)

# Break down time 
diff_y = Advectra.build_operator(Val(:diff_y), domain)
dudy = diff_y(get_fwd(domain) * u0)
t_b = -1 / (minimum(real(get_bwd(domain) * dudy)))

# Time span
tspan = [0, 0.8 * t_b]

diagnostics = @diagnostics [
    progress(; stride=10),
    cfl(; stride=10, velocity=:burger),
    sample_density(; stride=1000)
]

# Initialize problem
prob = SpectralODEProblem(Linear, NonLinear, u0, domain, tspan; p=parameters, dt=0.0001,
                          diagnostics=diagnostics,
                          additional_operators=[OperatorRecipe(:quadratic_term)])

# Initialize output
output_file_name = joinpath(@__DIR__, "output", "burgers equation gaussian.h5")
output = Output(prob; filename=output_file_name, simulation_name=:parameters)

## Solve problem
sol = spectral_solve(prob, MSS3(), output)

## TODO update plots below
using Plots, LaTeXStrings
plot(domain.y, sol.simulation["Density/data"][:, 1, end-1];
     label=L"U(" * "\$$(round(last(tspan),digits=2))\$" * L")")
plot!(domain.y, burgers_equation_analytical_solution(u0, domain, parameters, last(t_span));
      linestyle=:dash, label=L"u_a(" * "\$$(round(last(t_span),digits=2))\$" * L")",
      c=:yellow)
plot!(; xlabel=L"y", ylabel=L"u(y)", labelfontsize=10)
savefig("figures/burgers steepning gaussian.pdf")

using JLD
jldopen("output/burgers steepning gaussian.jld", "w") do file
    g = create_group(file, "data")
    g["analytical"] = burgers_equation_analytical_solution(u0, domain, parameters,
                                                           last(t_span))
    g["approximate"] = sol.u[end]
    g["t"] = last(t_span)
    #g["colors"] = "#".*hex.(getindex.(p.series_list[1:end], :seriescolor))
end

#plot(sol.u[end] - burgers_equation_analytical_solution(u0, domain, parameters, last(t_span)))

## Time convergence test
timesteps = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7]
_,
convergence1 = test_timestep_convergence(prob, burgers_equation_analytical_solution,
                                         timesteps, MSS1())
_,
convergence2 = test_timestep_convergence(prob, burgers_equation_analytical_solution,
                                         timesteps, MSS2())
_,
convergence3 = test_timestep_convergence(prob, burgers_equation_analytical_solution,
                                         timesteps, MSS3())
plot(timesteps, convergence1; xaxis=:log, yaxis=:log, label="MSS1")
plot!(timesteps, convergence2; xaxis=:log, yaxis=:log, label="MSS2", color="dark green")
plot!(timesteps, convergence3; xaxis=:log, yaxis=:log, label="MSS3", color="orange",
      xlabel="dt",
      ylabel=L"||U-u_a||", title="Timestep convergence, Burgers equation (N =$(domain.Ny))",
      xticks=timesteps)
savefig("figures/Timestep convergence, Burgers equation (N =$(domain.Ny)).pdf")

using JLD
jldopen("output/burgers gaussian timestep.jld", "w") do file
    g = create_group(file, "data")
    g["convergence1"] = convergence1
    g["convergence2"] = convergence2
    g["convergence3"] = convergence3
    g["timesteps"] = timesteps
    #g["colors"] = "#".*hex.(getindex.(p.series_list[1:end], :seriescolor))
end

## Resolution convergence test
resolutions = [2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
_,
convergence1 = test_resolution_convergence(prob, gaussianWallY,
                                           burgers_equation_analytical_solution,
                                           resolutions, MSS1(); oneDimensional=true)
_,
convergence2 = test_resolution_convergence(prob, gaussianWallY,
                                           burgers_equation_analytical_solution,
                                           resolutions, MSS2(); oneDimensional=true)
_,
convergence3 = test_resolution_convergence(prob, gaussianWallY,
                                           burgers_equation_analytical_solution,
                                           resolutions, MSS3(); oneDimensional=true)

plot(resolutions, convergence1; xaxis=:log2, yaxis=:log, label="MSS1")
plot!(resolutions, convergence2; xaxis=:log2, yaxis=:log, label="MSS2", color="dark green")
plot!(resolutions, convergence3; xaxis=:log2, yaxis=:log, label="MSS3", color="orange")
plot!(resolutions[1:(end-4)], 0.5 * exp.(-0.5 * resolutions)[1:(end-4)];
      label=L"\frac{1}{2}\exp\left(-\frac{N}{2}\right)", linestyle=:dash,
      xaxis=:log2, yaxis=:log, xticks=resolutions, xlabel=L"N_x \wedge N_y",
      ylabel=L"||U-u_a||/N_xN_y",
      title="Resolution convergence, Burgers equation (dt=$(prob.dt))")
savefig("figures/Resolution convergence, Burgers equation (dt=$(prob.dt)).pdf")

jldopen("output/burgers gaussian resolution.jld", "w") do file
    g = create_group(file, "data")
    g["convergence1"] = convergence1
    g["convergence2"] = convergence2
    g["convergence3"] = convergence3
    g["resolutions"] = resolutions
    #g["colors"] = "#".*hex.(getindex.(p.series_list[1:end], :seriescolor))
end

## ----------------------------------- Plot ------------------------------------------------

plot(domain.y, u0; xlabel=L"y", ylabel=L"u(y)", label="",
     title="Gaussian initial condition")
savefig("figures/Gaussian intial condition.pdf")
send_mail("Burgers test (dt=1e-4) finished!")
