## Run all (alt+enter)
using Advectra

## Run scheme test for Burgers equation
domain = Domain(1, 1, 1, 1; dealiased=false)
u0 = [1.0;;]

L(u, d, p, t) = p["lambda"] * u

N(u, d, p, t) = zero(u)

# Parameters
parameters = Dict("lambda" => -1.0)

t_span = [0, 10]

prob = SpectralODEProblem(L, N, u0, domain, t_span; p=parameters, dt=0.0001)

## Solve and plot
cd(relpath(@__DIR__, pwd()))
output = Output(prob, -1, []; store_hdf=false)

sol = spectral_solve(prob, MSS2(), output)

plot(stack(sol.u)[1, 1, :]; xlabel=L"X axis", ylabel=L"Y axis")

analytical_solution(u, domain, p, t) = [u0 * exp(p["lambda"] * t);;]

## Time convergence test
timesteps = [
    2^-3,
    2^-4,
    2^-5,
    2^-6,
    2^-7,
    2^-8,
    2^-9,
    2^-10,
    2^-11,
    2^-12,
    2^-13,
    2^-14,
    2^-15,
    2^-16,
    2^-17,
    2^-18,
    2^-19,
    2^-20
]
#timesteps = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6]
_, convergence1 = test_timestep_convergence(prob, analytical_solution, timesteps, MSS1())
_, convergence2 = test_timestep_convergence(prob, analytical_solution, timesteps, MSS2())
_, convergence3 = test_timestep_convergence(prob, analytical_solution, timesteps, MSS3())
p = plot(timesteps, convergence1; xaxis=:log, yaxis=:log10, label=" MSS1")
plot!(timesteps, convergence2; xaxis=:log, yaxis=:log10, label=" MSS2", color="dark green")
plot!(timesteps, convergence3; xaxis=:log, yaxis=:log10, label=" MSS3", color="orange")
plot!(timesteps, 0.0002 * timesteps; linestyle=:dash, label=" " * L"\mathcal{O}(dt)")
plot!(timesteps, 0.0001 * timesteps .^ 2; linestyle=:dash, label=" " * L"\mathcal{O}(dt^2)")
plot!(timesteps, 0.0001 * timesteps .^ 3; linestyle=:dash, label=" " * L"\mathcal{O}(dt^3)")
plot!(; xlabel=L"dt", ylabel=L"||U-u_a||", title="Timestep convergence (Scheme test)",
      xticks=timesteps[1:2:end],
      xscale=:log2, legend_positions=:bottomright, size=(3.37 * 100, 2.5 * 100),
      minorticks=1, yticks=[1e-20, 1e-17, 1e-14, 1e-11, 1e-8, 1e-5], titlefontsize=10)
savefig("Timestep convergence, exponential test exact.png")

## Export data
using JLD

jldopen("schemetest exact ic.jld", "w") do file
    g = create_group(file, "data")
    g["convergence1"] = convergence1
    g["convergence2"] = convergence2
    g["convergence3"] = convergence3
    g["timesteps"] = timesteps
    g["colors"] = "#" .* hex.(getindex.(p.series_list[1:end], :seriescolor))
end
