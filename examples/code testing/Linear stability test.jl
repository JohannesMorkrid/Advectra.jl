## Run all (alt+enter)
using Advectra
using CUDA

# Run alternative linear stability test
domain = Domain(256, 256; Lx=48, Ly=48, MemoryType=CuArray)
ic = initial_condition(random_crossphased, domain; value=1e-3)

# Linear operator
function Linear(du, u, operators, p, t)
    @unpack laplacian, diff_x, diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ν, μ, ζ = p
    dη .= ν * laplacian(η) - ζ * diff_y(η) - 2 * ν * κ * diff_x(η)
    dΩ .= μ * laplacian(Ω)
end

# Non-linear operator
function NonLinear(du, u, operators, p, t)
    @unpack solve_phi, diff_x, diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, κ, ν = p
    ϕ = solve_phi(Ω)
    dη .= -(κ - ζ) * diff_y(ϕ) + σ * ϕ
    dΩ .= -ζ * diff_y(η) + σ * ϕ
end

# Parameters
κ = 1e-2
parameters = (κ=κ, ζ=1e-3, σ=1e-3, ν=1e-4, μ=1e-4)

# Time interval
tspan = [0.0, 10000.0]

# Array of diagnostics want
diagnostics = @diagnostics [
    progress(; stride=1000),
    cfl(; stride=100, component=:radial, silent=true),
    plot_density(; stride=10000),
    get_modes(; stride=100),
    get_log_modes(; stride=100, axis=:ky)
]

# The problem
prob = SpectralODEProblem(Linear, NonLinear, ic, domain, tspan; p=parameters, dt=1e-1,
                          diagnostics=diagnostics)

# The output
output_file_name = joinpath(@__DIR__, "output", "linear-stability.h5")
output = Output(prob; filename=output_file_name, simulation_name=:parameters,
                storage_limit="100 MB", resume=true)

# Solve
sol = spectral_solve(prob, MSS3(), output;)

## ----------------------------------- Mode Andalysis --------------------------------------

using LaTeXStrings
using Plots
log_modes = sol.simulation["Log modes/data"]
t = sol.simulation["Log modes/t"][:]
gamma = (log_modes[:, :, end] - log_modes[:, :, end-1]) / (diff(t)[end])
@unpack ζ, κ, σ, ν = parameters
w0 = sqrt(ζ * κ) * sqrt(1 - ζ / κ)
plot(gamma[:, 1] / w0; xaxis=:log, xlabel=L"k_y = k_x", ylabel=L"\gamma/\gamma_0",
     ylim=[-2, 1])
vline!([(σ / ν)^(1 / 4)])
