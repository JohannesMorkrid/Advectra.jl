## Run all (alt+enter)
using Advectra
using CUDA

# Run linear stability test
domain = Domain(256, 256; Lx=100, Ly=100, MemoryType=CuArray)
ic = initial_condition(random_crossphased, domain; value=1e-6)

# Linear operator (May not be static actually)
function Linear(du, u, operators, p, t)
    @unpack laplacian = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ν, μ = p
    dη .= ν * laplacian(η)
    dΩ .= μ * laplacian(Ω)
end

# Non-linear operator
function NonLinear(du, u, operators, p, t)
    @unpack solve_phi, diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack κ, g, σ = p
    ϕ = solve_phi(Ω)
    dη .= -(κ - g) * diff_y(ϕ) - g * diff_y(η) + σ * ϕ
    #dη .-= σ * η  # This is an additional term
    dΩ .= -g * diff_y(η) + σ * ϕ
end

# Parameters
parameters = (ν=1e-2, μ=1e-2, g=1e-3, σ=1e-5, κ=sqrt(1e-3))

# Time interval
tspan = [0.0, 500.0]

# Array of diagnostics want
diagnostics = @diagnostics [
    progress(; stride=100),
    cfl(; silent=true),
    plot_density(; stride=1000),
    get_modes(; stride=100),
    get_log_modes(; stride=100, axis=:ky)
]

# The problem
prob = SpectralODEProblem(Linear, NonLinear, ic, domain, tspan; p=parameters, dt=1e-2,
                          diagnostics=diagnostics)

# The output
output_file_name = joinpath(@__DIR__, "output", "linear-stability test Synne.h5")
output = Output(prob; filename=output_file_name, simulation_name=:parameters, resume=true)

# Solve and plot
sol = spectral_solve(prob, MSS3(), output;)

## ----------------------------------- Mode Andalysis --------------------------------------

using LaTeXStrings
using Plots
log_modes = sol.simulation["Log modes/data"]
t = sol.simulation["Log modes/t"][:]
gamma = (log_modes[:, :, end] - log_modes[:, :, end-1]) / (diff(t)[end])
@unpack g, κ = parameters
w0 = sqrt(g * κ) * sqrt(1 - g / κ)
plot(gamma[:, 1] / w0; xaxis=:log, xlabel=L"k_y = k_x", ylabel=L"\gamma/\gamma_0",
     ylim=[-2, 1])
vline!([κ^(-1 / 4)])
