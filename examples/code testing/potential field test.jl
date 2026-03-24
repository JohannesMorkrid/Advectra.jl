#Potential field test
using Advectra

function analytical_solution(x, y; l=1)
    r2 = x^2 + y^2
    exp(-r2 / (2 * l^2))
end

function omega(x, y; l=1)
    r2 = x^2 + y^2
    (r2 - 2 * l^2) * exp(-r2 / (2 * l^2)) / (l^4)
end

domain = Domain(1024; L=1)
l = 0.01

Ω = initial_condition(omega, domain; l=l)
using Plots
surface(Ω)

Ω_hat = get_fwd(domain) * Ω
solve_phi = Advectra.build_operator(Val(:solve_phi), domain)
ϕ_hat = solve_phi(Ω_hat)
ϕ = get_bwd(domain) * ϕ_hat
#Remove the value at boundary since should be zero for periodicity
ϕ_n = ϕ .- ϕ[1]

ϕ_a = initial_condition(analytical_solution, domain; l=l)

using LinearAlgebra
println("Norm without correcting boundary: ", norm(ϕ .- ϕ_a))
println("Norm after correcting boundary: ", norm(ϕ_n .- ϕ_a))

surface(domain, ϕ_n; xlabel="x", ylabel="y", zlabel="ϕ")
surface(domain, ϕ_a)
