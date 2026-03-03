using Advectra
using CUDA

domain = Domain(1024, 1024; Lx=50, Ly=50, MemoryType=CuArray, precision=Float64)

# 1) N = cst => ϖ = N∇²ϕ

# Density
N = 10 * CUDA.ones(size(domain))
N_hat = spectral_transform(N, get_fwd(domain))

# Construct dipole vorticity
Ω = initial_condition(gaussian, domain) |> CuArray
Ω_hat = spectral_transform(Ω, get_fwd(domain))
diff_y = build_operator(:diff_y, domain)
ϖ_hat = diff_y(Ω_hat)

solve_phi_b = build_operator(:solve_phi, domain; boussinesq=true)

# Required to build non-Boussinesq operator
diff_x = build_operator(:diff_x, domain)
quadratic_term = build_operator(:quadratic_term, domain)

solve_phi_nb = build_operator(:solve_phi, domain; boussinesq=false, diff_x, diff_y,
                              quadratic_term, density=:linear)

# Solve for Potential                              
ϕ1_hat = solve_phi_b(ϖ_hat)
ϕ2_hat = solve_phi_nb(N_hat, ϖ_hat)

# Build Laplacian operator
laplacian = build_operator(:laplacian, domain)

# Compute ∇²ϕ
ϖ1_hat = laplacian(ϕ1_hat)
ϖ1 = spectral_transform(ϖ1_hat, get_bwd(domain))
ϖ2_hat = laplacian(ϕ2_hat)
ϖ2 = spectral_transform(ϖ2_hat, get_bwd(domain))

# To be compared against
ddϖ = spectral_transform(ϖ_hat, get_bwd(domain))

using Plots
heatmap(domain, Array(ϖ1); aspect_ratio=:equal, title="Boussinesq ∇²ϕ")
heatmap(domain, Array(N .* ϖ2); aspect_ratio=:equal, title="Non-Boussinesq N∇²ϕ")
heatmap(domain, Array(ddϖ); aspect_ratio=:equal, title="Reference ϖ")
heatmap(domain, Array(N); aspect_ratio=:equal, title="Density")

# 2) ϕ = f(x,y) => ϖ(x,y) = ∇⋅(N∇ϕ)

function vorticity(N, ϕ; diff_x, diff_y, quadratic_term)
    diff_x(quadratic_term(N, diff_x(ϕ))) + diff_y(quadratic_term(N, diff_y(ϕ)))
end

N = initial_condition(gaussian, domain; A=10, B=1) |> CuArray
spectral_transform!(N_hat, get_fwd(domain), N)
ϕ = CUDA.@allowscalar sin.(domain.ky[2] * domain.y) .+ 0 * domain.x' |> CuArray
ϕ_hat = spectral_transform(ϕ, get_fwd(domain))
ϖ_hat = vorticity(N_hat, ϕ_hat; diff_x=diff_x, diff_y, quadratic_term)

ϖ = spectral_transform(ϖ_hat, get_bwd(domain))
heatmap(Array(ϖ))

ϕ1_hat = solve_phi_b(ϖ_hat)
ϕ2_hat = solve_phi_nb(N_hat, ϖ_hat)
ϕ1 = spectral_transform(ϕ1_hat, get_bwd(domain))
ϕ2 = spectral_transform(ϕ2_hat, get_bwd(domain))
heatmap(Array(ϕ1))
heatmap(Array(ϕ2))
using LinearAlgebra
norm(ϕ2 - ϕ)