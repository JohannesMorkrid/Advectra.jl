# ------------------------------------------------------------------------------------------
#                            Gradient-Driven Sheath Interchange                             
# ------------------------------------------------------------------------------------------

# ---------------------------------- Bohm Normalization ------------------------------------

# ----------------------------------- Nonlinear Sheaths ------------------------------------

function SI_GD_B(du, u, operators, p, t)
    @unpack solve_phi, diff_x, diff_y = operators
    @unpack poisson_bracket, grad_dot_grad, spectral_expm1 = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, ν, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (κ - ζ) * diff_y(ϕ) - ζ * diff_y(η) +
          ν * grad_dot_grad(η, η) - 2ν * κ * diff_x(η) - σ * spectral_expm1(-ϕ)
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) - σ * spectral_expm1(-ϕ)

    CUDA.@allowscalar dη[1] = 0
    CUDA.@allowscalar dΩ[1] = 0
end

# ------------------------------- Different Damping Models ---------------------------------

function SI_GD_NPD_B(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, diff_x, diff_y, spectral_expm1 = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, ν, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (κ - ζ) * diff_y(ϕ) - ζ * diff_y(η) -
          2ν * κ * diff_x(η) - σ * spectral_expm1(-ϕ)
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) - σ * spectral_expm1(-ϕ)

    CUDA.@allowscalar dη[1] = 0
    CUDA.@allowscalar dΩ[1] = 0
end

function SI_GD_NX_B(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, grad_dot_grad, diff_y, spectral_expm1 = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, ν, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (κ - ζ) * diff_y(ϕ) - ζ * diff_y(η) +
          ν * grad_dot_grad(η, η) - σ * spectral_expm1(-ϕ)
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) - σ * spectral_expm1(-ϕ)

    CUDA.@allowscalar dη[1] = 0
    CUDA.@allowscalar dΩ[1] = 0
end

function SI_GD_SD_B(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, diff_y, spectral_expm1 = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (κ - ζ) * diff_y(ϕ) - ζ * diff_y(η) -
          σ * spectral_expm1(-ϕ)
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) - σ * spectral_expm1(-ϕ)

    CUDA.@allowscalar dη[1] = 0
    CUDA.@allowscalar dΩ[1] = 0
end

const SI_GD_NPD_NX_B = SI_GD_SD_B

# ------------------------------------- Linear Sheath --------------------------------------

function SI_GD_LS_B(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, grad_dot_grad, diff_x, diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, ν, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (κ - ζ) * diff_y(ϕ) - ζ * diff_y(η) +
          ν * grad_dot_grad(η, η) - 2ν * κ * diff_x(η) + σ * ϕ
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) + σ * ϕ
end

# ------------------------------- Different Damping Models ---------------------------------

function SI_GD_LS_NPD_B(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, diff_x, diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, ν, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (κ - ζ) * diff_y(ϕ) - ζ * diff_y(η) -
          2ν * κ * diff_x(η) + σ * ϕ
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) + σ * ϕ
end

function SI_GD_LS_NX_B(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, grad_dot_grad, diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, ν, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (κ - ζ) * diff_y(ϕ) - ζ * diff_y(η) +
          ν * grad_dot_grad(η, η) + σ * ϕ
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) + σ * ϕ
end

function SI_GD_LS_SD_B(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (κ - ζ) * diff_y(ϕ) - ζ * diff_y(η) + σ * ϕ
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) + σ * ϕ
end

const SI_GD_LS_NPD_NX_B = SI_GD_LS_SD_B

# -------------------------------- Gyro-Bohm Normalization ---------------------------------

# ----------------------------------- Nonlinear Sheaths ------------------------------------

function SI_GD_GB(du, u, operators, p, t)
    @unpack solve_phi, diff_x, diff_y = operators
    @unpack poisson_bracket, grad_dot_grad, spectral_expm1 = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, ν, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (1 - ζ) * diff_y(ϕ) - ζ * diff_y(η) +
          ν * κ * grad_dot_grad(η, η) - 2ν * κ * diff_x(η) -
          (σ / κ) * spectral_expm1(-κ * ϕ)
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) - (σ / κ) * spectral_expm1(-κ * ϕ)

    CUDA.@allowscalar dη[1] = 0
    CUDA.@allowscalar dΩ[1] = 0
end

# ------------------------------- Different Damping Models ---------------------------------

function SI_GD_NPD_GB(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, diff_x, diff_y, spectral_expm1 = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, ν, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (1 - ζ) * diff_y(ϕ) - ζ * diff_y(η) -
          2ν * κ * diff_x(η) - (σ / κ) * spectral_expm1(-κ * ϕ)
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) - (σ / κ) * spectral_expm1(-κ * ϕ)

    CUDA.@allowscalar dη[1] = 0
    CUDA.@allowscalar dΩ[1] = 0
end

function SI_GD_NX_GB(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, grad_dot_grad, diff_y, spectral_expm1 = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, ν, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (1 - ζ) * diff_y(ϕ) - ζ * diff_y(η) +
          ν * κ * grad_dot_grad(η, η) - (σ / κ) * spectral_expm1(-κ * ϕ)
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) - (σ / κ) * spectral_expm1(-κ * ϕ)

    CUDA.@allowscalar dη[1] = 0
    CUDA.@allowscalar dΩ[1] = 0
end

function SI_GD_SD_GB(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, diff_y, spectral_expm1 = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (1 - ζ) * diff_y(ϕ) - ζ * diff_y(η) -
          (σ / κ) * spectral_expm1(-κ * ϕ)
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) - (σ / κ) * spectral_expm1(-κ * ϕ)

    CUDA.@allowscalar dη[1] = 0
    CUDA.@allowscalar dΩ[1] = 0
end

const SI_GD_NPD_NX_GB = SI_GD_SD_GB

# ------------------------------------- Linear Sheath --------------------------------------

function SI_GD_LS_GB(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, grad_dot_grad, diff_x, diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, ν, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (1 - ζ) * diff_y(ϕ) - ζ * diff_y(η) +
          ν * κ * grad_dot_grad(η, η) - 2ν * κ * diff_x(η) + σ * ϕ
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) + σ * ϕ
end

# ------------------------------- Different Damping Models ---------------------------------

function SI_GD_LS_NPD_GB(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, diff_x, diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, ν, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (1 - ζ) * diff_y(ϕ) - ζ * diff_y(η) -
          2ν * κ * diff_x(η) + σ * ϕ
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) + σ * ϕ
end

function SI_GD_LS_NX_GB(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, grad_dot_grad, diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ, ν, κ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (1 - ζ) * diff_y(ϕ) - ζ * diff_y(η) +
          ν * κ * grad_dot_grad(η, η) + σ * ϕ
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) + σ * ϕ
end

function SI_GD_LS_SD_GB(du, u, operators, p, t)
    @unpack solve_phi, poisson_bracket, diff_y = operators
    η, Ω = eachslice(u; dims=3)
    dη, dΩ = eachslice(du; dims=3)
    @unpack ζ, σ = p
    ϕ = solve_phi(η, Ω)

    dη .= poisson_bracket(η, ϕ) - (1 - ζ) * diff_y(ϕ) - ζ * diff_y(η) + σ * ϕ
    dΩ .= poisson_bracket(Ω, ϕ) - ζ * diff_y(η) + σ * ϕ
end

const SI_GD_LS_NPD_NX_GB = SI_GD_LS_SD_GB