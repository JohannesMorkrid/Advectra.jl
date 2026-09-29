# ------------------------------------------------------------------------------------------
#                                         Profiles
# ------------------------------------------------------------------------------------------

# ---------------------------------------- Helpers -----------------------------------------

"""
    radial_profile(field)

  Average `field` over the poloidal (y) direction: `f₀(x) = 1/L_y∫_0^L_y f(x,y)dy`.
"""
radial_profile(field::AbstractArray) = vec(mean(field; dims=1))

"""
    poloidal_profile(field)

  Average `field` over the radial (x) direction: `f₀(y) = 1/L_x∫_0^L_x f(x,y)dx`.
"""
poloidal_profile(field::AbstractArray) = vec(mean(field; dims=2))

# ---------------------------------------- Density -----------------------------------------

function radial_density_profile(state, prob, time)
    radial_profile(selectdim(state, ndims(prob.domain) + 1, 1))
end

function build_diagnostic(::Val{:radial_density_profile}; kwargs...)
    Diagnostic(; name="Radial density profile",
               method=radial_density_profile,
               metadata="Poloidally averaged density n₀(x)")
end

function poloidal_density_profile(state, prob, time)
    poloidal_profile(selectdim(state, ndims(prob.domain) + 1, 1))
end

function build_diagnostic(::Val{:poloidal_density_profile}; kwargs...)
    Diagnostic(; name="Poloidal density profile",
               method=poloidal_density_profile,
               metadata="Radially averaged density n₀(y)")
end

# --------------------------------------- Vorticity ----------------------------------------

function radial_vorticity_profile(state, prob, time)
    radial_profile(selectdim(state, ndims(prob.domain) + 1, 2))
end

function build_diagnostic(::Val{:radial_vorticity_profile}; kwargs...)
    Diagnostic(; name="Radial vorticity profile",
               method=radial_vorticity_profile,
               metadata="Poloidally averaged vorticity Ω₀(x)")
end

function poloidal_vorticity_profile(state, prob, time)
    poloidal_profile(selectdim(state, ndims(prob.domain) + 1, 2))
end

function build_diagnostic(::Val{:poloidal_vorticity_profile}; kwargs...)
    Diagnostic(; name="Poloidal vorticity profile",
               method=poloidal_vorticity_profile,
               metadata="Radially averaged vorticity Ω₀(y)")
end

# ------------------------------------------ Flux ------------------------------------------

"""
    radial_flux_profile(state_hat, prob, time)

  Computes the poloidally averaged radial particle flux `Γ₀(x) = 1/L_y∫_0^L_y nv_x dy`,
  where `v_x = -∂ϕ/∂y`.
"""
function radial_flux_profile(state_hat, prob, time)
    @unpack domain, operators = prob
    @unpack solve_phi, diff_y = operators

    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    n = bwd_plan(domain) * n_hat
    v_x = bwd_plan(domain) * -diff_y(solve_phi(n_hat, Ω_hat))

    radial_profile(n .* v_x)
end

function requires_operator(::Val{:radial_flux_profile}; kwargs...)
    [OperatorRecipe(:solve_phi), OperatorRecipe(:diff_y)]
end

function build_diagnostic(::Val{:radial_flux_profile}; kwargs...)
    Diagnostic(; name="Radial flux profile",
               method=radial_flux_profile,
               metadata="Poloidally averaged radial particle flux Γ₀(x)",
               assumes_spectral_state=true)
end
