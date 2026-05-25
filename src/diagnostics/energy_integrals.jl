# ------------------------------------------------------------------------------------------
#                                Energy Integral Diagnostics                                
# ------------------------------------------------------------------------------------------

# ----------------------------- Parseval's Theorem Utilities -------------------------------

"""
    parseval_integral(u_hat::AbstractArray, domain::Domain; average=true)
 
  Compute an integral using Parseval's theorem by working in spectral space.
  
  Uses spectral_sum to properly account for Hermitian symmetry in real transforms.
  
  Returns the integral ∫∫ |a|² dxdy if average=false,
  or the average (integral/area) if average=true.
"""
function parseval_integral(u_hat::AbstractArray, domain::Domain; average=true)
    integral = spectral_sum(abs2, u_hat, domain) / length(domain)^2
    return average ? integral : area(domain) * integral
end

"""
    parseval_integral(a_hat::AbstractArray, b_hat::AbstractArray, domain::Domain; 
                      average=true)
 
  Compute the inner product using Parseval's theorem:
  ∫∫ a(x,y) · b̄(x,y) dxdy = ∫∫ â(k) · b̂*(k) dk
  
  Uses spectral_sum to properly account for Hermitian symmetry in real transforms.
  
  Returns the integral ∫∫ a·conj(b) dxdy if average=false,
  or the average (integral/area) if average=true.
"""
function parseval_integral(a_hat::AbstractArray, b_hat::AbstractArray, domain::Domain;
                           average=true)
    integral = spectral_sum(a_hat .* conj(b_hat), domain) / length(domain)^2
    return average ? integral : area(domain) * integral
end

const parsevals_theorem = parseval_integral

# ------------------------------ Integral Of Quadratic Term --------------------------------

# TODO move or remove
function integral_of_quadratic_term(u, v, domain, quadratic_term; average=true)
    @unpack transforms, U, V, up, vp, padded, dealiasing_coefficient = quadratic_term
    mul!(U, bwd(transforms), padded ? pad!(up, u, typeof(transforms)) : u)
    mul!(V, bwd(transforms), padded ? pad!(vp, v, typeof(transforms)) : v)
    @. U *= V
    # ∫∫ UV dxdy ≈ ∑ UV dxdy
    integral = dealiasing_coefficient .* sum(U) * differential_area(domain)
    return average ? integral / area(domain) : integral
end

# ----------------------------------- Energy Integrals -------------------------------------

# ------------------------------- Potential Energy Integral --------------------------------

# P(t) = ∫dx 1/2n^2
function potential_energy_integral(state_hat, prob, time)
    @unpack domain = prob
    n_hat = selectdim(state_hat, ndims(state_hat), 1)
    parseval_integral(n_hat, domain) / 2
end

function build_diagnostic(::Val{:potential_energy_integral}; kwargs...)
    Diagnostic(; name="Potential energy integral",
               method=potential_energy_integral,
               metadata="Potential energy density.",
               assumes_spectral_state=true)
end

# -------------------------------- Kinetic Energy Integral ---------------------------------

# K(t) = ∫1/2(∇_⟂Φ)^2 = ∫dx1/2 U_E^2
function kinetic_energy_integral(state_hat, prob, time)
    @unpack domain, operators = prob
    @unpack solve_phi, diff_x, diff_y = operators
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    ϕ_hat = solve_phi(n_hat, Ω_hat)
    K = parseval_integral(diff_x(ϕ_hat), domain) .+ parseval_integral(diff_y(ϕ_hat), domain)
    return K / 2
end

function requires_operator(::Val{:kinetic_energy_integral}; kwargs...)
    [OperatorRecipe(:solve_phi), OperatorRecipe(:diff_x), OperatorRecipe(:diff_y)]
end

function build_diagnostic(::Val{:kinetic_energy_integral}; kwargs...)
    Diagnostic(; name="Kinetic energy integral",
               method=kinetic_energy_integral,
               metadata="Kinetic energy density.",
               assumes_spectral_state=true)
end

# --------------------------------- Total Energy Integral ----------------------------------

# E(t) = P(T) + K(T)
function total_energy_integral(state_hat, prob, time)
    potential_energy_integral(state_hat, prob, time) .+
    kinetic_energy_integral(state_hat, prob, time)
end

function build_diagnostic(::Val{:total_energy_integral}; kwargs...)
    Diagnostic(; name="Total energy integral",
               method=total_energy_integral,
               metadata="Total energy density.",
               assumes_spectral_state=true)
end

# ------------------------------- Enstrophy Energy Integral --------------------------------

# U(t) = ∫1/2(∇_⟂^2Φ)^2 = ∫dx1/2 Ω^2
function enstrophy_energy_integral(state_hat, prob, time)
    @unpack domain = prob
    Ω_hat = selectdim(state_hat, ndims(state_hat), 2)
    parseval_integral(Ω_hat, domain) / 2
end

function build_diagnostic(::Val{:enstrophy_energy_integral}; kwargs...)
    Diagnostic(; name="Enstrophy energy integral",
               method=enstrophy_energy_integral,
               metadata="Enstrophy energy density.",
               assumes_spectral_state=true)
end

# -------------------------------------- Radial Flux ---------------------------------------

"""
    radial_flux(state_hat, prob, time)
  
  Computes the radial (x-direction) particle flux:
  ```math
    \\Gamma_x(t) = \\frac{1}{L_xL_y}\\int_0^{L_x}\\int_0^{L_y}nv_x dxdy,
  ```

  Uses Parseval's theorem to avoid FFTs, by computing directly in spectral space.
"""
function radial_flux(state_hat::AbstractArray, prob, time)
    @unpack domain, operators = prob
    @unpack solve_phi, diff_y = operators

    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    dϕ_hat = -solve_phi(n_hat, Ω_hat)
    diff_y(dϕ_hat, dϕ_hat)

    parseval_integral(n_hat, dϕ_hat, domain)
end

function requires_operator(::Val{:radial_flux}; kwargs...)
    [OperatorRecipe(:solve_phi), OperatorRecipe(:diff_y)]
end

function build_diagnostic(::Val{:radial_flux}; kwargs...)
    Diagnostic(; name="Radial flux",
               method=radial_flux,
               metadata="Average radial particle flux",
               assumes_spectral_state=true)
end

# ------------------------------------- Poloidal Flux --------------------------------------

"""
    poloidal_flux(state_hat, prob, time)

  Computes the poloidal (y-direction) particle flux:
  ```math
    \\Gamma_y(t) = \\frac{1}{L_xL_y}\\int_0^{L_x}\\int_0^{L_y}nv_y dxdy,
  ```

  Uses Parseval's theorem to avoid FFTs, by computing directly in spectral space.
"""
function poloidal_flux(state_hat::AbstractArray, prob, time)
    @unpack domain, operators = prob
    @unpack solve_phi, diff_x = operators

    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    dϕ_hat = solve_phi(n_hat, Ω_hat)
    diff_x(dϕ_hat, dϕ_hat)

    parseval_integral(n_hat, dϕ_hat, domain)
end

function requires_operator(::Val{:poloidal_flux}; kwargs...)
    [OperatorRecipe(:solve_phi), OperatorRecipe(:diff_x)]
end

function build_diagnostic(::Val{:poloidal_flux}; kwargs...)
    Diagnostic(; name="Poloidal flux",
               method=poloidal_flux,
               metadata="Average poloidal particle flux",
               assumes_spectral_state=true)
end

# ------------------------------------- Flux Magnitude -------------------------------------

"""
    flux_magnitude(state_hat, prob, time)

  Computes the magnitude of the total flux:
  ```math
    |\\Gamma|(t) = \\sqrt{\\Gamma_x^2 + \\Gamma_y^2}
  ```
"""
function flux_magnitude(state_hat::AbstractArray, prob, time)
    Γ_x = radial_flux(state_hat, prob, time)
    Γ_y = poloidal_flux(state_hat, prob, time)
    return sqrt(Γ_x^2 + Γ_y^2)
end

function requires_operator(::Val{:flux_magnitude}; kwargs...)
    vcat(requires_operator(Val(:radial_flux); kwargs...),
         requires_operator(Val(:poloidal_flux); kwargs...))
end

function build_diagnostic(::Val{:flux_magnitude}; kwargs...)
    Diagnostic(; name="Flux magnitude",
               method=flux_magnitude,
               metadata="Magnitude of total particle flux",
               assumes_spectral_state=true)
end

# ----------------------------- Dissipative Energy Integrals -------------------------------

# ---------------------------- Resistive Dissipation Integral ------------------------------

# Γ_c(t) = C∫(n-ϕ)^2
function resistive_dissipation_integral(state_hat, prob, time; adiabaticity_symbol=:C)
    @unpack domain, operators, p = prob
    @unpack solve_phi = operators
    C = getfield(p, adiabaticity_symbol)
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    h_hat = n_hat .- solve_phi(n_hat, Ω_hat)
    return C * parseval_integral(h_hat, domain)
end

function requires_operator(::Val{:resistive_dissipation_integral}; kwargs...)
    [OperatorRecipe(:solve_phi)]
end

function build_diagnostic(::Val{:resistive_dissipation_integral}; adiabaticity_symbol=:C,
                          kwargs...)
    diagnostic_kwargs = (; adiabaticity_symbol=adiabaticity_symbol)
    Diagnostic(; name="Resistive dissipation integral",
               method=resistive_dissipation_integral,
               metadata="Resistive dissipation energy density.",
               assumes_spectral_state=true,
               kwargs=diagnostic_kwargs)
end

# ---------------------------- Potential Dissipation Integral ------------------------------

# D^E_N(t) = ν∫n∇⁶_⟂n
function potential_dissipation_integral(state_hat, prob, time; diffusivity_symbol=:ν)
    @unpack domain, p, operators = prob
    @unpack hyper_laplacian = operators
    ν = getfield(p, diffusivity_symbol)
    n_hat = eachslice(state_hat; dims=ndims(state_hat))[1]
    ν * parseval_integral(n_hat, hyper_laplacian(n_hat), domain)
end

function requires_operator(::Val{:potential_dissipation_integral}; order=3, kwargs...)
    [OperatorRecipe(:laplacian; order=order, alias=:hyper_laplacian)]
end

function build_diagnostic(::Val{:potential_dissipation_integral}; diffusivity_symbol=:ν,
                          kwargs...)
    diagnostic_kwargs = (; diffusivity_symbol=diffusivity_symbol)
    Diagnostic(; name="Potential dissipation integral",
               method=potential_dissipation_integral,
               metadata="Potential energy dissipation density.",
               assumes_spectral_state=true,
               kwargs=diagnostic_kwargs)
end

# ----------------------------- Kinetic Dissipation Integral -------------------------------

# D^E_V(t) = μ∫ϕ∇⁶_⟂Ω = μ∫(∇²_⟂Ω)² 
function kinetic_dissipation_integral(state_hat, prob, time; viscosity_symbol=:μ)
    @unpack domain, p, operators = prob
    @unpack solve_phi, hyper_laplacian = operators
    μ = getfield(p, viscosity_symbol)
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    ϕ_hat = solve_phi(n_hat, Ω_hat)
    μ * parseval_integral(ϕ_hat, hyper_laplacian(Ω_hat), domain)
end

function requires_operator(::Val{:kinetic_dissipation_integral}; order=3, kwargs...)
    [OperatorRecipe(:solve_phi),
     OperatorRecipe(:laplacian; order=order, alias=:hyper_laplacian)]
end

function build_diagnostic(::Val{:kinetic_dissipation_integral}; viscosity_symbol=:μ,
                          kwargs...)
    diagnostic_kwargs = (; viscosity_symbol=viscosity_symbol)
    Diagnostic(; name="Kinetic dissipation integral",
               method=kinetic_dissipation_integral,
               metadata="Kinetic energy dissipation density.",
               assumes_spectral_state=true,
               kwargs=diagnostic_kwargs)
end

# ----------------------------- Viscous Dissipation Integral -------------------------------

# D^E(t) = D^E_N(t) + D^E_V(t) 
function viscous_dissipation_integral(state_hat, prob, time; diffusivity_symbol=:ν,
                                      viscosity_symbol=:μ)
    potential_dissipation_integral(state_hat, prob, time; diffusivity_symbol) .+
    kinetic_dissipation_integral(state_hat, prob, time; viscosity_symbol)
end

function requires_operator(::Val{:viscous_dissipation_integral}; order=3, kwargs...)
    vcat(requires_operator(Val(:potential_dissipation_integral); order=order, kwargs...),
         requires_operator(Val(:kinetic_dissipation_integral); order=order, kwargs...))
end

function build_diagnostic(::Val{:viscous_dissipation_integral}; diffusivity_symbol=:ν,
                          viscosity_symbol=:μ, kwargs...)
    diagnostic_kwargs = (; diffusivity_symbol=diffusivity_symbol,
                         viscosity_symbol=viscosity_symbol)
    Diagnostic(; name="Viscous dissipation integral",
               method=viscous_dissipation_integral,
               metadata="Viscous energy dissipation density.",
               assumes_spectral_state=true,
               kwargs=diagnostic_kwargs)
end

# ---------------------------- Enstrophy Dissipation Integral ------------------------------

# D^U(t) = ∫(n-Ω)(ν∇⁶_⟂n - μ∇⁶_⟂Ω)
function enstrophy_dissipation_integral(state_hat, prob, time; diffusivity_symbol=:ν,
                                        viscosity_symbol=:μ, kwargs...)
    @unpack domain, p, operators = prob
    @unpack hyper_laplacian = operators
    ν = getfield(p, diffusivity_symbol)
    μ = getfield(p, viscosity_symbol)
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    h_hat = n_hat - Ω_hat
    diffusive_terms_hat = ν * hyper_laplacian(n_hat) - μ * hyper_laplacian(Ω_hat)
    parseval_integral(h_hat, diffusive_terms_hat, domain)
end

function requires_operator(::Val{:enstrophy_dissipation_integral}; kwargs...)
    [OperatorRecipe(:laplacian; order=3, alias=:hyper_laplacian)]
end

function build_diagnostic(::Val{:enstrophy_dissipation_integral}; diffusivity_symbol=:ν,
                          viscosity_symbol=:μ, kwargs...)
    diagnostic_kwargs = (; diffusivity_symbol=diffusivity_symbol,
                         viscosity_symbol=viscosity_symbol)
    Diagnostic(; name="Enstrophy dissipation integral",
               method=enstrophy_dissipation_integral,
               metadata="Enstrophy energy dissipation density.",
               assumes_spectral_state=true,
               kwargs=diagnostic_kwargs)
end

# ---------------------------------- Evolution Integrals -----------------------------------

# ------------------------------- Energy Evolution Integral --------------------------------

# dE/dt(t) = Γ_n - Γ_c - D^E 
function energy_evolution_integral(state_hat, prob, time; adiabaticity_symbol=:C,
                                   diffusivity_symbol=:ν, viscosity_symbol=:μ)
    radial_flux(state_hat, prob, time) .-
    resistive_dissipation_integral(state_hat, prob, time; adiabaticity_symbol) .-
    viscous_dissipation_integral(state_hat, prob, time; diffusivity_symbol,
                                 viscosity_symbol)
end

function requires_operator(::Val{:energy_evolution_integral}; order=3, kwargs...)
    vcat(requires_operator(Val(:radial_flux); kwargs...),
         requires_operator(Val(:resistive_dissipation_integral); kwargs...),
         requires_operator(Val(:viscous_dissipation_integral); order=order, kwargs...))
end

function build_diagnostic(::Val{:energy_evolution_integral}; adiabaticity_symbol=:C,
                          diffusivity_symbol=:ν,
                          viscosity_symbol=:μ, kwargs...)
    diagnostic_kwargs = (; adiabaticity_symbol=adiabaticity_symbol,
                         diffusivity_symbol=diffusivity_symbol,
                         viscosity_symbol=viscosity_symbol)
    Diagnostic(; name="Energy evolution integral",
               method=energy_evolution_integral,
               metadata="Energy density evolution.",
               assumes_spectral_state=true,
               kwargs=diagnostic_kwargs)
end

# ------------------------------- Enstrophy Energy Integral --------------------------------

#dU / dt(t) = Γ_n - D^U
function enstrophy_evolution_integral(state_hat, prob, time; diffusivity_symbol=:ν,
                                      viscosity_symbol=:μ)
    radial_flux(state_hat, prob, time) .-
    enstrophy_dissipation_integral(state_hat, prob, time; diffusivity_symbol,
                                   viscosity_symbol)
end

function requires_operator(::Val{:enstrophy_evolution_integral}; order=3, kwargs...)
    vcat(requires_operator(Val(:radial_flux); kwargs...),
         requires_operator(Val(:enstrophy_dissipation_integral); order=order, kwargs...))
end

function build_diagnostic(::Val{:enstrophy_evolution_integral}; diffusivity_symbol=:ν,
                          viscosity_symbol=:μ, kwargs...)
    diagnostic_kwargs = (; diffusivity_symbol=diffusivity_symbol,
                         viscosity_symbol=viscosity_symbol)
    Diagnostic(; name="Enstrophy evolution integral",
               method=enstrophy_evolution_integral,
               metadata="Enstrophy density evolution.",
               assumes_spectral_state=true,
               kwargs=diagnostic_kwargs)
end