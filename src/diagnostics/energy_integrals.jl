# ------------------------------------------------------------------------------------------
#                                Energy Integral Diagnostics                                
# ------------------------------------------------------------------------------------------

# ----------------------------------- Parsevals Theorem ------------------------------------

# For the most part uses Parsevals theorem, these does not correct for aliasing
function parsevals_theorem(coeffs::AbstractArray, domain::Domain; compute_density=true)
    _parsevals_theorem(coeffs, domain, Val(domain.real_transform); compute_density)
end

function _parsevals_theorem(coeffs::AbstractArray, domain::Domain,
                            real_transform::Val{true}; compute_density=true)
    integral = @views (2 * sum(abs2.(coeffs)) .- sum(abs2.(coeffs[1,:]))) *
                      differential_area(domain)) / length(domain)
    return compute_density ? integral / area(domain) : integral
end

function _parsevals_theorem(coeffs::AbstractArray, domain::Domain,
                            real_transform::Val{false}; compute_density=true)
    integral = @views (sum(abs2.(coeffs))) * differential_area(domain) / length(domain)
    return compute_density ? integral / area(domain) : integral
end

# ------------------------------ Integral Of Quadratic Term --------------------------------

function integral_of_quadratic_term(u, v, domain, quadratic_term; compute_density=true)
    @unpack transforms, U, V, up, vp, padded, dealiasing_coefficient = quadratic_term
    mul!(U, bwd(transforms), padded ? pad!(up, u, typeof(transforms)) : u)
    mul!(V, bwd(transforms), padded ? pad!(vp, v, typeof(transforms)) : v)
    @. U *= V
    # ∫∫Udxdy ≈ ∑Udxdy
    integral = dealiasing_coefficient * sum(U) * differential_area(domain)
    return compute_density ? integral / area(domain) : integral
end

# ------------------------------- Potential Energy Integral --------------------------------

# P(t) = ∫dx 1/2n^2
function potential_energy_integral(state_hat, prob, time; quadrature=nothing)
    @unpack domain = prob
    n_hat = selectdim(state_hat, ndims(domain) + 1, 1)
    parsevals_theorem(n_hat, domain) / 2
end

function build_diagnostic(::Val{:potential_energy_integral}; kwargs...)
    Diagnostic(; name="Potential energy integral",
               method=potential_energy_integral,
               metadata="Potential energy density.",
               assumes_spectral_state=true)
end

# ----------------------------------- Energy Integrals -------------------------------------

# -------------------------------- Kinetic Energy Integral ---------------------------------

# K(t) = ∫1/2(∇_⟂Φ)^2 = ∫dx1/2 U_E^2
function kinetic_energy_integral(state_hat, prob, time; quadrature=nothing)
    @unpack domain, operators = prob
    @unpack solve_phi, diff_x, diff_y = operators
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    ϕ_hat = solve_phi(n_hat, Ω_hat)
    parsevals_theorem(diff_x(ϕ_hat), domain) / 2 + parsevals_theorem(diff_y(ϕ_hat), domain) / 2
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

# --------------------------------- Total Enegy Integral -----------------------------------

# E(t) = P(T) + K(T)
function total_energy_integral(state, prob, time; quadrature=nothing)
    potential_energy_integral(state, prob, time; quadrature=quadrature) .+
    kinetic_energy_integral(state, prob, time; quadrature=quadrature)
end

function build_diagnostic(::Val{:total_energy_integral}; kwargs...)
    Diagnostic(; name="Total energy integral",
               method=total_energy_integral,
               metadata="Total energy density.",
               assumes_spectral_state=true)
end

# ------------------------------- Enstropy Energy Integral ---------------------------------

# U(t) = ∫1/2(∇_⟂^2Φ)^2 = ∫dx1/2 Ω^2
function enstropy_energy_integral(state_hat, prob, time; quadrature=nothing)
    @unpack domain = prob
    Ω_hat = selectdim(state_hat, ndims(domain) + 1, 2)
    parsevals_theorem(Ω_hat, domain) / 2
end

function build_diagnostic(::Val{:enstropy_energy_integral}; kwargs...)
    Diagnostic(; name="Enstropy energy integral",
               method=enstropy_energy_integral,
               metadata="Enstropy energy density.",
               assumes_spectral_state=true)
end

# ----------------------------- Dissipative Energy Integrals -------------------------------

# ---------------------------- Resistive Dissipation Integral ------------------------------

# Γ_c(t) = C∫(n-ϕ)^2
function resistive_dissipation_integral(state_hat, prob, time; adiabaticity_symbol=:C,
                                        quadrature=nothing)
    @unpack domain, operators, p = prob
    @unpack solve_phi, quadratic_term = operators
    C = getfield(p, adiabaticity_symbol)
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    h_hat = n_hat .- solve_phi(n_hat, Ω_hat)
    return C * parsevals_theorem(h_hat, domain)
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
function potential_dissipation_integral(state, prob, time; diffusivity_symbol=:ν,
                                        quadrature=nothing)
    @unpack domain, p, operators = prob
    @unpack hyper_laplacian, quadratic_term = operators
    ν = getfield(p, diffusivity_symbol)
    slices = eachslice(state; dims=ndims(state))
    n = slices[1]
    ν * integral_of_quadratic_term(n, hyper_laplacian(n), domain, quadratic_term)
end

function requires_operator(::Val{:potential_dissipation_integral}; kwargs...)
    [OperatorRecipe(:laplacian; order=3, alias=:hyper_laplacian),
     OperatorRecipe(:quadratic_term)]
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

# D^E_V(t) = μ∫ϕ∇⁶_⟂Ω = μ∫Ω∇⁴_⟂Ω 
function kinetic_dissipation_integral(state, prob, time; viscosity_symbol=:μ,
                                      quadrature=nothing)
    @unpack domain, p, operators = prob
    @unpack solve_phi, hyper_laplacian, quadratic_term = operators
    μ = getfield(p, viscosity_symbol)
    slices = eachslice(state; dims=ndims(state))
    n = slices[1]
    Ω = slices[2]
    ϕ = solve_phi(n, Ω)
    μ * integral_of_quadratic_term(ϕ, hyper_laplacian(Ω), domain, quadratic_term)
end

function requires_operator(::Val{:kinetic_dissipation_integral}; kwargs...)
    [OperatorRecipe(:solve_phi), OperatorRecipe(:quadratic_term),
     OperatorRecipe(:laplacian; order=3, alias=:hyper_laplacian)]
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
function viscous_dissipation_integral(state, prob, time; diffusivity_symbol=:ν,
                                      viscosity_symbol=:μ, quadrature=nothing)
    potential_dissipation_integral(state, prob, time; diffusivity_symbol=diffusivity_symbol,
                                   quadrature=quadrature) .+
    kinetic_dissipation_integral(state, prob, time; viscosity_symbol=viscosity_symbol,
                                 quadrature=quadrature)
end

function requires_operator(::Val{:viscous_dissipation_integral}; kwargs...)
    vcat(requires_operator(Val(:potential_dissipation_integral); kwargs...),
         requires_operator(Val(:kinetic_dissipation_integral); kwargs...))
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

# ----------------------------- Enstropy Dissipation Integral ------------------------------

# D^U(t) = ∫(n-Ω)(ν∇⁶_⟂n - μ∇⁶_⟂Ω)
function enstropy_dissipation_integral(state, prob, time; diffusivity_symbol=:ν,
                                       viscosity_symbol=:μ, kwargs...)
    @unpack domain, p, operators = prob
    @unpack hyper_laplacian, quadratic_term = operators
    ν = getfield(p, diffusivity_symbol)
    μ = getfield(p, viscosity_symbol)
    slices = eachslice(state; dims=ndims(state))
    n = slices[1]
    Ω = slices[2]
    h = n - Ω
    diffusive_terms = ν * hyper_laplacian(n) - μ * hyper_laplacian(Ω)
    integral_of_quadratic_term(h, diffusive_terms, domain, quadratic_term)
end

function requires_operator(::Val{:enstropy_dissipation_integral}; kwargs...)
    [OperatorRecipe(:laplacian; order=3, alias=:hyper_laplacian),
     OperatorRecipe(:quadratic_term)]
end

function build_diagnostic(::Val{:enstropy_dissipation_integral}; diffusivity_symbol=:ν,
                          viscosity_symbol=:μ, kwargs...)
    diagnostic_kwargs = (; diffusivity_symbol=diffusivity_symbol,
                         viscosity_symbol=viscosity_symbol)
    Diagnostic(; name="Enstropy dissipation integral",
               method=enstropy_dissipation_integral,
               metadata="Enstropy energy dissipation density.",
               assumes_spectral_state=true,
               kwargs=diagnostic_kwargs)
end

# ---------------------------------- Evolution Integrals -----------------------------------

# ------------------------------- Energy Evolution Integral --------------------------------

# dE/dt(t) = Γ_n - Γ_c - D^E 
function energy_evolution_integral(state, prob, time; adiabaticity_symbol=:C,
                                   diffusivity_symbol=:ν, viscosity_symbol=:μ,
                                   quadrature=nothing)
    radial_flux(state, prob, time; quadrature=quadrature) .-
    resistive_dissipation_integral(state, prob, time;
                                   adiabaticity_symbol=adiabaticity_symbol,
                                   quadrature=quadrature) .-
    viscous_dissipation_integral(state, prob, time; diffusivity_symbol=diffusivity_symbol,
                                 viscosity_symbol=viscosity_symbol, quadrature=quadrature)
end

function requires_operator(::Val{:energy_evolution_integral}; kwargs...)
    vcat(requires_operator(Val(:radial_flux); kwargs...),
         requires_operator(Val(:resistive_dissipation_integral); kwargs...),
         requires_operator(Val(:viscous_dissipation_integral); kwargs...))
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

# ------------------------------- Enstropy Energy Integral ---------------------------------

#dU / dt(t) = Γ_n - D^U
function enstropy_evolution_integral(state, prob, time; diffusivity_symbol=:ν,
                                     viscosity_symbol=:μ, quadrature=nothing)
    radial_flux(state, prob, time; quadrature=quadrature) .-
    enstropy_dissipation_integral(state, prob, time; diffusivity_symbol=diffusivity_symbol,
                                  viscosity_symbol=viscosity_symbol, quadrature=quadrature)
end

function requires_operator(::Val{:enstropy_evolution_integral}; kwargs...)
    vcat(requires_operator(Val(:radial_flux); kwargs...),
         requires_operator(Val(:enstropy_dissipation_integral); kwargs...))
end

function build_diagnostic(::Val{:enstropy_evolution_integral}; diffusivity_symbol=:ν,
                          viscosity_symbol=:μ, kwargs...)
    diagnostic_kwargs = (; diffusivity_symbol=diffusivity_symbol,
                         viscosity_symbol=viscosity_symbol)
    Diagnostic(; name="Enstropy evolution integral",
               method=enstropy_evolution_integral,
               metadata="Enstropy density evolution.",
               assumes_spectral_state=true,
               kwargs=diagnostic_kwargs)
end
