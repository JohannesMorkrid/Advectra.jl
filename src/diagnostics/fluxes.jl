# ------------------------------------------------------------------------------------------
#                                     Flux Diagnostics                                      
# ------------------------------------------------------------------------------------------

# ------------------------------------ Helper Function -------------------------------------

"""
    physical_flux(n_hat::T, dϕ_hat::T, domain)
  
  Compute the flux in physical space given a 'spectal gradient' of `ϕ` and the spectral `n`.
  Dispatches on `AbstractArray` type `T` to differentiate between CPU and GPU.
"""
function physical_flux(n_hat::T, dϕ_hat::T, domain) where {T<:AbstractArray}
    # Allocate arrays
    v = zeros(size(domain)) # TODO cache these perhaps?
    n = similar(v)

    # Bwd transforms
    task_v = Threads.@spawn mul!(v, get_bwd(prob), dϕ_hat)
    task_n = Threads.@spawn mul!(n, get_bwd(prob), n_hat)
    wait(task_v)
    wait(task_n)

    # Compute Γ
    @threads for i in eachindex(n)
        @inbounds v[i] *= n[i]
    end
    return v
end

function physical_flux(n_hat::T, dϕ_hat::T, domain) where {T<:AbstractGPUArray}
    # Allocate arrays
    v = zeros(size(domain)) |> memory_type(domain) # TODO cache these perhaps?
    n = similar(v)

    mul!(n, get_bwd(domain), n_hat)
    mul!(v, get_bwd(domain), dϕ_hat)
    v .*= n
end

# -------------------------------------- Radial Flux ---------------------------------------

"""
    radial_flux(state_hat, prob, time)
  Computes 
  ```math
    \\Gamma_0(t) = \\frac{1}{L_xL_y}\\int_0^{L_x}\\int_0^{L_y} nv_x dydx,
  ```
  does not take into acount dealiasing.
"""
function radial_flux(state_hat::AbstractGPUArray, prob, time; quadrature=nothing)
    @unpack domain, operators = prob
    @unpack solve_phi, diff_y = operators
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    # Compute dϕdy 
    dϕ_hat = -solve_phi(n_hat, Ω_hat)
    diff_y(dϕ_hat, dϕ_hat)

    # Compute flux average
    Γ = physical_flux(n_hat, dϕ_hat, domain)
    return sum(Γ) / area(domain)
end

function requires_operator(::Val{:radial_flux}; kwargs...)
    [OperatorRecipe(:solve_phi), OperatorRecipe(:diff_y)]
end

function build_diagnostic(::Val{:radial_flux}; kwargs...)
    Diagnostic(; name="Radial flux",
               method=radial_flux,
               metadata="Average radial flux",
               assumes_spectral_state=true)
end

# ------------------------------------- Poloidal Flux --------------------------------------
"""
    poloidal_flux(state_hat, prob, time)

  Computes Γ_0(t) = 1/(L_xL_y)∫_0^L_x∫_0^L_y nv_y dydx, does not take into acount dealiasing.
"""
function poloidal_flux(state_hat::AbstractGPUArray, prob, time; quadrature=nothing)
    @unpack domain, operators = prob
    @unpack solve_phi, diff_x = operators
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    # Compute dϕdy 
    dϕ_hat = solve_phi(n_hat, Ω_hat)
    diff_x(dϕ_hat, dϕ_hat)

    # Compute flux average
    Γ = physical_flux(n_hat, dϕ_hat, domain)
    return sum(Γ) / area(domain)
end

function requires_operator(::Val{:poloidal_flux}; kwargs...)
    [OperatorRecipe(:solve_phi), OperatorRecipe(:diff_x)]
end

function build_diagnostic(::Val{:poloidal_flux}; kwargs...)
    Diagnostic(; name="Poloidal flux",
               method=poloidal_flux,
               metadata="Average poloidal flux",
               assumes_spectral_state=true)
end