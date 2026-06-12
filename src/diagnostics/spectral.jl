# ------------------------------------------------------------------------------------------
#                                   Spectral Diagnostics                                    
# ------------------------------------------------------------------------------------------

# ----------------------------------- Raw Spectral Data ------------------------------------

"""
    get_modes(state_hat, prob, time; axis::Symbol=:both)

  Return `state_hat`, along an `axis`.
  
  ### `axis` options (`Symbol`):
  - `kx`: gets modes along the `kx` axis in spectral space.
  - `ky`: gets modes along the `ky` axis in spectral space.
  - `both`: gets all the modes in spectral space.
  - `kx`: gets modes along the `kx=ky` line in spectral space. 
""" # Interface
function get_modes(state_hat::AbstractArray, prob, time; axis::Symbol=:both)
    get_modes(state_hat, prob, time, Val(axis))
end

# Catch-all (invalid axis)
function get_modes(state_hat::AbstractArray, prob, time, axis::Val{T}) where {T}
    throw(ArgumentError("axis has to be either :kx, :ky, :both or :diag, instead :$axis was given."))
end

# Specializations
get_modes(state_hat::AbstractArray, prob, time, ::Val{:kx}) = selectdim(state_hat, 1, 1)
get_modes(state_hat::AbstractArray, prob, time, ::Val{:ky}) = selectdim(state_hat, 2, 1)
get_modes(state_hat::AbstractArray, prob, time, ::Val{:both}) = @views state_hat

function get_modes(state_hat::AbstractArray, prob, time, ::Val{:diag})
    if ndims(state_hat) > 2
        return stack(diag.(eachslice(state_hat; dims=ndims(prob.domain) + 1)))
    elseif ndims(state_hat) == 2
        return diag(state_hat)
    else
        error("axis=:diag is not supported for $(ndims(state_hat))D-Arrays.")
    end
end

function build_diagnostic(::Val{:get_modes}; axis=:both, kwargs...)
    Diagnostic(; name="Modes",
               method=get_modes,
               metadata="Modes (Complex) along $axis axis",
               assumes_spectral_state=true,
               args=(Val(axis),))
end

"""
    get_log_modes(state_hat, prob, time; axis::Symbol=:diag)
  
  Return log(|`state_hat`|) along an `axis`. See [`get_modes`](@ref) for the `axis` options.
"""
function get_log_modes(state_hat, prob, time; axis::Symbol=:diag)
    get_log_modes(state_hat, prob, time, Val(axis))
end

function get_log_modes(state_hat, prob, time, axis::Val{T}) where {T}
    modes = get_modes(state_hat, prob, time, axis)
    log.(abs.(modes))
end

function build_diagnostic(::Val{:get_log_modes}; axis=:diag, kwargs...)
    Diagnostic(; name="Log modes",
               method=get_log_modes,
               metadata="log(|modes|) along $axis axis",
               assumes_spectral_state=true,
               args=(Val(axis),))
end

# ---------------------------------- Spectral Utilities ------------------------------------

_spectral_sum(f, A::AbstractArray, dims, ::FFTPlans) = sum(f, A; dims)

function _spectral_sum(f, A::AbstractArray, dims::Colon, plan::rFFTPlans)
    S = 2 * sum(f, A) - sum(f, selectdim(A, 1, 1:1))
    remove_nyquist = iseven(size(fwd(plan), 1)) && size(A, 1) > 1
    remove_nyquist ? S - sum(f, selectdim(A, 1, size(A, 1):size(A, 1))) : S
end

function _spectral_sum(f, A::AbstractArray, dims, plan::rFFTPlans)
    !(1 in dims) && return sum(f, A; dims)
    S = 2 .* sum(f, A; dims) .- sum(f, selectdim(A, 1, 1:1); dims)
    remove_nyquist = iseven(size(fwd(plan), 1)) && size(A, 1) > 1
    remove_nyquist ? S .- sum(f, selectdim(A, 1, size(A, 1):size(A, 1)); dims) : S
end

# User interface
function spectral_sum(f, A::AbstractArray, trait::AbstractTransformPlans; dims=:)
    _spectral_sum(f, A, dims, trait)
end

function spectral_sum(f, A::AbstractArray, domain::AbstractDomain; dims=:)
    _spectral_sum(f, A, dims, get_transform_plans(domain))
end

# Default f = identity
function spectral_sum(A::AbstractArray, trait::AbstractTransformPlans; dims=:)
    spectral_sum(identity, A, trait; dims)
end

function spectral_sum(A::AbstractArray, domain::AbstractDomain; dims=:)
    spectral_sum(identity, A, domain; dims)
end

export spectral_sum

# ------------------------------------ Energy Spectra --------------------------------------

"""
    energy_spectrum(power_spectrum::AbstractArray, prob, time, ::Val{spectrum})

  Computes the energy spectrum `E(k)`, based on the `spectrum` argument.

  ### `spectrum` options:
  - `:radial`: radial (kx) spectrum, averaged over the poloidal direction.
  - `:poloidal`: poloidal (ky) spectrum, averaged over the radial direction. 
  - `:wavenumber`: wavenumber (k) spectrum, averaged over wavenumber magnitude |k|.
"""
function energy_spectrum(power_spectrum::AbstractArray, prob, time, ::Val{:radial})
    @unpack domain = prob
    return vec(spectral_sum(power_spectrum, domain; dims=1)) * differential_area(domain) /
           (2π * domain.Ly * length(domain)) # domain.kx,
end

function energy_spectrum(power_spectrum::AbstractArray, prob, time, ::Val{:poloidal})
    @unpack domain = prob
    return vec(spectral_sum(power_spectrum, domain; dims=2)) * differential_area(domain) /
           (2π * domain.Lx * length(domain)) # prob.domain.ky, 
end

function energy_spectrum(power_spectrum::AbstractArray{<:Number,2},
                         prob, time, ::Val{:wavenumber})
    @unpack domain = prob

    # Determine dk for binning [k-dk, k+dk] inspired by Camargo & Durran
    dk = 0.5 * max(2π / domain.Lx, 2π / domain.Ly)
    # Compute magnitudes
    k_magnitude = hypot.(domain.kx', domain.ky)
    # Determine masks
    k_values = reshape((0:cld(maximum(k_magnitude), dk)) .* dk, 1, 1, :)
    k_magnitude = reshape(k_magnitude, size(k_magnitude)..., 1)
    masks = @. (k_values - dk ≤ k_magnitude) & (k_magnitude < k_values + dk)
    # Reshape ps for GPU optimization   
    ps = reshape(power_spectrum, size(power_spectrum)..., 1)
    # Compute energy spectrum (S(k)k/2π = (∫S(k, θ)dθ, θ∈[0, 2π])/2π)
    k_sums = vec(sum(k_magnitude .* masks; dims=(1, 2)))
    ps_sums = vec(sum(ps .* masks; dims=(1, 2)))
    denominator = vec(sum(masks; dims=(1, 2))) .^ 2

    # Serves to average both k and power_spectrum
    E = @. k_sums * ps_sums / max(denominator, 1)
    # Return spectrum alongside wavenumbers
    return E * (differential_area(domain) / (2π * length(domain))) #k_values, 
end

"""
    wavenumber_metadata(::Val{spectrum}) where spectrum<:Symbol

  Return human readable metadata about which wavenumber is stored.
"""
wavenumber_metadata(::Val{:radial}) = "Radial wavenumber (kx);"
wavenumber_metadata(::Val{:poloidal}) = "Poloidal wavenumber (ky);"
wavenumber_metadata(::Val{:wavenumber}) = "Wavenumbers (k);"

# Catch-all
function wavenumber_metadata(::Val{spectrum}) where {spectrum}
    throw(ArgumentError("spectrum has to be either :radial, :poloidal or :wavenumber, :" *
                        string(spectrum) * " was given."))
end

# ----------------------------------- Potential/density ------------------------------------

"""
    potential_energy_spectrum(state_hat, prob, time, spectrum=Val(:radial))
  
  Computes energy spectrum of the potential power spectrum |̂n(k)|², based on `spectrum` type.
  
  See [`energy_spectrum`](@ref) for `spectrum` type options.
"""
function potential_energy_spectrum(state_hat, prob, time, spectrum=Val(:radial))
    @unpack domain = prob
    n_hat = selectdim(state_hat, ndims(state_hat), 1)
    energy_spectrum(abs2.(n_hat), prob, time, spectrum) / 2
end

function potential_energy_spectrum(state_hat::AbstractArray, prob, time;
                                   spectrum::Symbol=:radial)
    potential_energy_spectrum(state_hat::AbstractArray, prob, time, Val(spectrum))
end

const density_energy_spectrum = potential_energy_spectrum

function build_diagnostic(::Union{Val{:potential_energy_spectrum},
                                  Val{:density_energy_spectrum}};
                          spectrum::Symbol=:wavenumber, kwargs...)
    start = spectrum == :wavenumber ? "Potential" :
            titlecase(string(spectrum)) * " potential"
    metadata = wavenumber_metadata(Val(spectrum)) * " Potential energy spectrum ($spectrum)"
    Diagnostic(; name=start * " energy spectrum",
               method=potential_energy_spectrum,
               metadata=metadata,
               assumes_spectral_state=true,
               args=(Val(spectrum),))
end

# ---------------------------------------- Kinetic -----------------------------------------

"""
    kinetic_energy_spectrum(state_hat, prob, time, spectrum=Val(:radial))
  
  Computes energy spectrum of the kinetic power spectrum |̂Ω(k)|², based on `spectrum` type.
  
  See [`energy_spectrum`](@ref) for `spectrum` type options.
"""
function kinetic_energy_spectrum(state_hat::AbstractArray, prob, time,
                                 spectrum=Val{:radial})
    @unpack domain, operators = prob
    @unpack solve_phi, diff_x, diff_y = operators
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    ϕ_hat = solve_phi(n_hat, Ω_hat)

    energy_spectrum(abs2.(diff_x(ϕ_hat)) + abs2.(diff_y(ϕ_hat)), prob, time, spectrum) / 2
end

function kinetic_energy_spectrum(state_hat, prob, time; spectrum::Symbol=:radial)
    kinetic_energy_spectrum(state_hat, prob, time, Val(spectrum))
end

function requires_operator(::Val{:kinetic_energy_spectrum}; kwargs...)
    [OperatorRecipe(:solve_phi), OperatorRecipe(:diff_x), OperatorRecipe(:diff_y)]
end

function build_diagnostic(::Val{:kinetic_energy_spectrum}; spectrum::Symbol=:wavenumber,
                          kwargs...)
    start = spectrum == :wavenumber ? "Kinetic" : titlecase(string(spectrum)) * " kinetic"
    metadata = wavenumber_metadata(Val(spectrum)) * " Kinetic energy spectrum ($spectrum)"
    Diagnostic(; name=start * " energy spectrum",
               method=kinetic_energy_spectrum,
               metadata=metadata,
               assumes_spectral_state=true,
               args=(Val(spectrum),))
end

# ----------------------------------------- Flux -------------------------------------------

"""
    flux_spectrum(state_hat, prob, time, spectrum=Val(:poloidal))
  
  Computes flux decomposition spectrum of the radial flux Γ=nvₓ, based on `spectrum` type.
  
  See [`energy_spectrum`](@ref) for `spectrum` type options.
"""
function flux_spectrum(state_hat::AbstractArray, prob, time, spectrum=Val{:poloidal})
    @unpack domain, operators = prob
    @unpack solve_phi, diff_y = operators
    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    ϕ_hat = solve_phi(n_hat, Ω_hat)
    vx_hat = -diff_y(ϕ_hat)
    energy_spectrum(real(n_hat .* conj.(vx_hat)), prob, time, spectrum)
end

function flux_spectrum(state_hat, prob, time; spectrum::Symbol=:poloidal)
    flux_spectrum(state_hat, prob, time, Val(spectrum))
end

function requires_operator(::Val{:flux_spectrum}; kwargs...)
    [OperatorRecipe(:solve_phi), OperatorRecipe(:diff_y)]
end

function build_diagnostic(::Val{:flux_spectrum}; spectrum::Symbol=:poloidal, kwargs...)
    start = spectrum == :wavenumber ? "Flux" : titlecase(string(spectrum)) * " flux"
    metadata = wavenumber_metadata(Val(spectrum)) * " Flux decomposition ($spectrum)"
    Diagnostic(; name=start * " spectrum",
               method=flux_spectrum,
               metadata=metadata,
               assumes_spectral_state=true,
               args=(Val(spectrum),))
end

# --------------------------------------- Enstrophy ----------------------------------------

"""
    enstrophy_spectrum(state_hat, prob, time, spectrum=Val(:radial))
  
  Computes energy spectrum of the enstrophy power spectrum |̂Ω(k)|², based on `spectrum` type.
  
  See [`energy_spectrum`](@ref) for `spectrum` type options.
"""
function enstrophy_spectrum(state_hat::AbstractArray, prob, time, spectrum=Val{:radial})
    @unpack domain = prob
    Ω_hat = eachslice(state_hat; dims=ndims(state_hat))[2]
    energy_spectrum(abs2.(Ω_hat), prob, time, spectrum) / 2
end

function enstrophy_spectrum(state_hat, prob, time; spectrum::Symbol=:radial)
    enstrophy_spectrum(state_hat, prob, time, Val(spectrum))
end

function build_diagnostic(::Val{:enstrophy_spectrum}; spectrum::Symbol=:wavenumber,
                          kwargs...)
    start = spectrum == :wavenumber ? "Enstrophy" :
            titlecase(string(spectrum)) * " enstrophy"
    metadata = wavenumber_metadata(Val(spectrum)) * " Enstrophy spectrum ($spectrum)"
    Diagnostic(; name=start * " spectrum",
               method=enstrophy_spectrum,
               metadata=metadata,
               assumes_spectral_state=true,
               args=(Val(spectrum),))
end

# -------------------------------- Electrostatic Potential ---------------------------------

"""
    electrostatic_potential_spectrum(state_hat, prob, time, spectrum=Val(:radial))
  
  Computes energy spectrum of the electrostatic potential |̂ϕ(k)|², based on `spectrum` type.
  
  See [`energy_spectrum`](@ref) for `spectrum` type options.
"""
function electrostatic_potential_spectrum(state_hat::AbstractArray, prob, time,
                                          spectrum=Val{:radial})
    @unpack domain, operators = prob
    @unpack solve_phi = operators

    slices = eachslice(state_hat; dims=ndims(state_hat))
    n_hat = slices[1]
    Ω_hat = slices[2]
    ϕ_hat = solve_phi(n_hat, Ω_hat)
    energy_spectrum(abs2.(ϕ_hat), prob, time, spectrum) / 2
end

function electrostatic_potential_spectrum(state_hat, prob, time; spectrum::Symbol=:radial)
    electrostatic_potential_spectrum(state_hat, prob, time, Val(spectrum))
end

function requires_operator(::Val{:electrostatic_potential_spectrum}; kwargs...)
    [OperatorRecipe(:solve_phi)]
end

function build_diagnostic(::Val{:electrostatic_potential_spectrum};
                          spectrum::Symbol=:wavenumber, kwargs...)
    start = spectrum == :wavenumber ? "Electrostatic potential" :
            titlecase(string(spectrum)) * " electrostatic potential"
    metadata = wavenumber_metadata(Val(spectrum)) *
               " Electrostatic potential spectrum ($spectrum)"
    Diagnostic(; name=start * " spectrum",
               method=electrostatic_potential_spectrum,
               metadata=metadata,
               assumes_spectral_state=true,
               args=(Val(spectrum),))
end

# -------------------------------------- Temperature ---------------------------------------
# TODO implement temperature energy spectrum