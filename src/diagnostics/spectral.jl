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

# ------------------------------------ Energy Spectra --------------------------------------

"""
    energy_spectrum(power_spectrum::AbstractArray, prob, time, ::Val{spectrum})

  Computes the energy spectrum `E(k)`, based on the `spectrum` argument.

  ### `spectrum` options:
  - `:radial`: radial (kx) spectrum, averaged over the poloidal direction.
  - `:poloidal`: poloidal (ky) spectrum, averaged over the radial direction. 
  - `:wavenumber`: wavenumber (k) spectrum, averaged over wavenumber magnitude |k|.

  ## Returns:
    (wavenumbers, E) where E = E(k) (`Tuple`).
"""
function energy_spectrum(power_spectrum::AbstractArray, prob, time, ::Val{:radial})
    return vec(sum(power_spectrum; dims=1)) * (2 * pi / prob.domain.Ly)
end

function energy_spectrum(power_spectrum::AbstractArray, prob, time, ::Val{:poloidal})
    return vec(sum(power_spectrum; dims=2)) * (2 * pi / prob.domain.Lx)
end

# TODO add windowed option?
function energy_spectrum(power_spectrum::AbstractArray{<:Number,2},
                         prob, time, ::Val{:wavenumber})
    @unpack domain = prob

    # Determine dk for binning [k-dk, k+dk] inspired by Camargo
    dk = 0.5 * min(2 * pi / domain.Lx, 2 * pi / domain.Ly)
    # Compute magnitudes
    k_magnitude = hypot.(domain.kx', domain.ky)
    # Determine bins
    nbins = ceil(Int, maximum(k_magnitude) / dk) # Or = max(size(domain)...)
    k_values = (0:nbins) .* dk
    # Compute energy spectrum (S(k) = 2π*(∫S(k, θ)dθ, θ∈[0, 2π])/2π)
    E = 2pi .* [mean(power_spectrum[k-dk.<=k_magnitude.<k+dk]) for k in k_values]
    # Return spectrum alongside wavenumbers
    return E
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

# --------------------------------------- Potential ----------------------------------------

"""
    potential_energy_spectrum(state_hat, prob, time, spectrum=Val(:radial))
  
  Computes energy spectrum of the potential power spectrum |̂n(k)|², based on `spectrum` type.
  
  See [`energy_spectrum`](@ref) for `spectrum` type options.
"""
function potential_energy_spectrum(state_hat, prob, time, spectrum=Val(:radial))
    @unpack domain = prob
    n_hat = selectdim(state_hat, ndims(state_hat), 1)
    energy_spectrum(abs2.(n_hat), prob, time, spectrum)
end

function potential_energy_spectrum(state_hat::AbstractArray, prob, time;
                                   spectrum::Symbol=:radial)
    potential_energy_spectrum(state_hat::AbstractArray, prob, time, Val(spectrum))
end

function build_diagnostic(::Val{:potential_energy_spectrum}; spectrum::Symbol=:wavenumber,
                          kwargs...)
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
    @unpack domain = prob
    Ω_hat = selectdim(state_hat, ndims(state_hat), 1)
    energy_spectrum(abs2.(Ω_hat), prob, time, spectrum)
end

function kinetic_energy_spectrum(state_hat, prob, time; spectrum=:radial)
    kinetic_energy_spectrum(state_hat, prob, time, Val(spectrum))
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
