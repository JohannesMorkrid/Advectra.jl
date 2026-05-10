# ------------------------------------------------------------------------------------------
#                                       FFT Utilities                                       
# ------------------------------------------------------------------------------------------

# ---------------------------------- FFT Transform Plans -----------------------------------

abstract type AbstractTransformPlans end

abstract type FourierTransformPlans <: AbstractTransformPlans end

"""
    FFTPlans{FT<:FFTW.Plan, iFT<:FFTW.Plan} 
  Collection of transform plans using standard FFT.
"""
struct FFTPlans{FWD<:FFTW.Plan,BWD<:FFTW.Plan} <: FourierTransformPlans
    FT::FWD
    iFT::BWD
end

"""
    rFFTPlans{FT<:FFTW.Plan, iFT<:FFTW.Plan}
  Collection of transform plans using real FFT (rFFT), utilizing hermitian symmetry.
"""
struct rFFTPlans{FWD<:FFTW.Plan,BWD<:FFTW.Plan} <: FourierTransformPlans
    FT::FWD
    iFT::BWD
end

# ------------------------------ Spectral Transform Methods --------------------------------

"""
    spectral_transform!(U<:AbstractArray, transformplan<:FFTW.Plan)
    spectral_transform!(U<:Union{Tuple,Vector}, transformplan<:FFTW.Plan)
  Spectral transform, applies transform plan p to u in-place returning du.
"""
spectral_transform!(du, p::P, u) where {P<:FFTW.Plan} = _spectral_transform!(du, p, u)

function _spectral_transform!(du, p::P, u::AbstractArray{<:Number}) where {P<:FFTW.Plan}
    idx = ntuple(_ -> :, ndims(p))
    for i in axes(u, ndims(p) + 1)
        mul!(view(du, idx..., i), p, u[idx..., i])
    end
end

function _spectral_transform!(du, p::P,
                              u::AbstractArray{<:AbstractArray}) where {P<:FFTW.Plan}
    for i in eachindex(u)
        _spectral_transform!(du[i], p, u[i])
    end
end

# TODO optimize spectral_transform methods, maybe use spectral_transform! TODO or remove
"""
    spectral_transform(U<:AbstractArray, transformplan<:FFTW.Plan)
    spectral_transform(U<:Union{Tuple,Vector}, transformplan<:FFTW.Plan)
  Spectral transform out-of-place, applies transform plan p to U. 
"""
function spectral_transform(U::T, p::P) where {T<:AbstractArray,P<:FFTW.Plan}
    mapslices(u -> p * u, U; dims=(1, 2))
end

function spectral_transform(U::T, p::P) where {T<:Union{Tuple,Vector},P<:FFTW.Plan}
    map(u -> p * u, U)
end

# function spectral_transform(U, p::P) where {P<:FFTW.Plan}
#     allocate_coefficients()
#     spectral_transform!
# end

# ---------------------------------------- Helpers -----------------------------------------

get_fwd(transformplans::FourierTransformPlans) = transformplans.FT
const fwd = get_fwd

get_bwd(transformplans::FourierTransformPlans) = transformplans.iFT
const bwd = get_bwd

"""
    Base.show(io::IO, transformplans::AbstractTransformPlans)

  Pretty-print `AbstractTransformPlans`.
"""
function Base.show(io::IO, transformplans::AbstractTransformPlans)
    typename = nameof(typeof(transformplans))

    fwd = get_fwd(transformplans)
    bwd = get_bwd(transformplans)

    # TODO perhaps do a better check if the plan is real
    if transformplans isa rFFTPlans
        println(io, typename)
        println(io, "  forward (real→complex): ", fwd)
        println(io, "  backward (complex→real): ", bwd)
    else
        println(io, typename)
        println(io, "  forward: ", fwd)
        println(io, "  backward:", bwd)
    end
end

"""
    Base.show(io::IO, ::MIME"text/plain", transformplans::AbstractTransformPlans)

  Compact one-line show of AbstractTransformPlans for use in arrays and etc.
"""
function Base.show(io::IO, ::MIME"text/plain", transformplans::AbstractTransformPlans)
    typename = nameof(typeof(transformplans))
    fwd = get_fwd(transformplans)
    bwd = get_bwd(transformplans)

    print(io, typename, "(fwd: ", fwd, ", bwd: ", bwd, ")")
end