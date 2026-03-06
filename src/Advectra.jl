module Advectra

using FFTW, HDF5, H5Zblosc, LinearAlgebra, LaTeXStrings, MuladdMacro, UnPack, Base.Threads,
      Dates, Printf, GPUArrays, Adapt, Statistics
export @unpack

# TODO make ext
using Plots

include("operators/fftutilities.jl")
export spectral_transform, spectral_transform!, get_fwd, get_bwd

include("domains/domain.jl")
export Domain, wave_vectors, get_points, spectral_size, spectral_length, get_transform_plans

include("operators/spectralOperators.jl")
export OperatorRecipe, build_operators, build_operator # TODO perhaps remove and swap with @op
# reciprocal, spectral_exp, spectral_expm1,
# spectral_log, hyper_diffusion

using ProgressMeter, Interpolations
include("diagnostics/diagnostics.jl")
export radial_density_profile, poloidal_density_profile, radial_vorticity_profile,
       poloidal_vorticity_profile, poloidal_vorticity_profile, @diagnostics

include("spectralODEProblem.jl")
export SpectralODEProblem

include("schemes.jl")
export MSS1, MSS2, MSS3

include("outputer.jl")
export Output

include("spectralSolve.jl")
export spectral_solve

include("utilities.jl")
export initial_condition, gaussian, log_gaussian, sinusoidal, sinusoidalX, sinusoidalY,
       gaussianWallX, gaussianWallY, quadratic_function, exponential_background, randomIC,
       random_phase, random_crossphased, isolated_blob, frequencies, remove_zonal_modes!,
       remove_streamer_modes!, remove_asymmetric_modes!, remove_nothing, add_constant!,
       add_constant, send_mail, logspace
end