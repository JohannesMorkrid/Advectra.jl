module Advectra

using FFTW, HDF5, H5Zblosc, LinearAlgebra, LaTeXStrings, MuladdMacro, UnPack, Base.Threads,
      Dates, GPUArrays, Adapt, Statistics
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
       poloidal_vorticity_profile, poloidal_vorticity_profile, ProgressDiagnostic,
       plot_frequencies, @diagnostics, DiagnosticRecipe
export cfl, radial_COM, plot_density, plot_vorticity, plot_potential,
       potential_energy_integral, kinetic_energy_integral, total_energy_integral,
       enstropy_energy_integral, resistive_dissipation_integral,
       potential_dissipation_integral, kinetic_dissipation_integral,
       viscous_dissipation_integral, enstropy_dissipation_integral,
       energy_evolution_integral, enstropy_evolution_integral, radial_flux, poloidal_flux,
       probe_density, probe_vorticity, probe_potential, probe_radial_velocity, probe_all,
       progress, get_modes, get_log_modes, potential_energy_spectrum,
       kinetic_energy_spectrum

include("spectralODEProblem.jl")
export SpectralODEProblem

include("schemes.jl")
export MSS1, MSS2, MSS3

include("outputer.jl")
export Output

include("spectralSolve.jl")
export spectral_solve

include("utilities.jl")
export initial_condition, @nobroadcast, gaussian, log_gaussian, gaussian_x, gaussian_y,
       sinusoidal, sinusoidal_x, sinusoidal_y, quadratic_y, exponential_x, white_noise,
       random_phase, random_crossphased, isolated_blob, isolated_temperature_blob, frequencies, remove_zonal_modes!,
       remove_streamer_modes!, remove_asymmetric_modes!, remove_nothing, add_constant!,
       add_constant, send_mail, logspace
end
