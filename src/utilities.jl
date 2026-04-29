# ------------------------------ Initial condition helpers ---------------------------------

# Default trait
broadcastable_ic(::Function) = Val(true)

# TODO add @nobroadcast macro to make function broadcastable_ic(::typeof(func)) = Val(false)

function initial_condition(f::Function, domain::AbstractDomain; kwargs...)
    initial_condition(broadcastable_ic(f), f, domain; kwargs...)
end

function initial_condition(::Val{true}, f::Function, domain::AbstractDomain; kwargs...)
    f.(domain.x', domain.y; kwargs...)
end

function initial_condition(::Val{false}, f::Function, domain::AbstractDomain; kwargs...)
    f(domain; kwargs...)
end

# ------------------------------- Initial conditions ---------------------------------------

function gaussian(x, y; A=1, B=0, l=1, lx=l, ly=l, x0=0, y0=0)
    B + A * exp(-(x - x0)^2 / (2 * lx^2) - (y - y0)^2 / (2 * ly^2))
end

function log_gaussian(x, y; A=1, B=1, l=1, lx=l, ly=l, x0=0, y0=0)
    log(gaussian(x, y; A=A, B=B, lx=lx, ly=ly, x0=x0, y0=y0))
end

sinusoidal(x, y; Lx=1, Ly=1, n=1, m=1) = sin(2 * π * n * x / Lx) * cos(2 * π * m * y / Ly)

sinusoidalX(x, y; L=1, N=1) = sin(2 * π * N * x / L)

sinusoidalY(x, y; L=1, N=1) = sin(2 * π * N * y / L)

gaussianWallX(x, y; A=1, l=1) = A * exp(-x^2 / (2 * l^2))

gaussianWallY(x, y; A=1, l=1) = A * exp(-y^2 / (2 * l^2))

exponential_background(x, y; kappa=1) = exp(-kappa * x)

quadratic_function(x, y; kappa=1) = abs(y) <= 1 ? 1 - y .^ 2 : 0.0

randomIC(x, y) = rand()

function random_phase(domain::AbstractDomain; value=10^-6)
    θ = 2 * π * rand(domain.Ny, domain.Nx)
    n_hat = value * ones(domain.Ny, domain.Nx) .* exp.(im * θ)
    n_hat[:, 1] .= 0
    n_hat[1, :] .= 0
    real(ifft(n_hat))
end

function random_crossphased(domain::AbstractDomain; value=10^-6, cross_phase=pi / 2)
    θ = 2 * π * rand(domain.Ny, domain.Nx)
    phi_hat = value * ones(domain.Ny, domain.Nx) .* exp.(im * θ)
    n_hat = phi_hat .* exp(im * cross_phase)
    [real(ifft(n_hat));;; real(ifft(phi_hat))]
end

function isolated_blob(domain::AbstractDomain; density::Symbol=:lin, ndims=2, kwargs...)
    isolated_blob(domain, Val(density); ndims, kwargs...)
end

function isolated_blob(domain::AbstractDomain, ::Val{:lin}; ndims=2, kwargs...)
    u0 = initial_condition(gaussian, domain; kwargs...)
    ic = zeros(size(u0)..., ndims)
    ic[:, :, 1] .= u0
    return ic
end

function isolated_blob(domain::AbstractDomain, ::Val{:log}; ndims=2, kwargs...)
    u0 = initial_condition(log_gaussian, domain; kwargs...)
    ic = zeros(size(u0)..., ndims)
    ic[:, :, 1] .= u0
    return ic
end

broadcastable_ic(::typeof(random_phase)) = Val(false)
broadcastable_ic(::typeof(random_crossphased)) = Val(false)
broadcastable_ic(::typeof(isolated_blob)) = Val(false)

function isolated_temperature_blob(domain::AbstractDomain; density::Symbol=:lin, ndims=3,
                                   kwargs...)
    isolated_temperature_blob(domain, Val(density); ndims, kwargs...)
end

function isolated_temperature_blob(domain::AbstractDomain, ::Val{:lin}; ndims=3, kwargs...)
    u0 = initial_condition(gaussian, domain; kwargs...)
    ic = zeros(size(u0)..., ndims)
    ic[:, :, 1] .= 1.0
    ic[:, :, 3] .= u0
    return ic
end

function isolated_temperature_blob(domain::AbstractDomain, ::Val{:log}; ndims=3, kwargs...)
    u0 = initial_condition(log_gaussian, domain; kwargs...)
    ic = zeros(size(u0)..., ndims)
    ic[:, :, 1] .= 0.0
    ic[:, :, 3] .= u0
    return ic
end

broadcastable_ic(::typeof(isolated_temperature_blob)) = Val(false)

# ---------------------- Inverse functions / transforms ------------------------------------

expTransform(u::AbstractArray) = [exp.(u[:, :, 1]);;; u[:, :, 2]]

# ------------------------------------- Mode Related ---------------------------------------

# In-place method
function add_constant!(out::AbstractArray, field::AbstractArray, val::Number)
    out .= field
    @allowscalar out[1] += val
    return out
end

# In-place but overwrites the field
function add_constant!(field::AbstractArray, val::Number)
    @allowscalar field[1] += val
    return field
end

# Out-of place method
add_constant(field::AbstractArray, val::Number) = add_constant!(similar(field), field, val)

#------------------------------ Removal of modes -------------------------------------------

function remove_zonal_modes!(u::U, d::D) where {U<:AbstractArray,D<:AbstractDomain}
    @inbounds u[1, :, :] .= 0
end

function remove_streamer_modes!(u::U, d::D) where {U<:AbstractArray,D<:AbstractDomain}
    @inbounds u[:, 1, :] .= 0
end

function remove_asymmetric_modes!(u::U,
                                  domain::D) where {U<:AbstractArray,
                                                    D<:AbstractDomain}
    if domain.Nx % 2 == 0
        @inbounds u[:, domain.Nx÷2+1, :] .= 0
    end
    if Ny % 2 == 0
        @inbounds u[domain.Ny÷2+1, :, :] .= 0
    end
end

remove_nothing(u::U, d::D) where {U<:AbstractArray,D<:AbstractDomain} = nothing

# ------------------------------------ Other -----------------------------------------------

# For parameter scans
logspace(start, stop, length) = 10 .^ range(start, stop, length)

# TODO move to ext
# Extend plotting to allow domain as input
import Plots.plot
function plot(domain::AbstractDomain, args...; kwargs...)
    plot(domain.x, domain.y, args...; kwargs...)
end

"""
frequencies(state)

  Displays a heatmap of the mode-amplitudes using log scale.
"""
frequencies(state::AbstractArray) = heatmap(log10.(abs.(state)); title="Frequencies")

# --------------------------------------- Mailing ------------------------------------------

function send_mail(subject; attachment="")
    if length(methods(send_mail)) == 1
        error("SMTPClient is not loaded. Please add SMTPClient.jl and configure the .env file.")
    else
        throw(MethodError(send_mail, Tuple{typeof(subject)}))
    end
end
