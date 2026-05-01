# ------------------------------ Initial condition helpers ---------------------------------

"""
    broadcastable_ic(::Function) -> Val{Bool}

Trait indicating whether an initial condition function should be broadcasted over spatial coordinates.

Default `Val(true)`.
"""
broadcastable_ic(::Function) = Val(true)

"""
    @nobroadcast f(domain::AbstractDomain; kwargs...) = ...

Mark a domain-level initial condition function as non-broadcastable for use with
[`initial_condition`](@ref). By default, [`initial_condition`](@ref) broadcasts `f` pointwise over
the grid; this macro overrides that behaviour so the domain is passed to `f` instead.

Internally, this expands to:

```julia
broadcastable_ic(::typeof(f)) = Val(false)
```

Works with both `function f(...) ... end` and `f(args...) = ...` syntax.
"""
macro nobroadcast(expr)
    fname = expr
    while fname isa Expr
        fname = fname.head === :(.) ? fname.args[end].value : fname.args[1]
    end

    fname isa Symbol || throw(ArgumentError("Unsupported function signature"))

    trait = quote
        broadcastable_ic(::typeof($fname)) = Val(false)
    end

    return quote
        Base.@__doc__ $(esc(expr))
        $(esc(trait))
    end
end

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

"""
    gaussian(x, y; A=1, B=0, l=1, lx=l, ly=l, x0=0, y0=0)

2D Gaussian function:
```math 
f(x,y) = B + A\\exp{\\left(-\\frac{(x-x_0)^2}{2l_x^2}-\\frac{(y-y_0)^2}{2l_y^2}\\right)}
```

# Arguments
- `A`: Amplitude (default `1`)
- `B`: Background level (default `0`)
- `lx`, `ly`: Length scales in `x` and `y` (default (`1`, `1`)), jointly controlled via `l`
- `x0`, `y0`: centre coordinates (default (`0`,`0`))
"""
function gaussian(x, y; A=1, B=0, l=1, lx=l, ly=l, x0=0, y0=0)
    B + A * exp(-(x - x0)^2 / (2 * lx^2) - (y - y0)^2 / (2 * ly^2))
end

"""
  log_gaussian(x, y; A=1, B=1, l=1, lx=l, ly=l, x0=0, y0=0)

Logarithm of the 2D Gaussian function:
```math 
f(x,y) = \\log\\left(B + A\\exp{\\left(-\\frac{(x-x_0)^2}{2l_x^2}-\\frac{(y-y_0)^2}{2l_y^2}\\right)}\\right)
```

# Arguments
- `A`: Amplitude (default `1`)
- `B`: Background level (default `1`)
- `lx`, `ly`: Length scales in `x` and `y` (default (`1`, `1`)), jointly controlled via `l`
- `x0`, `y0`: centre coordinates (default (`0`,`0`))
"""
function log_gaussian(x, y; A=1, B=1, l=1, lx=l, ly=l, x0=0, y0=0)
    log(gaussian(x, y; A=A, B=B, lx=lx, ly=ly, x0=x0, y0=y0))
end

"""
    gaussian_x(x, y; A=1, l=1, x0=0)

2D Gaussian function varying in `x` only:
```math
f(x,y) = A\\exp{\\left(-\\frac{(x-x_0)^2}{2l^2}\\right)}
```
# Arguments
- `A`: Amplitude (default `1`)
- `l`: Length scale in `x` (default `1`)
- `x0`: Centre coordinate in `x` (default `0`)
"""
gaussian_x(x, y; A=1, l=1, x0=0) = A * exp(-(x - x0)^2 / (2 * l^2))

"""
    gaussian_y(x, y; A=1, l=1, x0=0)

2D Gaussian function varying in `y` only:
```math
f(x,y) = A\\exp{\\left(-\\frac{(y-y_0)^2}{2l^2}\\right)}
```
# Arguments
- `A`: Amplitude (default `1`)
- `l`: Length scale in `y` (default `1`)
- `y0`: Centre coordinate in `y` (default `0`)
"""
gaussian_y(x, y; A=1, l=1, y0=0) = A * exp(-(y - y0)^2 / (2 * l^2))

"""
    sinusoidal(x, y; Lx=1, Ly=1, n=1, m=1)

2D Sinusoidal function:
```math
f(x,y) = \\sin\\left(\\frac{2\\pi N_x x}{l_x}\\right)\\cos\\left(\\frac{2\\pi N_y y}{l_y}\\right)
```

# Arguments
- `lx`, `ly`: Domain sizes in x and y (default `(1, 1)`)
- `Nx`, `Ny`: Wave-numbers in x and y (default `(1, 1)`)
"""
function sinusoidal(x, y; lx=1, ly=1, Nx=1, Ny=1)
    sin(2 * π * Nx * x / lx) * cos(2 * π * Ny * y / ly)
end
"""
    sinusoidal_x(x, y; L=1, N=1)

2D Sinusoidal function varying in `x` only:
```math
f(x,y) = \\sin{\\left(\\frac{2\\pi Nx}{L}\\right)}
```

# Arguments
- `L`: Domain size in x (default `1`)
- `N`: Wave-number in x (default `1`)
"""
sinusoidal_x(x, y; L=1, N=1) = sin(2 * π * N * x / L)

"""
    sinusoidal_y(x, y; L=1, N=1)

2D Sinusoidal function varying in `y` only:
```math
f(x,y) = \\sin{\\left(\\frac{2\\pi Ny}{L}\\right)}
```

# Arguments
- `L`: Domain size in y (default `1`)
- `N`: Wave-number in y (default `1`)
"""
sinusoidal_y(x, y; L=1, N=1) = sin(2 * π * N * y / L)

"""
    exponential_x(x, y; κ=1)

2D Decaying exponential in `x` only:

```math
f(x,y) = \\exp(-\\kappa x)
```
# Arguments
- `κ`: decay rate (default `1`)
"""
exponential_x(x, y; κ=1) = exp(-κ * x)

"""
    quadratic_y(x, y)

2D Piecewise quadratic profile in `y` only:

```math
f(x,y) = \\begin{cases}
            1 - y^2, & |y| \\le 1 \\\\
            0, & |y| > 1
         \\end{cases}
```
"""
quadratic_y(x, y) = abs(y) <= 1 ? 1 - y .^ 2 : 0.0

"""
    white_noise(x, y; σ=1)

Complex-valued spatial white noise field.

```math
f(x,y) = \\sigma (a + ib), \\quad a,b \\sim \\mathcal{N}(0,1)
```
"""
white_noise(x, y; σ=1) = σ * randn(ComplexF64)

"""
    random_phase(domain::AbstractDomain; value=1e-6, ndims=1)

Generate random noise from randomly phased spectral modes, with zonal and streamer-modes removed.

```math
\\hat{u} = A\\exp(i\\theta), \\quad \\theta \\sim U[0,2\\pi)
```

# Arguments
- `value`: Amplitude of the spectral coefficients (default `1e-6`)
- `ndims`: Number of fields (default `1`)
"""
@nobroadcast function random_phase(domain::AbstractDomain; value=10^-6, ndims=1)
    θ = 2 * π * rand(spectral_size(domain)..., (ndims == 1 ? () : (ndims,))...)
    u_hat = value .* exp.(im * θ)
    u_hat[:, 1, :] .= 0
    u_hat[1, :, :] .= 0
    spectral_transform(u_hat, bwd(domain))
end

"""
    random_crossphased(domain::AbstractDomain; value=1e-6, cross_phase=π/2)

Generate random noise for the density n and vorticity Ω fields from randomly phased spectral 
    modes, with zonal and streamer modes removed.

The density nₖ is related to the potential ϕₖ via a fixed cross-phase θ:
```math
\\hat{n} = \\hat{\\phi}\\exp(i\\theta)
```
and the vorticity is derived from the potential via:
```math
\\Omega = \\nabla_\\perp^2 \\phi
```

# Arguments
- `value`: Amplitude of the spectral coefficients (default `1e-6`)
- `cross_phase`: Phase offset between ϕₖ and nₖ modes (default `π/2`)
"""
@nobroadcast function random_crossphased(domain::AbstractDomain; value=1e-6,
                                         cross_phase=π / 2)
    θ = 2 * π * rand(spectral_size(domain)...) |> domain.MemoryType
    ϕ_hat = value .* exp.(im * θ)
    n_hat = ϕ_hat .* exp(im * cross_phase)
    laplacian = build_operator(:laplacian, domain)
    Ω_hat = laplacian(ϕ_hat)
    u_hat = cat(n_hat, Ω_hat; dims=3)

    u_hat[1, :, :] .= 0
    u_hat[:, 1, :] .= 0

    spectral_transform(u_hat, bwd(domain))
end

"""
    isolated_blob(domain::AbstractDomain; density::Symbol=:lin, ndims=2, kwargs...)

Generate a 2D isolated Gaussian density blob with remaining fields set to zero.

The `density` keyword controls whether the blob is represented on a linear or
logarithmic scale. Additional keyword arguments are passed to the underlying
[`gaussian`](@ref) or [`log_gaussian`](@ref) function.

# Arguments
- `density`: Scale of the blob, either `:lin` (default) or `:log`
- `ndims`: Number of fields (default `2`)
- `kwargs...`: Passed to [`gaussian`](@ref) or [`log_gaussian`](@ref), typically `A`, `B`, `l`
"""
@nobroadcast function isolated_blob(domain::AbstractDomain; density::Symbol=:lin, kwargs...)
    isolated_blob(domain, Val(density); kwargs...)
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

"""
    isolated_temperature_blob(domain::AbstractDomain; density::Symbol=:lin, ndims=2, kwargs...)

Generate a 2D isolated Gaussian temperature blob with unity density and remaining fields set to zero.

The `density` keyword controls whether the blob is represented on a linear or
logarithmic scale. Additional keyword arguments are passed to the underlying
[`gaussian`](@ref) or [`log_gaussian`](@ref) function.

# Arguments
- `density`: Scale of the blob, either `:lin` (default) or `:log`
- `ndims`: Number of fields (default `2`)
- `kwargs...`: Passed to [`gaussian`](@ref) or [`log_gaussian`](@ref), typically `A`, `B`, `l`
"""
@nobroadcast function isolated_temperature_blob(domain::AbstractDomain; density::Symbol=:lin, ndims=3,
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
