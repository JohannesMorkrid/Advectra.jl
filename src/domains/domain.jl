# ------------------------------------- Domain ---------------------------------------------

abstract type AbstractDomain end

"""
Two dimensional `Domain` assuming bi-spectral boundary conditions, equiped with \
transform_plans` for how to transform between physical and spectral 'space'.

The `Domain` struct is a collection of hints for how to construct the [`SpectralOperator`]()s \
used in the `Linear` and `NonLinear` right hand sides (rhs) used in the [`SpectralODEProblem`](). \
In addition the `Domain` store info about the [`wave_vectors`]() associated with the `Domain``.

## Constructors
`Domain(N; L=1, kwargs...)`

`Domain(Nx, Ny; Lx=1, Ly=1, 
        MemoryType=Array, 
        precision=Float64, 
        real_transform=true,
        dealiased=true, 
        x0=-Lx / 2, y0=-Ly / 2)`
  
### Positional Arguments

#### For 'square' `Domain` specifically:

- `N`: number of points along both x- and y-direction (`Integer`).

#### For general 2D `Domain` specifically:

- `Nx`: number of points along the x-direction (`Integer`).
- `Ny`: number of points along the y-direction (`Integer`).

### Keyword Arguments

#### For 'square' `Domain` specifically:

- `L`: length of domain along both directions (`Number`).

#### For general 2D `Domain` specifically:

- `Lx`: length of domain along the x-directions (`Number`).
- `Ly`: length of domain along the y-directions (`Number`)

#### General keyword arguments

- `MemoryType`: memory `Type` used to store states, fields, spectral coefficients, etc., also \
used to configure transform plans. `Array` is used by default but all `<:AbstractArray` \
types should be supported as long as the corresponding `MemoryType` package is loaded.
- `precision`: precision `Type` used in numerical calculations, also applies to precision of \
type and output. Defaults to `Float64`. Should be a numerical `DataType`.
- `real_transform`: boolean flag to tell the program whether or not to use `rfft` and `irfft` \
methods to transform between physical and spectral space. This halves the number of spectral \
coefficient needed to be stored. Defaults to `true`.
- `dealiased`: whether or not to dealias the [Pseudo spectral](https://munin.uit.no/bitstream/handle/10037/37597/no.uit%3awiseflow%3a7267480%3a61779855.pdf#page=45)
method. Defaults to `true`.
- `x0`: x-position of the lower left corner of the `Domain`. Defaults to ``x0 = -Lx / 2``.
- `y0`: y-position of the lower left corner of the `Domain`. Defaults to ``y0 = -Ly / 2``.

## Examples

```jldoctest
julia> Domain(128; L = 10, precision = Float32)
Domain(Nx:128, Ny:128, Lx:10.0, Ly:10.0, real_transform:true, dealiased:true, MemoryType:Array{Float32}) offset by (-5.0, -5.0)

julia> Domain(128, 256; Ly = 2, real_transform = true, x0 = 0, y0 = 10)
Domain(Nx:128, Ny:256, Lx:1.0, Ly:2.0, real_transform:true, dealiased:true, MemoryType:Array{Float64}) offset by (0.0, 10.0)
```

## Fields
- `dx`: spatial resolution/grid spacing, ```dx = 2Lx÷(Nx-1)```, in x-direction.
- `dy`: spatial resolution/grid spacing, ```dy = 2Ly÷(Nx-1)```, in y-direction.
- `x`: spatial position of each point along the x-axis. Uniformly distributed by default.
- `y`: spatial position of each point along the y-axis. Uniformly distributed by default.
- `kx`: wave vector components along x-axis.
- `ky`: wave vector components along y-axis.
- `transforms`: collection of fwd and bwd transform methods.

!!! warning
    Restricted to 2D for the time being.
"""
struct Domain{X<:AbstractArray,Y<:AbstractArray,KX<:AbstractArray,
              KY<:AbstractArray,TP<:AbstractTransformPlans,
              T<:AbstractFloat} <: AbstractDomain
    Nx::Int
    Ny::Int
    Lx::T
    Ly::T
    dx::T
    dy::T
    x::X
    y::Y
    kx::KX
    ky::KY
    real_transform::Bool
    dealiased::Bool
    transforms::TP
    MemoryType::Type
    precision::DataType

    Domain(N::Integer; L::Number=1, kwargs...) = Domain(N, N; Lx=L, Ly=L, kwargs...)
    function Domain(Nx::Integer, Ny::Integer;
                    Lx::Number=1,
                    Ly::Number=1,
                    MemoryType::Type{<:AbstractArray}=Array,
                    precision::DataType=Float64,
                    real_transform::Bool=true,
                    dealiased::Bool=true,
                    x0::Number=-Lx / 2,
                    y0::Number=-Ly / 2)

        # Ensure MemoryType is not parameterized
        if MemoryType.var.name != :T
            throw(ArgumentError("MemoryType should not include type parameters (e.g., \
            Array{Float64}). Pass precision separately."))
        end

        # Ensure precision is a numeric type
        if !(precision <: Number)
            throw(ArgumentError("Precision must be a numeric type (e.g., Float64, Int)."))
        end

        # Compute step sizes
        dx = Lx / Nx
        dy = Ly / Ny

        # dx and dy is subtracted at the end, because of periodic boundary conditions
        x = LinRange{precision}(x0, x0 + Lx - dx, Nx)
        y = LinRange{precision}(y0, y0 + Ly - dy, Ny)

        # Prepare frequencies
        kx, ky = prepare_frequencies(Nx, Ny, dx, dy, MemoryType, precision, real_transform)

        # Prepare transform plans
        transform_plans = prepare_transform_plans(Nx, Ny, MemoryType, precision,
                                                  real_transform)

        new{typeof(x),typeof(y),typeof(kx),typeof(ky),typeof(transform_plans),
            precision}(Nx, Ny, Lx, Ly, dx, dy, x, y, kx, ky, real_transform, dealiased,
                       transform_plans, MemoryType, precision)
    end
end

# Helpers
function prepare_frequencies(Nx, Ny, dx, dy, MemoryType, precision, real_transform)

    # Compute frequencies, impose Hermitian symetry if real_transform
    kx = 2 * π * fftfreq(Nx, 1 / dx) |> MemoryType{precision}
    ky = 2 * π * (real_transform ? rfftfreq(Ny, 1 / dy) : fftfreq(Ny, 1 / dy)) |>
         MemoryType{precision}

    return kx, ky
end

"""
    prepare_transform_plans(Nx, Ny, use_cuda, precision, real_transform)
    
Prepare transform plan by preparing a domain using dispatching to call the right \
construction method.
"""
function prepare_transform_plans(Nx, Ny, MemoryType, precision, real_transform)

    # Temporarly create an array to create the transform plan
    utmp = zeros(precision, Ny, Nx) |> MemoryType

    # Dispatch on Domain and real_transform
    construct_transform_plans(utmp, Domain, Val(real_transform))
end

"""
    construct_transform_plans(utmp, ::Type{Domain}, ::Val{false})
    construct_transform_plans(utmp, ::Type{Domain}, ::Val{true})

Constructs transform plans based on Domain type and if Real or Complex valued fields.
"""
function construct_transform_plans(utmp, ::Type{Domain}, ::Val{true})
    FT = plan_rfft(utmp)
    iFT = plan_irfft(FT * utmp, first(size(utmp)))
    return rFFTPlans(FT, iFT)
end

function construct_transform_plans(utmp, ::Type{Domain}, ::Val{false})
    FFTPlans(plan_fft(utmp), plan_ifft(utmp))
end

# ----------------------------------- Interface --------------------------------------------

function Base.show(io::IO, m::MIME"text/plain", domain::AbstractDomain)
    typename = nameof(typeof(domain))

    # Assumes 2D domain
    labels = reverse(get_labels(domain))
    Ns = reverse(size(domain))
    Ls = reverse(lengths(domain))
    points = reverse(first.(get_points(domain)))

    print(io, typename, "(")
    print(io, join("N" .* string.(labels) .* ":" .* string.(Ns), ", "))
    print(io, ", ")
    print(io, join("L" .* string.(labels) .* ":" .* string.(Ls), ", "))

    if get(io, :compact, false)
        print(io, ")")
    else # Prints additional kwargs
        kwargs = domain_kwargs(domain)
        for (key, value) in pairs(kwargs)
            print(io, ", ", key, ":", value)
        end
        print(io, ", MemoryType:", memory_type(domain), ")")
        if any(points .!= 0)
            print(io, " offset by (", join(points, ", "), ")")
        end
    end
end

"""
lengths(domain::Domain)

Return a tuple of lengths for each axis, default (Ly, Lx).
"""
lengths(domain::Domain) = (domain.Ly, domain.Lx)

"""
    differential_elements()

Return a tuple of differential elements, the grid spacing along each axis, default (dy, dx).
"""
differential_elements(domain::Domain) = (domain.dy, domain.dx)

"""
    get_points(domain::Domain)

Return a tuple of points along each axis, default (y, x).
"""
get_points(domain::Domain) = (domain.y, domain.x)

"""
    get_label(domain::Domain)

Return a tuple of labels (Symbols) for each axis, default (:y, :x).
"""
get_labels(domain::Domain) = (:y, :x)

"""
    wave_vectors(domain::Domain)

Return a tuple of wave vectors for each axis, default (ky, kx).
"""
wave_vectors(domain::Domain) = (domain.ky, domain.kx)

"""
    domain_kwargs(domain::AbstractArray)

Return the domain specific keyword arguments, depending on the type of AbstractDomain.
"""
domain_kwargs(domain::Domain) = (; real_transform=domain.real_transform,
                                 dealiased=domain.dealiased)

"""
    spectral_size(domain::AbstractDomain)

Return a tuple containing the size of the spectral coefficient Array (size in spectral space).
"""
spectral_size(domain::AbstractDomain) = size(get_bwd(domain))

"""
    spectral_length(domain::AbstractDomain)

Return the number of spectral coefficients. (length in spectral space).
"""
spectral_length(domain::AbstractDomain) = prod(spectral_size(domain))

"""
    area(domain::AbstractDomain)

Compute the area of the domain. By default use prod(lengths(domain)).
"""
area(domain::AbstractDomain) = prod(lengths(domain))

"""
    differential_area(domain::AbstractDomain)

Compute the differential area of the domain. By default use prod(differential_area(domain)).
"""
differential_area(domain::AbstractDomain) = prod(differential_elements(domain))

# Getters
get_transform_plans(domain::AbstractDomain) = domain.transforms
get_fwd(domain::AbstractDomain) = get_fwd(get_transform_plans(domain))
get_bwd(domain::AbstractDomain) = get_bwd(get_transform_plans(domain))
get_precision(domain::AbstractDomain) = domain.precision
memory_type(domain::AbstractDomain) = domain.MemoryType{domain.precision}

# Overloading
Base.size(domain::AbstractDomain) = (domain.Ny, domain.Nx)
Base.length(domain::AbstractDomain) = prod(size(domain))
Base.ndims(domain::AbstractDomain) = length(size(domain))