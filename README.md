# Advectra.jl

[![Build Status](https://github.com/JohannesMorkrid/Advectra.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JohannesMorkrid/Advectra.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JohannesMorkrid.github.io/Advectra.jl/dev)

_Advectra_ is a Two-Dimensional Bi-Spectral Advection Diffusion Solver written in Julia. The code solves generic partial differential-equations (PDEs) of the form:

$$ \frac{\partial u}{\partial t} = \mathcal{L}(u, p, t) + \mathcal{N}(u, p, t),$$

where $u$ is the state at time $t$, $p$ are additional parameters, $\mathcal{L}$ is a linear operator usually associated with a diffusion process and $\mathcal{N}$ is a non-linear operator for the advective terms.

The package attempts to make solving the differential equations in spectral space as
trivial as possible through the use of [`SpectralOperators`]() so that the user does not
need to translate the equations to their spectral counterpart themselves. In addition,
the package support solving multiple nested PDEs simultanously. While the code is
specialized towards solving plasma fluid equations, it is also well suited for generic
advection diffusion problems.

The code features:

- A bi-periodic [`Domain`]()
- [`SpectralOperators`]() to compute spatial derivatives in spectral space
- Mixed Stiffly-Stable ([`MSS`]()) time integrators; up to third order
- [`HDF5`](https://github.com/JuliaIO/HDF5.jl) data output for binary format storage with [`Blosc`](https://github.com/JuliaIO/HDF5.jl/tree/master/filters/H5Zblosc) compression
- Pseudospectral methods for non-linear terms using FFTs ([`FFTW`](https://github.com/JuliaMath/FFTW.jl))
- 3/2-[`dealiasing`]() of non-linear operators
- [`Diagnostic`]()'s for sampling at high frequencies with minimal storage
- GPU support ([`CUDA`](https://github.com/JuliaGPU/CUDA.jl), [`AMD`](https://github.com/JuliaGPU/AMDGPU.jl))
- Easy construction of canonical [`initial conditions`]() for PDEs
- Option to [`remove modes`]() of interest

## Installation

Advectra.jl will soon be installable through the Julia package manager. From the Julia REPL,
type `]` to enter the Pkg REPL mode and run:

```
pkg> add Advectra
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("Advectra")
```

## Example usage

Say you want to evolve the following system of coupled partial differential-equations (model from [Garcia et al.](https://doi.org/10.1063/1.2044487)):

$$ \frac{\partial n}{\partial t} + \\{\phi, n\\} = \nu\nabla^2 n $$

$$ \frac{\partial\Omega}{\partial t} + \\{\phi,\Omega\\} + \frac{\partial n}{\partial y} = \mu\nabla^2 \Omega $$

where $n$ is the density field, $\Omega = \nabla^2\phi$ is the voriticy field, $\phi$ is
the potential field, $\\{f, g\\} = \frac{\partial f}{\partial x}\frac{\partial g}{\partial y} - \frac{\partial f}{\partial y}\frac{\partial g}{\partial x}$ denotes the non-linear [Poisson bracket](https://en.wikipedia.org/wiki/Poisson_bracket#Definition_in_canonical_coordinates) operator and $\nu$ and $\mu$ are damping coefficients.

The diffusive terms lead to the following `Linear` operator:

```julia
function Linear(du, u, operators, p, t)
    @unpack ν, μ = p
    n, Ω = eachslice(u; dims=3)
    dn, dΩ = eachslice(du; dims=3)
    @unpack laplacian = operators

    dn .= ν .* laplacian(n)
    dΩ .= μ .* laplacian(Ω)
end
```

where most of the function is just unpacking, while the actual computations happen at the
last two lines. To compute the Laplacian operator ($\nabla^2$) it is as trivial as calling
the `laplacian` method for each field. Note the use of [broadcasting](https://docs.julialang.org/en/v1/manual/arrays/#Broadcasting) for writing [in-place](https://docs.sciml.ai/DiffEqDocs/stable/basics/problem/#In-place-vs-Out-of-Place-Function-Definition-Forms).

Similarly, the advective terms lead to the following `NonLinear` operator:

```julia
function NonLinear(du, u, operators, p, t)
    n, Ω = eachslice(u; dims=3)
    dn, dΩ = eachslice(du; dims=3)
    @unpack diff_y, poisson_bracket, solve_phi = operators
    ϕ = solve_phi(n, Ω)
    dn .= poisson_bracket(n, ϕ)
    dΩ .= poisson_bracket(Ω, ϕ) - diff_y(n)
end
```

which is a bit more complicated, as the laplacian has to be inversed using the `solve_phi`
method. Other than that the derivative in y is computed using `diff_y` and the Poisson
bracket is computed using the `poisson_bracket`method.

The right hand sides needs to be collected in a [`SpectralODEProblem`]():
```julia
prob = SpectralODEProblem(Linear, NonLinear, u0, domain, time_span; p=parameters, dt=2.5e-3, diagnostics=diagnostics)
```
alongside the initial state `u0`, the simulation [`Domain`](), the `time_span` to integrate 
over, the `parameters` of the system and a `Vector` of [`Diagnostic`]()s to be performed.

To solve the system, use the following method:
```julia
sol = spectral_solve(prob, MSS3(), output;)
```
where `MMS3` is the time integration [`scheme`]() and `output` is a constructed [`Output`](). 

See the [example file](examples/Garcia%202005%20Pop/Garcia%202005%20PoP.jl) for more details.

### Results:

![Alt Text](assets/blob.gif)

## Submodules

Advectra.jl attempts to supports the use of all `AbstractArray` types, but can only confirm
that the following third-party types are supported:

- `CuArray` ([CUDA.jl](https://github.com/JuliaGPU/CUDA.jl))
- `ROCArray` ([AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl))
- `ComponentArrays` ([ComponentArrays.jl](https://github.com/SciML/ComponentArrays.jl))

See the [documentation]() for how to use these Array types.

In addition the code supports the [sending of mails]() throught the [`SMTPClient`](https://github.com/aviks/SMTPClient.jl) extension:
`send_mail(subject::AbstractString; attachment="")`,
which is enabled by including the package

```julia
using SMTPClient
```

## Contributing

Issues and contributions through pull requests are welcome. Please consult the
[contributor guide]() before submitting a pull request.

## Citation

If you use Advectra.jl in your research, teaching, or other activities, please cite this
repository using the following:

```bibtex
@software{Advectra,
  author       = {Johannes Mørkrid and contributors},
  title        = {Advectra},
  year         = {2026},
  url          = {https://github.com/JohannesMorkrid/Advectra.jl},
  version      = {0.1.0},
  license      = {MIT},
  note         = {GitHub repository},
}
```

## Copyright and license

Copyright (c) 2026 Johannes Mørkrid (johannes.e.morkrid@uit.no) and contributors for Advectra.jl

Software licensed under the [MIT License](LICENSE).


