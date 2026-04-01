# Advectra.jl: Advection-Diffusion Spectral Solver

_Advectra_ is a Two-Dimensional Bi-Spectral Advection Diffusion Solver package written in Julia (available on [GitHub](https://github.com/JohannesMorkrid/Advectra.jl)). 

## Overview
The package solves generic partial differential-equations (PDEs) of the form:

```math
\frac{\partial u}{\partial t} = \mathcal{L}(u, p, t) + \mathcal{N}(u, p, t),
```

where $u$ is the state at time $t$, $p$ are additional parameters, $\mathcal{L}$ is a linear operator usually associated with a diffusion process and $\mathcal{N}$ is a non-linear operator associated with the advective terms.

The package attempts to make solving the differential equations in spectral space as
trivial as possible through the use of [`SpectralOperators`](SpectralOperators.md) so that the user does not
need to translate the equations to their spectral counterpart themselves. In addition,
the package support solving multiple nested PDEs simultanously. While the code is
specialized towards solving plasma fluid equations, it is also well suited for generic
advection diffusion problems.

The code features:

- A bi-periodic [`Domain`](Domain.md)
- [`SpectralOperators`](SpectralOperators.md) to compute spatial derivatives in spectral space
- Mixed Stiffly-Stable ([`MSS`]()) time integrators; up to third order
- [`HDF5`](https://github.com/JuliaIO/HDF5.jl) data output for binary format storage with [`Blosc`](https://github.com/JuliaIO/HDF5.jl/tree/master/filters/H5Zblosc) compression
- Pseudospectral methods for non-linear terms using FFTs ([`FFTW`](https://github.com/JuliaMath/FFTW.jl))
- 3/2-[`dealiasing`]() of non-linear operators
- [`Diagnostic`](Diagnostics.md)'s for sampling at high frequencies with minimal storage
- GPU support ([`CUDA`](https://github.com/JuliaGPU/CUDA.jl), [`AMD`](https://github.com/JuliaGPU/AMDGPU.jl))
- Easy construction of canonical [`initial conditions`]() for PDEs
- Option to [`remove modes`]() of interest

<!--## Installation

Advectra.jl will soon be installable through the Julia package manager. From the Julia REPL,
type `]` to enter the Pkg REPL mode and run:

```
pkg> add Advectra
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("Advectra")
```
-->

## Contributing

Issues and contributions through pull requests are welcome. Please consult the
[contributor guide](contributor-guide.md) before submitting a pull request.

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

Software licensed under the [MIT License](https://github.com/JohannesMorkrid/Advectra.jl/blob/main/LICENSE).