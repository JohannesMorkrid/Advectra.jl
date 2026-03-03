# Advectra

[![Build Status](https://github.com/JohannesMorkrid/Advectra/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JohannesMorkrid/Advectra/actions/workflows/CI.yml?query=branch%3Amain)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://JohannesMorkrid.github.io/Advectra.jl/dev)

Advectra is a Bi-Spectral Advection Diffusion Solver written in Julia. The code solves generic differential-equations of the form:
$$ \frac{\partial u}{\partial t} = \mathcal{L}(u, p, t) + \mathcal{N}(u, p, t),$$
where $u$ is the state at time $t$, $p$ are additional parameters, $\mathcal{L}$ is a linear operator usually associated with a diffusion process and $\mathcal{N+}$ is a non-linear operator for the advective terms.

The code attempts to write the differ

And supports multiple fields

The code atempts to be modular and generalizable to be able to solve other spectral problems. 

The code features:
* Biperiodic domain (perpendicular to $\textbf{B}$)
* Fast Fourier Transform for spatial derivatives (FFTW)
* Third order Stiffly-Stable time integrator
* HDF data output for binary format storage with Blosc compression
* 2/3 Antialiasing on quadratic terms and non-linear functions
* Diagnostic modules
* GPU support (CUDA, AMD)

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

## How to use

to simulate sheath-interchange instability and resistive drift-wave turbulence in magnetized plasmas.
The resisitve drift-wave turbulence is described by the Hasegawa-Wakatani model

$$ \frac{\partial n}{\partial t} + \\{\phi, n\\} + \kappa\frac{\partial\phi}{\partial y} = \alpha(\phi-n) + D_n\nabla^2_\perp n $$

$$ \frac{\partial\Omega}{\partial t} + \\{\phi,\Omega\\} = \alpha(\phi-n) + D_\Omega\nabla^2_\perp\Omega $$

where $D_n$ and $D_\Omega$ may include higher order damping operators,while sheat-interchange instabilities are described by the following equations

Results:
![Alt Text](assets/blob.gif)

## Submodules

Advectra.jl attempts to supports the use of all `AbstractArray` types, but can only confirm
that the following third-party types are supported:
* `CuArray` ([CUDA.jl](https://github.com/JuliaGPU/CUDA.jl))
* `ROCArray` ([AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl))
* `ComponentArrays` ([ComponentArrays.jl](https://github.com/SciML/ComponentArrays.jl))

See the documentation for how to use these Array types.

In addition the code supports the sending of mails throught the [`SMTPClient`](https://github.com/aviks/SMTPClient.jl) extensions 
function `send_mail(subject::AbstractString; attachment="")`.

## Contributing

Issues and contributions through pull requests are welcome. Please consult the 
[contributor guide]() before submitting a pull request.

Things want to add in future versions:
* Operators, remediscent of SciMLOperators
* Rosenbrock-Euler method for first step

## Citation

If you use Advectra.jl in research, teaching, or other activities, please cite this 
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