# Getting started

## Constructing a Domain

```julia
using Advectra
using CUDA

domain = Domain(Nx, Ny, Lx=2π, Ly=2π)
``` 

It is possible to toggle `dealiasing` and `real_transform`. See [Domain kwargs](Domain.md) for all available kwargs.

### Running on GPU

Deploying the code to a GPU (running sequentially with GPU speedups) is as simple as changing the `MemoryType` of the constructed `Domain`. The code uses [GPUArrays.jl](https://github.com/JuliaGPU/GPUArrays.jl) so it should in theory be able to be ran on an arbitrary GPU, however the code has been primarily developed using [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) and is known to work on AMD GPUs too.

#### CUDA

```julia
using Advectra
using CUDA

domain = Domain(Nx, Ny, Lx=2\pi, Ly=2\pi, MemoryType=CuArray)
``` 

#### AMDGPU

```julia
using Advectra
using AMDGPU

domain = Domain(Nx, Ny, Lx=2\pi, Ly=2\pi, MemoryType=ROCArray)
```

## Explain spectral and physical space