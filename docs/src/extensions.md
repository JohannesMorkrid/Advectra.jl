# Extensions

## Submodules

Advectra.jl attempts to supports the use of all `AbstractArray` types, but can only confirm
that the following third-party types are supported:

- `CuArray` ([CUDA.jl](https://github.com/JuliaGPU/CUDA.jl))
- `ROCArray` ([AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl))
- `ComponentArrays` ([ComponentArrays.jl](https://github.com/SciML/ComponentArrays.jl))

See the [documentation]() for how to use these Array types.

In addition the code supports the [sending of mails](@ref Advectra.send_mail) throught the [`SMTPClient`](https://github.com/aviks/SMTPClient.jl) extension:
`send_mail(subject::AbstractString; attachment="")`,
which is enabled by including the package

```julia
using SMTPClient
```