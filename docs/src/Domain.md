# Domain
Currently the only supported domain is the default Two-dimensional Bi-Spectral domain, and as such it is what is constructed when calling the `Domain` constructor. To construct a unit-square domain it is as simple as calling `Domain(N)`, with the kwarg `L` controlling the length of the square domain. 

The most general constructor takes the form:
```julia
    Domain(Nx::Integer, Ny::Integer; Lx::Number=1, Ly::Number=1, 
            MemoryType::Type{<:AbstractArray}=Array, precision::DataType=Float64, 
            real_transform::Bool=true, dealiased::Bool=true, x0::Number=-Lx / 2, y0::Number=-Ly / 2)
```
where $L_i$ and $N_i$ is the length and size respectively of the $i$-th dimension.

In addition to the parsed args the `Domain` struct also stores:

- `dx`: spatial resolution/grid spacing, ```dx = 2Lx÷(Nx-1)```, in x-direction.
- `dy`: spatial resolution/grid spacing, ```dy = 2Ly÷(Nx-1)```, in y-direction.
- `x`: spatial position of each point along the x-axis. Uniformly distributed by default.
- `y`: spatial position of each point along the y-axis. Uniformly distributed by default.
- `kx`: wave vector components along x-axis.
- `ky`: wave vector components along y-axis.
- `transforms`: collection of fwd and bwd transform methods.

The `Domain` struct is the backbone for allocating `Array`s and computing all the [`SpectralOperators`](SpectralOperators.jl) as well transforming between physical and spectral space. 

The remainder of the kwargs are explained in their own subsection.

## Real transform
Assuming the fields are `Real`, one can utilize the conjugate symmetry of the Fourier coefficients
```math
\hat{u}_{-k} = \hat{u}^*_{k},
```
hence only half the coefficients are needed to evolve the equations, which is exploited by the real Fast Fourier Transform (`rFFT`). The use of `rFFT` can be toggled using the `real_transform` kwarg (default: `true`), otherwise the standard Fast Fourier Transform (`FFT`) methods are used.

## Dealiasing
Transforming a non-linear operator acting on a field $u$ to spectral space results in a convolution between the spectral modes. For instance, the quadratic operator $u^2$ in Fourier space becomes a convolution $\hat{u}*\hat{u}$. While for a general non-linear operator the convolutions becomes apparent when using the Taylor approximation of the operator and then Fourier transforming. 

In the case of a Discrete Fourier Transform, the convolution can be written as a discretized sum
```math
    \hat{u}*\hat{u} = \sum_{k'} \hat{u}_{k} \hat{u}_{k-k'}
```
where we sum over all pairs of wave-vectors $k$ and $k-k'$ resolved on the domain. When the Discrete Fourier Transform is truncated the periodic nature of the domain results in a combination of low- and high-frequency un-physically influencing the $(k-k') \mod N = k_a$ wave-vectors. This is similar to what happen in integer overflow, where very large numbers "rolls-over and becomes negative".

To preserve all information and mitigate the aliasing, Orszag's 3/2-rule for dealiasing has been implemented, which utilizes zero-padding of the spectral modes such that the overflowing wave-numbers have zero amplitude. The optimal zero-padded spectral array takes the form:
```math
\tilde{u}_{ji} = \begin{cases} \hat{u}_{ji}, & -N_x/2 \leq i \leq N_x/2 - 1, & -N_y/2 \leq j \leq N_y/2 - 1, \\ 0, & \text{otherwise}. \end{cases}
```
with `size(ũ) = (3N_y÷2, 3N_x÷2)`. The use of Orszag's 3/2 rule can be toggled using the `dealiased` kwarg (default: `true`), otherwise aliasing will occur.

## Offset
The domain is usually centered at zero, however the position can be offset by specifying the location of the lower left point $(x_0,y_0)$ using the respective `x0` and `y0` kwargs.

## 1D-domain
The code technically supports one dimensional Arrays by setting either $Nx=1$ or $Ny=1$.

## Helper functions

These helper functions come in handy when constructing `Arrays` and for use in diagnostics.

```@docs
Advectra.size(::Advectra.AbstractDomain)
Advectra.spectral_size
Advectra.ndims(::Advectra.AbstractDomain)
Advectra.length(::Advectra.AbstractDomain)
Advectra.spectral_length
Advectra.lengths
Advectra.area
Advectra.get_transform_plans
Advectra.get_fwd
Advectra.get_bwd
Advectra.memory_type
Advectra.get_precision
Advectra.wave_vectors
Advectra.differential_elements
Advectra.get_points
Advectra.get_labels
Advectra.domain_kwargs
```