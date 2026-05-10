# Examples

## Example usage

Say you want to evolve the following system of coupled partial differential-equations (model from [Garcia et al.](https://doi.org/10.1063/1.2044487)):

```math 
\frac{\partial n}{\partial t} + \{\phi, n\} = \nu\nabla^2 n
```

```math 
\frac{\partial\Omega}{\partial t} + \{\phi,\Omega\} + \frac{\partial n}{\partial y} = \mu\nabla^2 \Omega
```

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

## Example files
In addition there are a lot of example files that solves different models used in the field of plasma physics

### Resistive drift-waves
The resisitve drift-wave turbulence is described by the Hasegawa-Wakatani model

```math 
    \frac{\partial n}{\partial t} + \{\phi, n\} + \kappa\frac{\partial\phi}{\partial y} = \alpha(\phi-n) + D_n\nabla^2_\perp n 
```

```math
    \frac{\partial\Omega}{\partial t} + \{\phi,\Omega\} = \alpha(\phi-n) + D_\Omega\nabla^2_\perp\Omega
```

where $D_n$ and $D_\Omega$ may include higher order damping operators.

### Sheath-interchange instabilities
Sheat-interchange instabilities are described by the following equations

```math
\frac{\partial n}{\partial t} + \{\phi, n\} - gn\frac{\partial\phi}{\partial y} + g\frac{\partial n}{\partial y} = D_n\nabla^2_\perp n - \sigma_nn\exp(\Lambda-\phi) + S_n
```

```math
\frac{\partial\Omega}{\partial t} + \{\phi,\Omega\} + g\frac{\partial\ln(n)}{\partial y} = D_\Omega\nabla^2_\perp\Omega + \sigma_\Omega[1-\exp(\Lambda-\phi)]
```