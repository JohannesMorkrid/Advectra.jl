# Spectral Operators

## Spatial Derivatives
One of the benefits of spectral codes is that computing spatial derivatives becomes trivial, as Linear-Operators becomes a multiplication, in Fourier space:
```math
    \mathcal{F}\left\{\frac{\partial u}{\partial x}\right\} = (ik_x)\hat{u}
```

**Available Operators:**
```julia
diff_x
diff_y
diff_xx
diff_yy
laplacian
grad_dot_grad
```
It is also possible to raise the derivative to power by using the `order` kwarg.

## Non-Linear Operators

On the contrary to Linear-Operators, Non-Linear Operators becomes a bit more tricky to compute, and are solved [Pseudospectrally](https://en.wikipedia.org/wiki/Pseudo-spectral_method)
```math
    \mathcal{F}\left\{\mathcal{N}[u]\right\} = ...
```

```julia
quadratic_term
poisson_bracket
spectral_exp
spectral_expm1
spectral_log
reciprocal
```

## Laplacian Inversion

```julia
solve_phi
``` 

- `SolvePhiSimplified` (`boussinesq` = `false`)
- `SolvePhiNonBoussinesq` (`boussinesq` = `true`, `relaxation` = `false`)
- `SolvePhiRelaxation` (`boussinesq` = `false`, `relaxation` = `true`)

Follows Angus.

## Source Terms

```julia
spectral_constant
source
``` 

These operators are still in development.