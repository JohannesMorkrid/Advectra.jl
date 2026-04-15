# Initial conditions

The recommended way to construct initial conditions is through the use of the `initial_condition` method:

```julia
ic = initial_condition(f::Function, domain::AbstractDomain)
```

The function `f` can either be a pointwise 2D function `f(x, y; kwargs...)`, broadcast
over the grid, or a domain-level function `f(domain; kwargs...)` for cases where a
pointwise evaluation is not sufficient, such as when the initial condition is constructed
in spectral space and transformed back to physical space. When using the
[Available initial conditions](@ref), the user does not need to worry about this
distinction — it is handled internally. For implementing custom domain-level initial
conditions, see [Creating a new non-broadcastable initial condition](@ref).

## Available initial conditions

```@docs
Advectra.isolated_blob
Advectra.isolated_thermal_blob
Advectra.gaussian
Advectra.log_gaussian
Advectra.random_crossphased
Advectra.gaussian_x
Advectra.gaussian_y
Advectra.sinusoidal
Advectra.sinusoidal_x
Advectra.sinusoidal_y
Advectra.exponential_x
Advectra.random_phase
Advectra.white_noise
Advectra.quadratic_y
```

## Creating a new non-broadcastable initial condition

To implement a domain-level initial condition, mark the function with the `@nobroadcast`
macro:

```julia
@nobroadcast function my_initial_condition(domain::AbstractDomain; kwargs...)
    u = randn(size(domain)...)
    u[1, :] .= 0
    return u
end
```

This tells `initial_condition` to pass the full domain to `f` rather than broadcasting
`f` pointwise over the grid.
