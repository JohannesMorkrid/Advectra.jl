# Style guide
The following signatures should be followed

### Order of problem related things
The first arguments should always be field/state related (`du` and `u`), then followed by 
metadata/parameters, and time should always be last. As the `Domain` is more often used it 
should be first in the case of meta data, followed by the `parameters=p`.

In-place:
```du, u, d, p, t```

Out-of-place:
```u, d, p, t```

If `prob::SpectralODEProblem` is in use:
```u, prob, t```¹

¹Because the `SpectralODEProblem` stores the domain and parameters.

### Exception to the field/state first
All arguments that are modified or related to the method goes first. For instance
```spectral_function(f::F, u::U, SC::SOC)```

### Order during construction
Should follow order of struct fields as close as possible. As for `Bool` flags, they should 
be in the order of speed-up to the code or (?). In the case of domain that amounts to:
```use_cuda, precision, real_transform, dealiased```

### Order of other arguments
Should follow order of struct construction and or order of struct fields.

## Adding TODOs
When wanting to add TODO comments, please also create a [GitHub Issue](https://github.com/JohannesMorkrid/Advectra.jl/issues)
 and add it to your TODO comment. Example

```
# TODO add suitable test [#4](https://github.com/JohannesMorkrid/Advectra.jl/issues/4)
function example(...) = ...
```

## DocStrings
Mainly follow the 11 steps outlined in the [manual](https://docs.julialang.org/en/v1/manual/documentation/#Writing-Documentation). 
However emphasis should be put on separiting `## Positional Arguments` and `## Keyword Arguments`
as done by SciML ([example](https://github.com/SciML/NeuralPDE.jl/blob/master/src/ode_solve.jl)).
Type should be included in the description instead of for the argument name.