# Installation

Advectra.jl will soon be installable through the Julia package manager. From the Julia REPL,
type `]` to enter the Pkg REPL mode and run:

```
pkg> add Advectra
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("Advectra")
```

## Developing Advectra
To develop the Advectra package locally from a fresh start run the following:

```julia
julia> Pkg.develop(url="https://github.com/JohannesMorkrid/Advectra.jl")
``` 

or in Pkg REPL mode 

```
pkg> dev https://github.com/JohannesMorkrid/Advectra.jl
```

To develop from local path do

```
pkg> dev <LOCAL PATH>
```

## Julia version compatibility
Advectra.jl requires Julia v1.9 or later. The package is tested on Julia 1.9, 1.11 and latest.