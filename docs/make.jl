using Documenter
using Advectra
DocMeta.setdocmeta!(Advectra, :DocTestSetup, :(using Advectra); recursive=true)

makedocs(; sitename="Advectra.jl",
         authors="Johannes Mørkrid",
         modules=[Advectra],
         warnonly=[:doctest,
             :missing_docs,
             :docs_block,
             :autodocs_block,
             :cross_references],
         repo="github.com/JohannesMorkrid/Advectra.jl",
         pages=["Advectra.jl: Advection-Diffusion Spectral Solver" => "index.md",
             "Installation" => "installation.md",
             "Getting started" => "getting-started.md",
             "Examples" => "examples.md",
             "Manual" => [
                 "Domain.md",
                 "initial_conditions.md",
                 "SpectralOperators.md",
                 "Diagnostics.md",
                 "boussinesq-vs-non-boussinesq.md",
                 "SpectralODEProblem.md",
                 "output.md",
                 "time-integration.md",
                 "Solver.md"
             ], # Examples
             "Extensions" => "extensions.md",
             "Contributing" => "contributor-guide.md",
             "API" => "API.md"])

deploydocs(; repo="github.com/JohannesMorkrid/Advectra.jl")