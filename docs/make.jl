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
             "Installation" => "Installation.md",
             "Getting started" => "getting-started.md",
             "Examples" => "examples.md",
             "Implementation details" => [
                 "Domain.md",
                 "SpectralOperators.md",
                 "Diagnostics.md",
                 "boussinesq-vs-non-boussinesq.md",
                 "time-integration.md",
             ], # Examples
             "Extensions" => "extensions.md",
             "Contributing" => "contributor-guide.md",
             "API" => "API.md",
             )

deploydocs(; repo="github.com/JohannesMorkrid/Advectra.jl")