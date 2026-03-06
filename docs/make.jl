using Documenter
using Advectra
DocMeta.setdocmeta!(Advectra, :DocTestSetup, :(using Advectra); recursive=true)

makedocs(; sitename="Advectra",
         authors="Johannes Mørkrid",
         modules=[Advectra],
         warnonly=[:doctest, :missing_docs])

deploydocs(; repo="https://github.com/JohannesMorkrid/Advectra.jl")