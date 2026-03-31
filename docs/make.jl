using Documenter
using Advectra
DocMeta.setdocmeta!(Advectra, :DocTestSetup, :(using Advectra); recursive=true)

makedocs(; sitename="Advectra",
         authors="Johannes Mørkrid",
         modules=[Advectra],
         warnonly=[:doctest,
             :missing_docs,
             :docs_block,
             :autodocs_block,
             :cross_references])

deploydocs(; repo="github.com/JohannesMorkrid/Advectra.jl")