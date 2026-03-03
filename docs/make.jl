using Advectra, Documenter
DocMeta.setdocmeta!(Advectra, :DocTestSetup, :(using Advectra);
                    recursive=true)

makedocs(; sitename="Advectra",
         authors="Johannes Mørkrid",
         modules=[Advectra],
         warnonly=[:doctest, :missing_docs])

deploydocs(; repo="github.com/JohannesMorkrid/Advectra.jl.git")