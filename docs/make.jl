using FiniteVolumeMethod1D
using Documenter

DocMeta.setdocmeta!(FiniteVolumeMethod1D, :DocTestSetup, :(using FiniteVolumeMethod1D); recursive=true)

makedocs(;
    modules=[FiniteVolumeMethod1D],
    authors="DanielVandH <danj.vandenheuvel@gmail.com> and contributors",
    repo="https://github.com/DanielVandH/FiniteVolumeMethod1D.jl/blob/{commit}{path}#{line}",
    sitename="FiniteVolumeMethod1D.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://DanielVandH.github.io/FiniteVolumeMethod1D.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Mathematical Details" => "math.md"
    ],
)

deploydocs(;
    repo="github.com/DanielVandH/FiniteVolumeMethod1D.jl",
    devbranch="main",
)
