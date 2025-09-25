using AstrochemicalYields
using Documenter

DocMeta.setdocmeta!(AstrochemicalYields, :DocTestSetup, :(using AstrochemicalYields); recursive=true)

makedocs(;
    modules=[AstrochemicalYields],
    authors="cgarling <chris.t.garling@gmail.com> and contributors",
    sitename="AstrochemicalYields.jl",
    format=Documenter.HTML(;
        canonical="https://cgarling.github.io/AstrochemicalYields.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    doctest=false,
    linkcheck=true,
    warnonly=[:missing_docs, :linkcheck],
)

deploydocs(;
    repo="github.com/cgarling/AstrochemicalYields.jl",
    devbranch="main",
)
