using CQDBasics
using Documenter

DocMeta.setdocmeta!(CQDBasics, :DocTestSetup, :(using CQDBasics); recursive=true)

makedocs(;
    modules=[CQDBasics],
    authors="Xukun Lin",
    sitename="CQDBasics.jl",
    format=Documenter.HTML(;
        canonical="https://XK-Lin.github.io/CQDBasics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Base" => "cqdbase.md",
        "Data Analysis Base" => "cqddataanalysisbase.md"
    ],
)

deploydocs(;
    repo="github.com/XK-Lin/CQDBasics.jl",
    devbranch="main",
)
