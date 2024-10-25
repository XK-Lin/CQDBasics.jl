using PkgTemplates

t = Template(
    user="XK-Lin",
    dir="/Users/lin/Desktop",
    authors="Xukun Lin",
    host="github.com",
    julia=v"1.11",
    plugins=[
        Git(; ignore=["**/.DS_Store"], manifest=false, jl=false),
        GitHubActions(),
        Documenter{GitHubActions}(),
    ]
)

t("CQDBasics.jl")