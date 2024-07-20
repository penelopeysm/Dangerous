using Documenter

makedocs(
    sitename = "Do not enter!",
    format = Documenter.HTML(
        assets = [
            asset("favicon.ico"),
        ]
    ),
)

deploydocs(
    repo = "github.com/penelopeysm/Dangerous.jl.git",
)
