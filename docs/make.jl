using Documenter, Dangerous

makedocs(
    sitename = "Do not enter!",
    format = Documenter.HTML(
        assets = [
            "assets/favicon.ico",
        ]
    ),
)

deploydocs(
    repo = "github.com/penelopeysm/Dangerous.jl.git",
)
