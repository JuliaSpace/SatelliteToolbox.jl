using Documenter
using SatelliteToolbox

makedocs(
    modules = [SatelliteToolbox],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical = "https://juliaspace.github.io/SatelliteToolbox.jl/stable/",
    ),
    sitename = "Satellite Toolbox",
    authors = "Ronan Arraes Jardim Chagas",
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "ISS Observation" => "tutorials/iss_observation.md",
        ],
        "Library" => "lib/library.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaSpace/SatelliteToolbox.jl.git",
    target = "build",
)
