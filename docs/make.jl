using Documenter
using SatelliteToolbox

makedocs(
    format = :html,
    modules = [SatelliteToolbox],
    sitename = "Satellite Toolbox",
    authors = "Ronan Arraes Jardim Chagas",
    pages = [
        "Home" => "index.md",
        "Earth Atmospheric Models" => "man/earth/atmospheric_models.md",
        "Library" => "lib/library.md",
    ],
    html_prettyurls = !("local" in ARGS),
)

deploydocs(
    repo = "github.com/JuliaSpace/SatelliteToolbox.jl.git",
    julia = "1.0",
    target = "build",
    deps = nothing,
    make = nothing,
)
