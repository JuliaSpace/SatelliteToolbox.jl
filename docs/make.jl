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
        "Earth" => Any[
            "Earth atmospheric models" => "man/earth/atmospheric_models.md",
            "Earth geomagnetic field models" => "man/earth/geomagnetic_field_models.md",
            "Space indices" => "man/earth/space_indices.md",
           ],
        "Orbit" => Any[
            "Orbit propagators" => "man/orbit/propagators.md",
            "TLE" => "man/orbit/tle.md",
        ],
        "Transformations" => Any[
            "ECEF and ECI" => "man/transformations/ecef_eci.md",
            "Geodetic and Geocentric" => "man/transformations/geodetic_geocentric.md",
        ],
        "Library" => "lib/library.md",
    ],
)

deploydocs(
    repo = "github.com/JuliaSpace/SatelliteToolbox.jl.git",
    target = "build",
)
