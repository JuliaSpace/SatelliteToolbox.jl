using Documenter
using SatelliteToolbox

makedocs(;
    modules = [SatelliteToolbox],
    format = Documenter.HTML(;
        prettyurls = !("local" in ARGS),
        canonical = "https://juliaspace.github.io/SatelliteToolbox.jl/stable/",
        size_threshold = 500 * 1024,
        size_threshold_warn = 300 * 1024,
    ),
    sitename = "Satellite Toolbox",
    authors = "Ronan Arraes Jardim Chagas",
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "ISS Observation"                     => "tutorials/iss_observation.md",
            "GPS State Vector: ITRF to GCRF"      => "tutorials/gps_itrf_to_gcrf.md",
            "Geomagnetic Field Along a LEO Orbit" => "tutorials/geomagnetic_field_leo.md",
            "Comparing Gravity Models"            => "tutorials/gravity_model_comparison.md",
            "Mean Elements Fit and Covariance"    => "tutorials/mean_elements_fit_covariance.md",
        ],
        "Library" => "lib/library.md",
    ],
)

deploydocs(; repo = "github.com/JuliaSpace/SatelliteToolbox.jl.git", target = "build")
