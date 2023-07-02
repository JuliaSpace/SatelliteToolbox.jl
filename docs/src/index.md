SatelliteToolbox.jl
===================

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

The **SatelliteToolbox.jl** contains a set of packages with functions to perform analysis
and build simulations related to satellites. It is used on a daily basis on projects at the
[Brazilian National Institute for Space Research (INPE)](http://www.gov.br/inpe).

The set of packages are listed bellow. All of them are loaded and reexported in this one.

| Package Name                                                                   | Description                                            | Build Status                                          | Coverage                                               |
|--------------------------------------------------------------------------------|--------------------------------------------------------|-------------------------------------------------------|--------------------------------------------------------|
| [SatelliteToolboxAtmosphericModels.jl][SatelliteToolboxAtmosphericModels-link] | Atmospheric models                                     | ![Build status][SatelliteToolboxAtmosphericModels-ci] | ![Build status][SatelliteToolboxAtmosphericModels-cov] |
| [SatelliteToolboxBase.jl][SatelliteToolboxBase-link]                           | Base functions and type definitions                    | ![Build status][SatelliteToolboxBase-ci]              | ![Coverace][SatelliteToolboxBase-cov]                  |
| [SatelliteToolboxCelestialBodies.jl][SatelliteToolboxCelestialBodies-link]     | Celestial bodies                                       | ![Build status][SatelliteToolboxCelestialBodies-ci]   | ![Coverage][SatelliteToolboxCelestialBodies-cov]       |
| [SatelliteToolboxGeomagneticField.jl][SatelliteToolboxGeomagneticField-link]   | Geomagnetic field models                               | ![Build status][SatelliteToolboxGeomagneticField-ci]  | ![Coverage][SatelliteToolboxGeomagneticField-cov]      |
| [SatelliteToolboxGravityModels.jl][SatelliteToolboxGravityModels-link]         | Gravity models                                         | ![Build status][SatelliteToolboxGravityModels-ci]     | ![Coverage][SatelliteToolboxGravityModels-cov]         |
| [SatelliteToolboxLegendre.jl][SatelliteToolboxLegendre-link]                   | Legendre associated functions and its time-derivatives | ![Build status][SatelliteToolboxLegendre-ci]          | ![Coverage][SatelliteToolboxLegendre-cov]              |
| [SatelliteToolboxPropagators.jl][SatelliteToolboxPropagators-link]             | Orbit propagators                                      | ![Build status][SatelliteToolboxPropagators-ci]       | ![Coverage][SatelliteToolboxPropagators-cov]           |
| [SatelliteToolboxSgp4.jl][SatelliteToolboxSgp4-link]                           | SGP4/SDP4 orbit propagator                             | ![Build status][SatelliteToolboxSgp4-ci]              | ![Coverage][SatelliteToolboxSgp4-cov]                  |
| [SatelliteToolboxTle.jl][SatelliteToolboxTle-link]                             | Creating, fetching, and parsing TLEs                   | ![Build status][SatelliteToolboxTle-ci]               | ![Coverage][SatelliteToolboxTle-cov]                   |
| [SatelliteToolboxTransformations.jl][SatelliteToolboxTransformations-link]     | Transformations (reference frames, time, etc.)         | ![Build status][SatelliteToolboxTransformations-ci]   | ![Coverage][SatelliteToolboxTransformations-cov]       |

## Installation

This package can be installed using:

```julia-repl
julia> using Pkg
julia> Pkg.add("SatelliteToolbox")
```
