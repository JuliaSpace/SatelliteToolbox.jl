<p align="center">
  <img src="./docs/src/assets/logo.png" width="150" title="SatelliteToolboxTransformations.jl"><br>
</p>

# SatelliteToolbox.jl

[![Build status](https://github.com/JuliaSpace/SatelliteToolbox.jl/workflows/CI/badge.svg)](https://github.com/JuliaSpace/SatelliteToolbox.jl/actions)
[![codecov](https://codecov.io/gh/JuliaSpace/SatelliteToolbox.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaSpace/SatelliteToolbox.jl)
[![docs-stable](https://img.shields.io/badge/docs-stable-blue.svg)][docs-stable-url]
[![docs-dev](https://img.shields.io/badge/docs-dev-blue.svg)][docs-dev-url]
[![DOI](https://zenodo.org/badge/133264084.svg)](https://zenodo.org/doi/10.5281/zenodo.10396045)

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

## Documentation

Please, see each package for the related documentation of the functions.

[docs-dev-url]: https://juliaspace.github.io/SatelliteToolbox.jl/dev
[docs-stable-url]: https://juliaspace.github.io/SatelliteToolbox.jl/stable
[SatelliteToolboxAtmosphericModels-link]: https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl
[SatelliteToolboxAtmosphericModels-cov]: https://codecov.io/gh/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/branch/main/graph/badge.svg?token=oQOhGnQmdG
[SatelliteToolboxAtmosphericModels-ci]: https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl/workflows/CI/badge.svg
[SatelliteToolboxBase-link]: https://github.com/JuliaSpace/SatelliteToolboxBase.jl
[SatelliteToolboxBase-cov]: https://codecov.io/gh/JuliaSpace/SatelliteToolboxBase.jl/branch/main/graph/badge.svg?token=YADU7IB8CT
[SatelliteToolboxBase-ci]: https://github.com/JuliaSpace/SatelliteToolboxBase.jl/workflows/CI/badge.svg
[SatelliteToolboxCelestialBodies-link]: https://github.com/JuliaSpace/SatelliteToolboxCelestialBodies.jl
[SatelliteToolboxCelestialBodies-cov]: https://codecov.io/gh/JuliaSpace/SatelliteToolboxCelestialBodies.jl/branch/main/graph/badge.svg?token=CONQMSI4JD
[SatelliteToolboxCelestialBodies-ci]: https://github.com/JuliaSpace/SatelliteToolboxCelestialBodies.jl/workflows/CI/badge.svg
[SatelliteToolboxGeomagneticField-link]: https://github.com/JuliaSpace/SatelliteToolboxGeomagneticField.jl
[SatelliteToolboxGeomagneticField-cov]: https://codecov.io/gh/JuliaSpace/SatelliteToolboxGeomagneticField.jl/branch/main/graph/badge.svg?token=HW2Y9NA0L5
[SatelliteToolboxGeomagneticField-ci]: https://github.com/JuliaSpace/SatelliteToolboxGeomagneticField.jl/workflows/CI/badge.svg
[SatelliteToolboxGravityModels-link]: https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl
[SatelliteToolboxGravityModels-cov]: https://codecov.io/gh/JuliaSpace/SatelliteToolboxGravityModels.jl/branch/main/graph/badge.svg?token=47G4OLV6PD
[SatelliteToolboxGravityModels-ci]: https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl/workflows/CI/badge.svg
[SatelliteToolboxLegendre-link]: https://github.com/JuliaSpace/SatelliteToolboxLegendre.jl
[SatelliteToolboxLegendre-cov]: https://codecov.io/gh/JuliaSpace/SatelliteToolboxLegendre.jl/branch/main/graph/badge.svg?token=AUE8ZZ5IXJ
[SatelliteToolboxLegendre-ci]: https://github.com/JuliaSpace/SatelliteToolboxLegendre.jl/workflows/CI/badge.svg
[SatelliteToolboxPropagators-link]: https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl
[SatelliteToolboxPropagators-cov]: https://codecov.io/gh/JuliaSpace/SatelliteToolboxPropagators.jl/branch/main/graph/badge.svg?token=WSVR7QYKOD
[SatelliteToolboxPropagators-ci]: https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl/workflows/CI/badge.svg
[SatelliteToolboxSgp4-link]: https://github.com/JuliaSpace/SatelliteToolboxSgp4.jl
[SatelliteToolboxSgp4-cov]: https://codecov.io/gh/JuliaSpace/SatelliteToolboxSgp4.jl/branch/main/graph/badge.svg?token=480UYDX6H5
[SatelliteToolboxSgp4-ci]: https://github.com/JuliaSpace/SatelliteToolboxSgp4.jl/workflows/CI/badge.svg
[SatelliteToolboxTle-link]: https://github.com/JuliaSpace/SatelliteToolboxTle.jl
[SatelliteToolboxTle-cov]: https://codecov.io/gh/JuliaSpace/SatelliteToolboxTle.jl/branch/main/graph/badge.svg?token=SPIKBIN3ES
[SatelliteToolboxTle-ci]: https://github.com/JuliaSpace/SatelliteToolboxTle.jl/workflows/CI/badge.svg
[SatelliteToolboxTransformations-link]: https://github.com/JuliaSpace/SatelliteToolboxTransformations.jl
[SatelliteToolboxTransformations-cov]: https://codecov.io/gh/JuliaSpace/SatelliteToolboxTransformations.jl/branch/main/graph/badge.svg?token=SH31IN1JXM
[SatelliteToolboxTransformations-ci]: https://github.com/JuliaSpace/SatelliteToolboxTransformations.jl/workflows/CI/badge.svg
