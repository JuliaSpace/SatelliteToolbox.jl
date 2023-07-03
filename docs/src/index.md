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

- [SatelliteToolboxAtmosphericModels.jl](https://github.com/JuliaSpace/SatelliteToolboxAtmosphericModels.jl): Atmospheric models.
- [SatelliteToolboxBase.jl](https://github.com/JuliaSpace/SatelliteToolboxBase.jl): Base functions and type definitions.
- [SatelliteToolboxCelestialBodies.jl](https://github.com/JuliaSpace/SatelliteToolboxCelestialBodies.jl): Celestial bodies.
- [SatelliteToolboxGeomagneticField.jl](https://github.com/JuliaSpace/SatelliteToolboxGeomagneticField.jl): Geomagnetic field models.
- [SatelliteToolboxGravityModels.jl](https://github.com/JuliaSpace/SatelliteToolboxGravityModels.jl): Gravity models.
- [SatelliteToolboxLegendre.jl](https://github.com/JuliaSpace/SatelliteToolboxLegendre.jl): Legendre associated functions and its time-derivatives.
- [SatelliteToolboxPropagators.jl](https://github.com/JuliaSpace/SatelliteToolboxPropagators.jl): Orbit propagators.
- [SatelliteToolboxSgp4.jl](https://github.com/JuliaSpace/SatelliteToolboxSgp4.jl): SGP4/SDP4 orbit propagator.
- [SatelliteToolboxTle.jl](https://github.com/JuliaSpace/SatelliteToolboxTle.jl): Creating, fetching, and parsing TLEs.
- [SatelliteToolboxTransformations.jl](https://github.com/JuliaSpace/SatelliteToolboxTransformations.jl): Transformations (reference frames, time, etc.).

## Installation

This package can be installed using:

```julia-repl
julia> using Pkg
julia> Pkg.add("SatelliteToolbox")
```

## Documentation

This page contains tutorials with examples of analyses that can be performed using the
**SatelliteToolbox.jl** ecosystem. For the documentation of the functions, please refer to
the related package documentation page.

!!! warning
    This documentation is under construction, and more tutorials will be added.
