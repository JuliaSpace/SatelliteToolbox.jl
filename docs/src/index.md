SatelliteToolbox.jl
===================

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

This package contains several functions to build simulations related with
satellites. It is used on a daily basis on projects at the [Brazilian National
Institute for Space Research (INPE)](http://www.inpe.br), and it is the engine
of the [FOrPlan Satellite Simulator](http://old.esaconferencebureau.com/docs/default-source/16c11-secesa-docs/39_chagas_presentation.pdf?sfvrsn=2).

## Requirements

* Julia >= 0.7
* HTTP >= 0.6.13
* Interpolations >= 0.8.0
* Parameters >= 0.10.1
* OptionalData >= 0.2.0
* ReferenceFrameRotations >= 0.5.0
* RemoteFiles >= 0.2.1
* StaticArrays >= 0.9.2

## Installation

This package can be installed using:

```julia-repl
julia> Pkg.update()
julia> Pkg.add("SatelliteToolbox")
```

## Manual outline

```@contents
Pages = [
    "./man/earth/atmospheric_models.md"
]
Depth = 2
```

## Library documentation

```@index
Pages = ["lib/library.md"]
```

