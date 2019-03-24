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

!!! warning

    This documentation is under construction. However, the functions are
    extensively documented using the Julia documentation system, which can be
    accessed by typing `?` followed by the function name in REPL.

## Requirements

* Julia >= 1.0
* Crayons >= 4.0.0
* Interpolations >= 0.11.2
* Parameters >= 0.10.3
* OptionalData >= 0.2.1
* PolynomialRoots >= 0.2.0
* ReferenceFrameRotations >= 0.5.1
* RemoteFiles >= 0.2.1
* StaticArrays >= 0.10.3

## Installation

This package can be installed using:

```julia-repl
julia> using Pkg
julia> Pkg.add("SatelliteToolbox")
```
