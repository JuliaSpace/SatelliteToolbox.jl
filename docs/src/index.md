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
of the [Forplan Satellite Simulator](https://journals.sagepub.com/doi/abs/10.1177/1063293X18804006).

!!! warning

    This documentation is under construction. However, the functions are
    extensively documented using the Julia documentation system, which can be
    accessed by typing `?` followed by the function name in REPL.

## Requirements

* Julia >= 1.0
* Crayons >= 4.0.0
* Interpolations >= 0.12
* OptionalData >= 0.3
* Parameters >= 0.12
* PolynomialRoots >= 0.2
* ReferenceFrameRotations >= 0.5
* RemoteFiles >= 0.3
* StaticArrays >= 0.12

## Installation

This package can be installed using:

```julia-repl
julia> using Pkg
julia> Pkg.add("SatelliteToolbox")
```
