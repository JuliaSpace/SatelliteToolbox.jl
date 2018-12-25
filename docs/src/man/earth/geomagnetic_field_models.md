Earth geomagnetic field models
==============================

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

Currently, there is only the IGFR model in this toolbox to compute the Earth
geomagnetic field.

## IGRF v12

There is a native Julia implementation of the [International Geomagnetic
Reference Field (IGRF) v12](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html). This
can be accessed by two functions:

* `igrf12syn`: This is the native Julia implementation of the original FORTRAN
  source-code with **the same** input parameters.
* `igrf12`: An independent (more readable) implementation of the IGRF model.
  However, it is not as fast as `igrf12syn` yet (~20% slower).

The `igrf12` function has the following signature:

```julia
function igrf12(date::Number, r::Number, λ::Number, Ω::Number, T; show_warns = true)
```

It computes the geomagnetic field vector [nT] at the date `date` [Year A.D.] and
position (`r`, `λ`, `Ω`).

The position representation is defined by `T`. If `T` is `Val{:geocentric}`,
then the input must be **geocentric** coordinates:

1. Distance from the Earth center `r` \[m];
1. Geocentric latitude `λ` (``-\pi/2``, ``+\pi/2``) \[rad]; and
2. Geocentric longitude `Ω` (``-\pi``, +``\pi``) \[rad].

If `T` is `Val{:geodetic}`, then the input must be **geodetic**
coordinates:

1. Altitude above the reference ellipsoid `h` (WGS-84) \[m];
2. Geodetic latitude `λ` (``-\pi/2``, ``+\pi/2``) \[rad]; and
3. Geodetic longitude `Ω` (``-\pi``, ``+\pi``) \[rad].

If `T` is omitted, then it defaults to `Val{:geocentric}`.

Notice that the output vector will be represented in the same reference system
selected by the parameter `T` (geocentric or geodetic). The Y-axis of the output
reference system always points East. In case of **geocentric coordinates**, the
Z-axis points toward the center of Earth and the X-axis completes a right-handed
coordinate system. In case of **geodetic coordinates**, the X-axis is tangent to
the ellipsoid at the selected location and points toward North, whereas the
Z-axis completes a right-hand coordinate system.

If the keyword `show_warns` is `true` (default), then warnings will be printed
to STDOUT.

!!! note

    The IGRF v12 implemented here can be used to compute the geomagnetic field
    from 1900 up to 2025. Notice, however, that for dates after 2020 the
    accuracy is reduced.

```jldoctest
julia> igrf12(2017.12313, 640e3, 50*pi/180, 25*pi/180, Val{:geodetic})
3-element StaticArrays.SArray{Tuple{3},Float64,1,3}:
 15374.385312809889
  1267.9325221604724
 34168.28074619655

julia> igrf12(2017.12313, 6371e3+640e3, 50*pi/180, 25*pi/180, Val{:geocentric})
3-element StaticArrays.SArray{Tuple{3},Float64,1,3}:
 15174.122905727732
  1262.6765496083972
 34210.301848156094
```

```julia
julia> igrf12(2022, 6371e3+640e3, 50*pi/180, 25*pi/180)
┌ Warning: The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2020.
└ @ SatelliteToolbox ~/.julia/dev/SatelliteToolbox/src/earth/geomagnetic_field_models/igrf.jl:99
3-element StaticArrays.SArray{Tuple{3},Float64,1,3}:
 15167.758261587729
  1423.7811941605075
 34342.17638944679

julia> igrf12(2022, 6371e3+640e3, 50*pi/180, 25*pi/180; show_warns=false)
3-element StaticArrays.SArray{Tuple{3},Float64,1,3}:
 15167.758261587729
  1423.7811941605075
 34342.17638944679
```
