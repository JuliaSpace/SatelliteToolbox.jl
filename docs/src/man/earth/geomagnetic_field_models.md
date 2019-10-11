Earth geomagnetic field models
==============================

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

Currently, there is two models in this toolbox to compute the Earth geomagnetic
field: the IGFR model and the simplified dipole model.

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

## Simplified dipole model

This model assumes that the Earth geomagnetic field is a perfect dipole. The
approximation is not good, but it can be sufficient for some analysis, such as
those carried out at the Pre-Phase A of a space mission, when the uncertainties
are very high.

The functions that can be used to compute the Earth geomagnetic field using this
simplified model are:

```julia
    function geomag_dipole(r_e::AbstractVector, pole_lat::Number, pole_lon::Number, m::Number)
    function geomag_dipole(r_e::AbstractVector, year::Number = 2019)
```

In the first, the geomagnetic field \[nT] is computed using the simplified
dipole model at position `r_e` (ECEF reference frame). This function considers
that the latitude of the South magnetic pole (which lies in the North
hemisphere) is `pole_lat` [rad] and the longitude is `pole_lon` [rad].
Furthermore, the dipole moment is considered to be `m` [A.m²].

In the second, the geomagnetic field \[nT] is computed using the simplified
dipole model at position `r_e` (ECEF reference frame). This function uses the
year `year` to obtain the position of the South magnetic pole (which lies in the
North hemisphere) and the dipole moment. If `year` is omitted, then it will be
considered as 2019.

!!! note

    In both functions, the output vector will be represented in the ECEF
    reference frame.

```jldoctest
julia> r_e = [0;0;R0+200e3];

julia> geomag_dipole(r_e)
3-element Array{Float64,1}:
   1286.02428617178
  -4232.804339060698
 -53444.68086319672

julia> geomag_dipole(r_e, 1986)
3-element Array{Float64,1}:
   1715.2656071053527
  -4964.598060841779
 -54246.30480714958

julia> r_e = [R0+200e3;0;0];

julia> geomag_dipole(r_e)
3-element Array{Float64,1}:
 -2572.04857234356
 -4232.804339060698
 26722.34043159836

julia> geomag_dipole(r_e, 1986)
3-element Array{Float64,1}:
 -3430.5312142107055
 -4964.598060841779
 27123.15240357479
```
