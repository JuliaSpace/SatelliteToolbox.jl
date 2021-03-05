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

## IGRF

There is a native Julia implementation of the [International Geomagnetic
Reference Field (IGRF) v13](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html). This
can be accessed by two functions:

* `igrf13syn`: This is the native Julia implementation of the original FORTRAN
  source-code with **the same** input parameters.
* `igrf`: An independent (more readable) implementation of the IGRF model.
  However, it is not as fast as `igrf13syn` yet (~20% slower).

The `igrf` function has the following signature:

```julia
function igrf(date::Number, [r,h]::Number, λ::Number, Ω::Number, T[, P, dP]; show_warns = true)
```

It computes the geomagnetic field vector [nT] at the date `date` [Year A.D.] and
position (`r`, `λ`, `Ω`).

The position representation is defined by `T`. If `T` is `Val(:geocentric)`,
then the input must be **geocentric** coordinates:

1. Distance from the Earth center `r` \[m];
1. Geocentric latitude `λ` (``-\pi/2``, ``+\pi/2``) \[rad]; and
2. Geocentric longitude `Ω` (``-\pi``, +``\pi``) \[rad].

If `T` is `Val(:geodetic)`, then the input must be **geodetic**
coordinates:

1. Altitude above the reference ellipsoid `h` (WGS-84) \[m];
2. Geodetic latitude `λ` (``-\pi/2``, ``+\pi/2``) \[rad]; and
3. Geodetic longitude `Ω` (``-\pi``, ``+\pi``) \[rad].

If `T` is omitted, then it defaults to `Val(:geocentric)`.

Notice that the output vector will be represented in the same reference system
selected by the parameter `T` (geocentric or geodetic). The Y-axis of the output
reference system always points East. In case of **geocentric coordinates**, the
Z-axis points toward the center of Earth and the X-axis completes a right-handed
coordinate system. In case of **geodetic coordinates**, the X-axis is tangent to
the ellipsoid at the selected location and points toward North, whereas the
Z-axis completes a right-hand coordinate system.

The optional arguments `P` and `dP` must be two matrices with at least 14x14
real numbers. If they are present, then they will be used to store the Legendre
coefficients and their derivatives. In this case, no allocation will be
performed when computing the magnetic field. If they are not present, then 2
allocations will happen to create them.

If the keyword `show_warns` is `true` (default), then warnings will be printed
to STDOUT.

The latitude `λ` and longitude `Ω` can be passed in degrees instead of radians
by using the function `igrfd`. All the other arguments and keywords of this
function are the same as those in the function `igrf`.

!!! note

    The IGRF v13 implemented here can be used to compute the geomagnetic field
    from 1900 up to 2030. Notice, however, that for dates after 2025 the
    accuracy is reduced.

```jldoctest
julia> igrf(2017.12313, 640e3, 50*pi/180, 25*pi/180, Val(:geodetic))
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 15365.787505205592
  1274.9958640697
 34201.21820333791

julia> igrfd(2017.12313, 640e3, 50, 25, Val(:geodetic))
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 15365.787505205592
  1274.9958640697
 34201.21820333791

julia> igrf(2017.12313, 6371e3+640e3, 50*pi/180, 25*pi/180, Val(:geocentric))
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 15165.486702524944
  1269.7264334427598
 34243.04928373083

julia> igrfd(2017.12313, 6371e3+640e3, 50, 25, Val(:geocentric))
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 15165.486702524944
  1269.7264334427598
 34243.04928373083
```

```julia
julia> igrf(2026, 6371e3+640e3, 50*pi/180, 25*pi/180)
┌ Warning: The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.
└ @ SatelliteToolbox ~/.julia/dev/SatelliteToolbox/src/earth/geomagnetic_field_models/igrf/igrf.jl:103
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 15118.591511098817
  1588.129544718571
 34668.84185460438

julia> igrfd(2026, 6371e3+640e3, 50, 25)
┌ Warning: The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.
└ @ SatelliteToolbox ~/.julia/dev/SatelliteToolbox/src/earth/geomagnetic_field_models/igrf/igrf.jl:103
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 15118.591511098817
  1588.129544718571
 34668.84185460438

julia> igrf(2026, 6371e3+640e3, 50*pi/180, 25*pi/180; show_warns = false)
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 15118.591511098817
  1588.129544718571
 34668.84185460438

julia> igrfd(2026, 6371e3+640e3, 50, 25; show_warns = false)
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 15118.591511098817
  1588.129544718571
 34668.84185460438
```

!!! info

    For compatibility reasons, the v12 of the IGRF model is still available
    using the function `igrf12syn`.

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
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
   1286.02428617178
  -4232.804339060698
 -53444.68086319672

julia> geomag_dipole(r_e, 1986)
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
   1715.2656071053527
  -4964.598060841779
 -54246.30480714958

julia> r_e = [R0+200e3;0;0];

julia> geomag_dipole(r_e)
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 -2572.04857234356
 -4232.804339060698
 26722.34043159836

julia> geomag_dipole(r_e, 1986)
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 -3430.5312142107055
 -4964.598060841779
 27123.15240357479
```
