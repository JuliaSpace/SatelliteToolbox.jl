# ECEF and ECI

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

This package currently provides the entire IAU-76/FK5 model to transform
reference systems. The following table lists the available coordinate frames and
how they can be referenced in the functions that will be described later on.

| Reference | Type |            Coordinate frame name            |
|-----------|------|---------------------------------------------|
| `ITRF()`  | ECEF | International Terrestrial Reference Frame   |
| `PEF()`   | ECEF | Pseudo-Earth Fixed reference frame          |
| `MOD()`   | ECI  | Mean-Of-Date reference frame                |
| `TOD()`   | ECI  | True-Of-Data reference frame                |
| `GCRF()`  | ECI  | Geocentric Celestial Reference Frame (GCRF) |
| `J2000()` | ECI  | J2000 reference frame                       |
| `TEME()`  | ECI  | True Equator, Mean Equinox reference frame  |

!!! note

    ECEF stands for Earth-Centered, Earth-Fixed whereas ECI stands for
    Earth-Centered Inertial.

## EOP Data

The conversions here sometimes requires additional data related to the Earth
orientation. This information is provided by [IERS](https://www.iers.org)
(International Earth Rotation and Reference Systems Service). The
SatelliteToolbox.jl has the capability to automatically download and parse the
IERS EOP (Earth Orientation Parameters) data.

The function that will automatically download the files, store them in the file
system, and parse the data is:

```julia
function get_iers_eop(data_type::Symbol = :IAU1980; force_download = false)
```

in which:

* `data_type` is a symbol that specify what kind of data is desired (`:IAU1980`
  for IAU1980 data and `:IAU2000A` for IAU2000A data). If omitted, then it
  defaults to `:IAU1980`.
* The files are obtained on a daily-basis by the package RemoteFiles.jl. If the
  user wants to force the download, then the keyword `force_download` should be
  set to `true`.
* This function returns an instance of the structure `EOPData_IAU1980` or
  `EOPData_IAU2000A` depending on the selection of `data_type`. The returned
  value should be passed to the reference frame conversion functions as
  described in the following.

!!! note

    Notice that, although we can fetch IAU2000A data, this IAU2000A theory is
    not implemented yet.

```jldoctest ECEF_ECI
julia> eop_IAU1980 = get_iers_eop();
[ Info: Downloading file 'EOP_IAU1980.TXT' from 'https://datacenter.iers.org/data/latestVersion/223_EOP_C04_14.62-NOW.IAU1980223.txt'.
```

## ECEF to ECEF

One ECEF frame can be converted to another one by the following function:

```julia
function rECEFtoECEF([T,] ECEFo, ECEFf, JD_UTC::Number, eop_data)
```

where it will be computed the rotation from the ECEF reference frame `ECEFo` to
the ECEF reference frame `ECEFf` at the Julian Day [UTC] `JD_UTC`. The rotation
description that will be used is given by `T`, which can be `DCM` or
`Quaternion`. If `T` is omitted, then it defaults to `DCM`. The EOP data
`eop_data` in this case is always necessary. Hence, the user must initialize it
as described in the section [EOP Data](@ref).

```jldoctest ECEF_ECI
julia> rECEFtoECEF(PEF(), ITRF(), DatetoJD(1986,6,19,21,35,0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
  1.0          0.0         -4.3531e-7
 -6.30011e-13  1.0         -1.44727e-6
  4.3531e-7    1.44727e-6   1.0

julia> rECEFtoECEF(Quaternion, PEF(), ITRF(), DatetoJD(1986,6,19,21,35,0), eop_IAU1980)
Quaternion{Float64}:
  + 0.9999999999997147 - 7.236343481310813e-7.i + 2.1765518308012794e-7.j + 0.0.k
```

## ECI to ECI

One ECI frame can be converted to another ECI frame by one of the following
functions:

```julia
function rECEFtoECI([T,] ECIo, ECIf, JD_UTC::Number [, eop_data])
function rECEFtoECI([T,] ECIo, JD_UTCo::Number, ECIf, JD_UTCf::Number [, eop_data])
```

where it will be computed compute the rotation from the ECI reference frame
`ECIo` to another ECI reference frame `ECIf`. If the origin and destination
frame contain only one *of date* frame, then the first signature is used and
`JD_UTC` is the epoch of this frame. On the other hand, if the origin and
destination frame contain two *of date* frame[^1], e.g. `TOD => MOD`, then the
second signature must be used in which `JD_UTCo` is the epoch of the origin
frame and `JD_UTCf` is the epoch of the destination frame. The rotation
description that will be used is given by `T`, which can be `DCM` or
`Quaternion`. If `T` is omitted, then it defaults to `DCM`. The EOP data
`eop_data`, as described in section [EOP Data](@ref), is required in some
conversions, as described in the following table.

[^1]: TEME is an *of date* frame.

|   Model    |   ECIo  |   ECIf  |    EOP Data   | Function Signature |
|:-----------|:--------|:--------|:--------------|:-------------------|
| IAU-76/FK5 | `GCRF`  | `J2000` | EOP IAU1980   | First              |
| IAU-76/FK5 | `GCRF`  | `MOD`   | EOP IAU1980   | First              |
| IAU-76/FK5 | `GCRF`  | `TOD`   | EOP IAU1980   | First              |
| IAU-76/FK5 | `GCRF`  | `TEME`  | EOP IAU1980   | First              |
| IAU-76/FK5 | `J2000` | `GCRF`  | EOP IAU1980   | First              |
| IAU-76/FK5 | `J2000` | `MOD`   | EOP IAU1980   | First              |
| IAU-76/FK5 | `J2000` | `TOD`   | EOP IAU1980   | First              |
| IAU-76/FK5 | `J2000` | `TEME`  | Not required  | First              |
| IAU-76/FK5 | `MOD`   | `GCRF`  | EOP IAU1980   | First              |
| IAU-76/FK5 | `MOD`   | `J2000` | EOP IAU1980   | First              |
| IAU-76/FK5 | `MOD`   | `TOD`   | EOP IAU1980   | Second             |
| IAU-76/FK5 | `MOD`   | `TEME`  | EOP IAU1980   | Second             |
| IAU-76/FK5 | `TOD`   | `GCRF`  | EOP IAU1980   | First              |
| IAU-76/FK5 | `TOD`   | `J2000` | EOP IAU1980   | First              |
| IAU-76/FK5 | `TOD`   | `MOD`   | EOP IAU1980   | Second             |
| IAU-76/FK5 | `TOD`   | `TEME`  | EOP IAU1980   | Second             |
| IAU-76/FK5 | `TEME`  | `GCRF`  | EOP IAU1980   | First              |
| IAU-76/FK5 | `TEME`  | `J2000` | Not requrired | First              |
| IAU-76/FK5 | `TEME`  | `MOD`   | EOP IAU1980   | Second             |
| IAU-76/FK5 | `TEME`  | `TOD`   | EOP IAU1980   | Second             |

!!! note

    In this function, MOD and TOD frames are always defined with IERS EOP
    corrections. Hence, if one wants to obtain the MOD and TOD frames according
    to the original IAU-76/FK5 theory, it is necessary to use the low-level
    functions in file `./src/transformations/fk5/fk5.jl`.

```jldoctest ECEF_ECI
julia> rECItoECI(DCM, GCRF(), J2000(), DatetoJD(1986, 6, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
  1.0          -2.45469e-12   4.56602e-10
  2.45466e-12   1.0          -1.84455e-9
 -4.56602e-10   1.84455e-9    1.0

julia> rECItoECI(Quaternion, TEME(), GCRF(), DatetoJD(1986, 6, 19, 21, 35, 0), eop_IAU1980)
Quaternion{Float64}:
  + 0.9999986335698654 + 1.8300414020900853e-5.i + 0.0006653038276169474.j - 0.0015132396749411375.k

julia> rECItoECI(TOD(), DatetoJD(1986,6,19,21,35,0), TOD(), DatetoJD(1987,5,19,3,0,0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 1.0          -0.000224087  -9.73784e-5
 0.000224086   1.0          -5.79859e-6
 9.73797e-5    5.77677e-6    1.0

julia> rECItoECI(Quaternion, TOD(), JD_J2000, MOD(), JD_J2000, eop_IAU1980)
Quaternion{Float64}:
  + 0.9999999993282687 - 1.400220690336851e-5.i + 1.3473593746216003e-5.j - 3.107834312843103e-5.k

julia> rECItoECI(J2000(), TEME(), DatetoJD(1986,6,19,21,35,0))
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
  0.999995    0.0030265    0.00133055
 -0.00302645  0.999995    -3.86125e-5
 -0.00133066  3.45854e-5   0.999999
```

## ECEF to ECI

One ECEF frame can be convert to one ECI frame using the following function:

```julia
function rECEFtoECI([T,] ECEF, ECI, JD_UTC::Number [, eop_data])
```

where it will be compute the rotation from the ECEF frame `ECEF` to the ECI
frame `ECI` at the Julian Day [UTC] `JD_UTC`. The rotation description that will
be used is given by `T`, which can be `DCM` or `Quaternion`. If it is omitted,
then it defaults to `DCM`. The EOP data `eop_data`, as described in section [EOP
Data](@ref), is required in some conversions, as described in the following
table.

|   Model    |  ECEF  |   ECI   |    EOP Data    |
|:-----------|:-------|:--------|:---------------|
| IAU-76/FK5 | `ITRF` | `GCRF`  | EOP IAU1980    |
| IAU-76/FK5 | `ITRF` | `J2000` | EOP IAU1980    |
| IAU-76/FK5 | `ITRF` | `MOD`   | EOP IAU1980    |
| IAU-76/FK5 | `ITRF` | `TOD`   | EOP IAU1980    |
| IAU-76/FK5 | `ITRF` | `TEME`  | EOP IAU1980    |
| IAU-76/FK5 | `PEF`  | `GCRF`  | EOP IAU1980    |
| IAU-76/FK5 | `PEF`  | `J2000` | Not required\* |
| IAU-76/FK5 | `PEF`  | `MOD`   | EOP IAU1980    |
| IAU-76/FK5 | `PEF`  | `TOD`   | EOP IAU1980    |
| IAU-76/FK5 | `PEF`  | `TEME`  | Not required\* |

`*`: In this case, the Julian Time UTC will be assumed equal to Julian Time UT1
to compute the Greenwich Mean Sidereal Time. This is an approximation, but
should be sufficiently accurate for some applications. Notice that, if EOP Data
is provided, the Julian Day UT1 will be accurately computed.

!!! note

    In this function, MOD and TOD frames are always defined with IERS EOP
    corrections. Hence, if one wants to obtain the MOD and TOD frames according
    to the original IAU-76/FK5 theory, it is necessary to use the low-level
    functions in file `./src/transformations/fk5/fk5.jl`.

```jldoctest ECEF_ECI
julia> rECEFtoECI(DCM, ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267      0.78518     -0.00132979
 -0.78518      -0.619267     3.33492e-5
 -0.000797313   0.00106478   0.999999

julia> rECEFtoECI(ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267      0.78518     -0.00132979
 -0.78518      -0.619267     3.33492e-5
 -0.000797313   0.00106478   0.999999

julia> rECEFtoECI(PEF(), J2000(), DatetoJD(1986, 06, 19, 21, 35, 0))
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619271      0.785176    -0.00133066
 -0.785177     -0.619272     3.45854e-5
 -0.000796885   0.00106622   0.999999

julia> rECEFtoECI(PEF(), J2000(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267      0.78518     -0.00133066
 -0.78518      -0.619267     3.45854e-5
 -0.000796879   0.00106623   0.999999

julia> rECEFtoECI(Quaternion, ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
Quaternion{Float64}:
  + 0.4363098936462618 - 0.0005909969666939257.i + 0.00030510511316206974.j + 0.8997962182293519.k
```

# ECI to ECEF

One ECI frame can be converted to one ECEF frame using the following function:

```julia
function rECItoECEF([T,] ECI, ECEF, JD_UTC::Number [, eop_data])
```

which has the same characteristics of the function `rECEFtoECI` described in
Section [ECEF to ECI](@ref), but with the inputs `ECI`  and `ECEF` swapped.

!!! note

    This function actually calls `rECEFtoECI` first and then uses
    `inv_rotation`. Hence, it has a slightly overhead on top of `rECEFtoECI`,
    which should be negligible for both rotation representations that are
    supported.

```jldoctest ECEF_ECI
julia> rECItoECEF(DCM, GCRF(), ITRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267    -0.78518     -0.000797313
  0.78518     -0.619267     0.00106478
 -0.00132979   3.33492e-5   0.999999

julia> rECItoECEF(GCRF(), ITRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267    -0.78518     -0.000797313
  0.78518     -0.619267     0.00106478
 -0.00132979   3.33492e-5   0.999999

julia> rECItoECEF(J2000(), PEF(), DatetoJD(1986, 06, 19, 21, 35, 0))
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619271    -0.785177    -0.000796885
  0.785176    -0.619272     0.00106622
 -0.00133066   3.45854e-5   0.999999

julia> rECItoECEF(J2000(), PEF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267    -0.78518     -0.000796879
  0.78518     -0.619267     0.00106623
 -0.00133066   3.45854e-5   0.999999

julia> rECItoECEF(Quaternion, GCRF(), ITRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
Quaternion{Float64}:
  + 0.4363098936462618 + 0.0005909969666939257.i - 0.00030510511316206974.j - 0.8997962182293519.k
```
