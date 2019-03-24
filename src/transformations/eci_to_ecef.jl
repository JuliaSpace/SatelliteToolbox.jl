#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Rotations from an Earth-Fixed Inertial (ECI) reference frame to an
#   Earth-Fixed, Earth-Centered (ECEF) reference frame.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export rECItoECEF

"""
    function rECItoECEF([T,] ECI, ECEF, JD_UTC::Number [, eop_data])

Compute the rotation from an Earth-Centered Inertial (`ECI`) reference frame to
an Earth-Centered, Earth-Fixed (`ECEF`) reference frame at the Julian Day [UTC]
`JD_UTC`. The rotation description that will be used is given by `T`, which can
be `DCM` or `Quaternion`. The ECI frame is selected by the input `ECI` and the
`ECEF` frame is selected by the input `ECEF`. The possible values are listed
below. The model used to compute the rotation is specified by the selection of
the origin and destination frames. Currently, there are two models supported:
IAU-76/FK5 and IAU-2006 with 2010 conventions (CIO approach only).

# Rotation description

The rotations that aligns the ECI with ECEF can be described by Direction Cosine
Matrices or Quaternions. This is selected by the parameter `T`. The possible
values are:

* `DCM`: The rotation will be described by a Direction Cosine Matrix.
* `Quaternion`: The rotation will be described by a Quaternion.

If no value is specified, then it falls back to `DCM`.

# Conversion model

The model that will be used to compute the rotation is automatically inferred
given the selection of the origin and destination frames. **Notice that mixing
IAU-76/FK5 and IAU-2006/2010 frames is not supported yet.**

# ECI Frame

The ECI frame is selected by the parameter `ECI`. The possible values are:

* `TEME()`: ECI will be selected as the True Equator Mean Equinox (TEME)
            reference frame.
* `TOD()`: ECI will be selected as the True of Date (TOD).
* `MOD()`: ECI will be selected as the Mean of Date (MOD).
* `J2000()`: ECI will be selected as the J2000 reference frame.
* `GCRF()`: ECI will be selected as the Geocentric Celestial Reference Frame
            (GCRF).
* `CIRS()`: ECEF will be selected as the Celestial Intermediate Reference System
            (CIRS).

# ECEF Frame

The ECEF frame is selected by the parameter `ECEF`. The possible values are:

* `ITRF()`: ECEF will be selected as the International Terrestrial Reference
            Frame (ITRF).
* `PEF()`: ECEF will be selected as the Pseudo-Earth Fixed (PEF) reference
           frame.
* `TIRS()`: ECEF will be selected as the Terrestrial Intermediate Reference
            System (TIRS).

# EOP Data

The conversion between the frames depends on EOP Data (see `get_iers_eop` and
`read_iers_eop`). If IAU-76/FK5 model is used, then the type of `eop_data` must
be `EOPData_IAU1980`. Otherwise, if IAU-2006/2010 model is used, then the type
of `eop_data` must be `EOPData_IAU2000A`. The following table shows the
requirements for EOP data given the selected frames.

|   Model       |   ECI   |  ECEF  |    EOP Data     |
|:--------------|:--------|:-------|:---------------=|
| IAU-76/FK5    | `GCRF`  | `ITRF` | EOP IAU1980     |
| IAU-76/FK5    | `J2000` | `ITRF` | EOP IAU1980     |
| IAU-76/FK5    | `MOD`   | `ITRF` | EOP IAU1980     |
| IAU-76/FK5    | `TOD`   | `ITRF` | EOP IAU1980     |
| IAU-76/FK5    | `TEME`  | `ITRF` | EOP IAU1980     |
| IAU-76/FK5    | `GCRF`  | `PEF`  | EOP IAU1980     |
| IAU-76/FK5    | `J2000` | `PEF`  | Not required¹   |
| IAU-76/FK5    | `MOD`   | `PEF`  | Not required¹   |
| IAU-76/FK5    | `TOD`   | `PEF`  | Not required¹   |
| IAU-76/FK5    | `TEME`  | `PEF`  | Not required¹   |
| IAU-2006/2010 | `CIRS`  | `ITRF` | EOP IAU2000A    |
| IAU-2006/2010 | `GCRF`  | `ITRF` | EOP IAU2000A    |
| IAU-2006/2010 | `CIRS`  | `TIRS` | Not required¹   |
| IAU-2006/2010 | `GCRF`  | `TIRS` | Not required¹ ² |

`¹`: In this case, the Julian Time UTC will be assumed equal to Julian Time UT1
to compute the Greenwich Mean Sidereal Time. This is an approximation, but
should be sufficiently accurate for some applications. Notice that, if EOP Data
is provided, the Julian Day UT1 will be accurately computed.

`²`: In this case, the terms that account for the free-core nutation and time
dependent effects of the Celestial Intermediate Pole (CIP) position with respect
to the GCRF will not be available, reducing the precision.
The conversion between the frames depends on EOP Data (see `get_iers_eop` and
`read_iers_eop`). If IAU-76/FK5 model is used, then the type of `eop_data` must
be `EOPData_IAU1980`. The following table shows the requirements for EOP data
given the selected frames.

## MOD and TOD

In this function, if EOP corrections are not provided, then MOD and TOD frames
will be computed considering the original IAU-76/FK5 theory. Otherwise, the
corrected frame will be used.

# Returns

The rotation description represented by `T` that rotates the ECI reference frame
into alignment with the ECEF reference frame.

# Examples

```julia-repl
julia> eop_IAU1980 = get_iers_eop(:IAU1980);

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

julia> eop_IAU2000A = get_iers_eop(:IAU2000A);

julia> rECItoECEF(GCRF(), ITRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU2000A)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267    -0.78518     -0.000797312
  0.78518     -0.619267     0.00106478
 -0.00132979   3.33502e-5   0.999999

julia> rECItoECEF(GCRF(), TIRS(), DatetoJD(1986, 06, 19, 21, 35, 0))
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619271    -0.785177    -0.000796885
  0.785176    -0.619272     0.00106623
 -0.00133066   3.45884e-5   0.999999

julia> rECItoECEF(Quaternion, GCRF(), ITRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU2000A)
Quaternion{Float64}:
  + 0.4363098936309669 + 0.000590996988144556.i - 0.0003051056555230158.j - 0.8997962182365703.k
```
"""
@inline rECItoECEF(T_ECI::T_ECIs,
                   T_ECEF::T_ECEFs,
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980) =
    rECItoECEF(DCM, T_ECI, T_ECEF, JD_UTC, eop_data)

@inline rECItoECEF(T::T_ROT,
                   T_ECI::T_ECIs,
                   T_ECEF::T_ECEFs,
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980) =
    inv_rotation(rECEFtoECI(T, T_ECEF, T_ECI, JD_UTC, eop_data))

@inline rECItoECEF(T_ECI::T_ECIs_IAU_2006,
                   T_ECEF::T_ECEFs_IAU_2006,
                   JD_UTC::Number,
                   eop_data::EOPData_IAU2000A) =
    rECItoECEF(DCM, T_ECI, T_ECEF, JD_UTC, eop_data)

@inline rECItoECEF(T::T_ROT,
                   T_ECI::T_ECIs_IAU_2006,
                   T_ECEF::T_ECEFs_IAU_2006,
                   JD_UTC::Number,
                   eop_data::EOPData_IAU2000A) =
    inv_rotation(rECEFtoECI(T, T_ECEF, T_ECI, JD_UTC, eop_data))

# Specializations for those cases that EOP Data is not needed.
@inline rECItoECEF(T_ECI::Union{Type{Val{:J2000}},Type{Val{:TOD}},
                                Type{Val{:MOD}},Type{Val{:TEME}}},
                   T_ECEF::Type{Val{:PEF}}, JD_UTC::Number) =
    rECItoECEF(DCM, T_ECI, T_ECEF, JD_UTC)

@inline rECItoECEF(T::T_ROT,
                   T_ECI::Union{Type{Val{:J2000}},Type{Val{:TOD}},
                                Type{Val{:MOD}},Type{Val{:TEME}}},
                   T_ECEF::Type{Val{:PEF}}, JD_UTC::Number) =
    inv_rotation(rECEFtoECI(T, T_ECEF, T_ECI, JD_UTC))

@inline rECItoECEF(T_ECI::Union{Type{Val{:CIRS}},Type{Val{:GCRF}}},
                   T_ECEF::Type{Val{:TIRS}},
                   JD_UTC::Number) =
    rECItoECEF(DCM, T_ECI, T_ECEF, JD_UTC)

@inline rECItoECEF(T::T_ROT,
                   T_ECI::Union{Type{Val{:CIRS}},Type{Val{:GCRF}}},
                   T_ECEF::Type{Val{:TIRS}},
                   JD_UTC::Number) =
    inv_rotation(rECEFtoECI(T, T_ECEF, T_ECI, JD_UTC))
