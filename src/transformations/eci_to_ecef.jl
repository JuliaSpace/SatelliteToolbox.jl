# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Rotations from an Earth-Fixed Inertial (ECI) reference frame to an
#   Earth-Fixed, Earth-Centered (ECEF) reference frame.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export r_eci_to_ecef

"""
    r_eci_to_ecef([T,] ECI, ECEF, JD_UTC::Number [, eop_data])

Compute the rotation from an Earth-Centered Inertial (`ECI`) reference frame to
an Earth-Centered, Earth-Fixed (`ECEF`) reference frame at the Julian Day [UTC]
`JD_UTC`. The rotation description that will be used is given by `T`, which can
be `DCM` or `Quaternion`. The ECI frame is selected by the input `ECI` and the
`ECEF` frame is selected by the input `ECEF`. The possible values are listed
below. The model used to compute the rotation is specified by the selection of
the origin and destination frames. Currently, there are two models supported:
IAU-76/FK5 and IAU-2006 with 2010 conventions (CIO and equinox approaches).

# Rotation description

The rotations that aligns the ECI with ECEF can be described by Direction Cosine
Matrices or Quaternions. This is selected by the parameter `T`. The possible
values are:

- `DCM`: The rotation will be described by a Direction Cosine Matrix.
- `Quaternion`: The rotation will be described by a Quaternion.

If no value is specified, then it falls back to `DCM`.

# Conversion model

The model that will be used to compute the rotation is automatically inferred
given the selection of the origin and destination frames. **Notice that mixing
IAU-76/FK5 and IAU-2006/2010 frames is not supported.**

# ECI Frame

The ECI frame is selected by the parameter `ECI`. The possible values are:

- `TEME()`: ECI will be selected as the True Equator Mean Equinox (TEME)
    reference frame.
- `TOD()`: ECI will be selected as the True of Date (TOD).
- `MOD()`: ECI will be selected as the Mean of Date (MOD).
- `J2000()`: ECI will be selected as the J2000 reference frame.
- `GCRF()`: ECI will be selected as the Geocentric Celestial Reference Frame
    (GCRF).
- `CIRS()`: ECEF will be selected as the Celestial Intermediate Reference System
    (CIRS).
- `ERS()`: ECI will be selected as the Earth Reference System (ERS).
- `MOD06()`: ECI will be selected as the Mean of Date (MOD) according to the
    definition in IAU-2006/2010 theory.
- `MJ2000()`: ECI will be selected as the J2000 mean equatorial frame (MJ2000).

!!! note
    The frames `MOD()` and `MOD06()` are virtually the same. However, we
    selected different names to make clear which theory are being used since
    mixing transformation between frames from IAU-76/FK5 and IAU-2006/2010 must
    be performed with caution.

# ECEF Frame

The ECEF frame is selected by the parameter `ECEF`. The possible values are:

- `ITRF()`: ECEF will be selected as the International Terrestrial Reference
    Frame (ITRF).
- `PEF()`: ECEF will be selected as the Pseudo-Earth Fixed (PEF) reference
    frame.
- `TIRS()`: ECEF will be selected as the Terrestrial Intermediate Reference
    System (TIRS).

# EOP Data

The conversion between the frames depends on EOP Data (see
[`get_iers_eop`](@ref) and [`read_iers_eop`](@ref)). If IAU-76/FK5 model is
used, then the type of `eop_data` must be [`EOPData_IAU1980`](@ref). Otherwise,
if IAU-2006/2010 model is used, then the type of `eop_data` must be
[`EOPData_IAU2000A`](@ref). The following table shows the requirements for EOP
data given the selected frames.

|   Model                     |   ECI    |  ECEF  |    EOP Data     |
|:----------------------------|:---------|:-------|:----------------|
| IAU-76/FK5                  | `GCRF`   | `ITRF` | EOP IAU1980     |
| IAU-76/FK5                  | `J2000`  | `ITRF` | EOP IAU1980     |
| IAU-76/FK5                  | `MOD`    | `ITRF` | EOP IAU1980     |
| IAU-76/FK5                  | `TOD`    | `ITRF` | EOP IAU1980     |
| IAU-76/FK5                  | `TEME`   | `ITRF` | EOP IAU1980     |
| IAU-76/FK5                  | `GCRF`   | `PEF`  | EOP IAU1980     |
| IAU-76/FK5                  | `J2000`  | `PEF`  | Not required¹   |
| IAU-76/FK5                  | `MOD`    | `PEF`  | Not required¹   |
| IAU-76/FK5                  | `TOD`    | `PEF`  | Not required¹   |
| IAU-76/FK5                  | `TEME`   | `PEF`  | Not required¹   |
| IAU-2006/2010 CIO-based     | `CIRS`   | `ITRF` | EOP IAU2000A    |
| IAU-2006/2010 CIO-based     | `GCRF`   | `ITRF` | EOP IAU2000A    |
| IAU-2006/2010 CIO-based     | `CIRS`   | `TIRS` | Not required¹   |
| IAU-2006/2010 CIO-based     | `GCRF`   | `TIRS` | Not required¹ ² |
| IAU-2006/2010 Equinox-based | `ERS`    | `TIRS` | EOP IAU2000A    |
| IAU-2006/2010 Equinox-based | `MOD06`  | `ITRF` | EOP IAU2000A    |
| IAU-2006/2010 Equinox-based | `MJ2000` | `ITRF` | EOP IAU2000A    |
| IAU-2006/2010 Equinox-based | `ERS`    | `TIRS` | Not required¹ ³ |
| IAU-2006/2010 Equinox-based | `MOD06`  | `TIRS` | Not required¹ ³ |
| IAU-2006/2010 Equinox-based | `MJ2000` | `TIRS` | Not required¹ ³ |

`¹`: In this case, UTC will be assumed equal to UT1 to compute the Greenwich
Mean Sidereal Time. This is an approximation, but should be sufficiently
accurate for some applications. Notice that, if EOP Data is provided, UT1 will
be accurately computed.

`²`: In this case, the terms that account for the free-core nutation and time
dependent effects of the Celestial Intermediate Pole (CIP) position with respect
to the GCRF will not be available, reducing the precision.

## MOD and TOD

In this function, if EOP corrections are not provided, then MOD and TOD frames
will be computed considering the original IAU-76/FK5 theory. Otherwise, the
corrected frame will be used.

# Returns

The rotation description represented by `T` that rotates the ECI reference frame
into alignment with the ECEF reference frame.

# Examples

```julia-repl
julia> eop_IAU1980 = get_iers_eop(Val(:IAU1980));

julia> r_eci_to_ecef(DCM, GCRF(), ITRF(), date_to_jd(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.619267    -0.78518     -0.000797312
  0.78518     -0.619267     0.00106478
 -0.00132979   3.33509e-5   0.999999

julia> r_eci_to_ecef(GCRF(), ITRF(), date_to_jd(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.619267    -0.78518     -0.000797312
  0.78518     -0.619267     0.00106478
 -0.00132979   3.33509e-5   0.999999

julia> r_eci_to_ecef(J2000(), PEF(), date_to_jd(1986, 06, 19, 21, 35, 0))
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.619271    -0.785177    -0.000796885
  0.785176    -0.619272     0.00106622
 -0.00133066   3.45854e-5   0.999999

julia> r_eci_to_ecef(J2000(), PEF(), date_to_jd(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.619267    -0.78518     -0.000796879
  0.78518     -0.619267     0.00106623
 -0.00133066   3.45854e-5   0.999999

julia> r_eci_to_ecef(Quaternion, GCRF(), ITRF(), date_to_jd(1986, 06, 19, 21, 35, 0), eop_IAU1980)
Quaternion{Float64}:
  + 0.43631 + 0.000590997⋅i - 0.000305106⋅j - 0.000305106⋅k

julia> eop_IAU2000A = get_iers_eop(Val(:IAU2000A));

julia> r_eci_to_ecef(GCRF(), ITRF(), date_to_jd(1986, 06, 19, 21, 35, 0), eop_IAU2000A)
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.619267    -0.78518     -0.000797311
  0.78518     -0.619267     0.00106478
 -0.00132979   3.33516e-5   0.999999

julia> r_eci_to_ecef(GCRF(), TIRS(), date_to_jd(1986, 06, 19, 21, 35, 0))
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.619271    -0.785177    -0.000796885
  0.785176    -0.619272     0.00106623
 -0.00133066   3.45884e-5   0.999999

julia> r_eci_to_ecef(Quaternion, GCRF(), ITRF(), date_to_jd(1986, 06, 19, 21, 35, 0), eop_IAU2000A)
Quaternion{Float64}:
  + 0.43631 + 0.000590997⋅i - 0.000305106⋅j - 0.000305106⋅k
```
"""
@inline function r_eci_to_ecef(T_ECI::T_ECIs,
    T_ECEF::T_ECEFs,
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    return r_eci_to_ecef(DCM, T_ECI, T_ECEF, JD_UTC, eop_data)
end

@inline function r_eci_to_ecef(
    T::T_ROT,
    T_ECI::T_ECIs,
    T_ECEF::T_ECEFs,
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_ecef_to_eci(T, T_ECEF, T_ECI, JD_UTC, eop_data))
end

@inline function r_eci_to_ecef(
    T_ECI::T_ECIs_IAU_2006,
    T_ECEF::T_ECEFs_IAU_2006,
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_ecef(DCM, T_ECI, T_ECEF, JD_UTC, eop_data)
end

@inline function r_eci_to_ecef(
    T::T_ROT,
    T_ECI::T_ECIs_IAU_2006,
    T_ECEF::T_ECEFs_IAU_2006,
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    return inv_rotation(r_ecef_to_eci(T, T_ECEF, T_ECI, JD_UTC, eop_data))
end

# Specializations for those cases that EOP Data is not needed.
@inline function r_eci_to_ecef(
    T_ECI::Union{Val{:J2000}, Val{:TOD}, Val{:MOD}, Val{:TEME}},
    T_ECEF::Val{:PEF},
    JD_UTC::Number
)
    return r_eci_to_ecef(DCM, T_ECI, T_ECEF, JD_UTC)
end

@inline function r_eci_to_ecef(
    T::T_ROT,
    T_ECI::Union{Val{:J2000}, Val{:TOD}, Val{:MOD}, Val{:TEME}},
    T_ECEF::Val{:PEF},
    JD_UTC::Number
)
    return inv_rotation(r_ecef_to_eci(T, T_ECEF, T_ECI, JD_UTC))
end

@inline function r_eci_to_ecef(
    T_ECI::T_ECIs_IAU_2006,
    T_ECEF::Val{:TIRS},
    JD_UTC::Number
)
    return r_eci_to_ecef(DCM, T_ECI, T_ECEF, JD_UTC)
end

@inline function r_eci_to_ecef(
    T::T_ROT,
    T_ECI::T_ECIs_IAU_2006,
    T_ECEF::Val{:TIRS},
    JD_UTC::Number
)
    return inv_rotation(r_ecef_to_eci(T, T_ECEF, T_ECI, JD_UTC))
end
