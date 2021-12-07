# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Rotations from an Earth-Fixed, Earth-Centered (ECEF) reference frame to an
#   Earth-Fixed Inertial (ECI) reference frame.
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

export r_ecef_to_eci

"""
    r_ecef_to_eci([T,] ECEF, ECI, jd_utc::Number [, eop_data])

Compute the rotation from an Earth-Centered, Earth-Fixed (`ECEF`) reference
frame to an Earth-Centered Inertial (`ECI`) reference frame at the Julian Day
[UTC] `jd_utc`. The rotation description that will be used is given by `T`,
which can be `DCM` or `Quaternion`. The ECEF frame is selected by the input
`ECEF` and the `ECI` frame is selected by the input `ECI`. The possible values
are listed below. The model used to compute the rotation is specified by the
selection of the origin and destination frames. Currently, there are two models
supported: IAU-76/FK5 and IAU-2006 with 2010 conventions (CIO and equinox
approaches).

# Rotation description

The rotations that aligns the ECEF with ECI can be described by Direction Cosine
Matrices or Quaternions. This is selected by the parameter `T`. The possible
values are:

- `DCM`: The rotation will be described by a Direction Cosine Matrix.
- `Quaternion`: The rotation will be described by a Quaternion.

If no value is specified, then it falls back to `DCM`.

# Conversion model

The model that will be used to compute the rotation is automatically inferred
given the selection of the origin and destination frames. **Notice that mixing
IAU-76/FK5 and IAU-2006/2010 frames is not supported.**

# ECEF Frame

The ECEF frame is selected by the parameter `ECEF`. The possible values are:

- `ITRF()`: ECEF will be selected as the International Terrestrial Reference
    Frame (ITRF).
- `PEF()`: ECEF will be selected as the Pseudo-Earth Fixed (PEF) reference
    frame.
- `TIRS()`: ECEF will be selected as the Terrestrial Intermediate Reference
    System (TIRS).

# ECI Frame

The ECI frame is selected by the parameter `ECI`. The possible values are:

- `TEME()`: ECI will be selected as the True Equator Mean Equinox (TEME)
    reference frame.
- `TOD()`: ECI will be selected as the True of Date (TOD).
- `MOD()`: ECI will be selected as the Mean of Date (MOD).
- `J2000()`: ECI will be selected as the J2000 reference frame.
- `GCRF()`: ECI will be selected as the Geocentric Celestial Reference Frame
    (GCRF).
- `CIRS()`: ECI will be selected as the Celestial Intermediate Reference System
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

# EOP Data

The conversion between the frames depends on EOP Data (see
[`get_iers_eop`](@ref) and [`read_iers_eop`](@ref)). If IAU-76/FK5 model is
used, then the type of `eop_data` must be [`EOPData_IAU1980`](@ref). Otherwise,
if IAU-2006/2010 model is used, then the type of `eop_data` must be
[`EOPData_IAU2000A`](@ref). The following table shows the requirements for EOP
data given the selected frames.

|   Model                     |  ECEF  |   ECI    |    EOP Data     |
|:----------------------------|:-------|:---------|:----------------|
| IAU-76/FK5                  | `ITRF` | `GCRF`   | EOP IAU1980     |
| IAU-76/FK5                  | `ITRF` | `J2000`  | EOP IAU1980     |
| IAU-76/FK5                  | `ITRF` | `MOD`    | EOP IAU1980     |
| IAU-76/FK5                  | `ITRF` | `TOD`    | EOP IAU1980     |
| IAU-76/FK5                  | `ITRF` | `TEME`   | EOP IAU1980     |
| IAU-76/FK5                  | `PEF`  | `GCRF`   | EOP IAU1980     |
| IAU-76/FK5                  | `PEF`  | `J2000`  | Not required¹   |
| IAU-76/FK5                  | `PEF`  | `MOD`    | Not required¹   |
| IAU-76/FK5                  | `PEF`  | `TOD`    | Not required¹   |
| IAU-76/FK5                  | `PEF`  | `TEME`   | Not required¹   |
| IAU-2006/2010 CIO-based     | `ITRF` | `CIRS`   | EOP IAU2000A    |
| IAU-2006/2010 CIO-based     | `ITRF` | `GCRF`   | EOP IAU2000A    |
| IAU-2006/2010 CIO-based     | `TIRS` | `CIRS`   | Not required¹   |
| IAU-2006/2010 CIO-based     | `TIRS` | `GCRF`   | Not required¹ ² |
| IAU-2006/2010 Equinox-based | `ITRF` | `ERS`    | EOP IAU2000A    |
| IAU-2006/2010 Equinox-based | `ITRF` | `MOD06`  | EOP IAU2000A    |
| IAU-2006/2010 Equinox-based | `ITRF` | `MJ2000` | EOP IAU2000A    |
| IAU-2006/2010 Equinox-based | `TIRS` | `ERS`    | Not required¹ ³ |
| IAU-2006/2010 Equinox-based | `TIRS` | `MOD06`  | Not required¹ ³ |
| IAU-2006/2010 Equinox-based | `TIRS` | `MJ2000` | Not required¹ ³ |

`¹`: In this case, UTC will be assumed equal to UT1 to compute the Greenwich
Mean Sidereal Time. This is an approximation, but should be sufficiently
accurate for some applications. Notice that, if EOP Data is provided, UT1 will
be accurately computed.

`²`: In this case, the terms that account for the free core nutation and time
dependent effects of the Celestial Intermediate Pole (CIP) position with respect
to the GCRF will not be available, reducing the precision.

`³`: In this case, the terms that corrects the nutation in obliquity and in
longitude due to the free core nutation will not be available, reducing the
precision.

## MOD and TOD

In this function, if EOP corrections are not provided, then MOD and TOD frames
will be computed considering the original IAU-76/FK5 theory. Otherwise, the
corrected frame will be used.

# Returns

The rotation description represented by `T` that rotates the ECEF reference
frame into alignment with the ECI reference frame.

# Examples

```julia-repl
julia> eop_IAU1980 = get_iers_eop(Val(:IAU1980));

julia> r_ecef_to_eci(DCM, ITRF(), GCRF(), date_to_jd(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.619267      0.78518     -0.00132979
 -0.78518      -0.619267     3.33509e-5
 -0.000797312   0.00106478   0.999999

julia> r_ecef_to_eci(ITRF(), GCRF(), date_to_jd(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.619267      0.78518     -0.00132979
 -0.78518      -0.619267     3.33509e-5
 -0.000797312   0.00106478   0.999999

julia> r_ecef_to_eci(PEF(), J2000(), date_to_jd(1986, 06, 19, 21, 35, 0))
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.619271      0.785176    -0.00133066
 -0.785177     -0.619272     3.45854e-5
 -0.000796885   0.00106622   0.999999

julia> r_ecef_to_eci(PEF(), J2000(), date_to_jd(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.619267      0.78518     -0.00133066
 -0.78518      -0.619267     3.45854e-5
 -0.000796879   0.00106623   0.999999

julia> r_ecef_to_eci(Quaternion, ITRF(), GCRF(), date_to_jd(1986, 06, 19, 21, 35, 0), eop_IAU1980)
Quaternion{Float64}:
  + 0.43631 - 0.000590997⋅i + 0.000305106⋅j + 0.000305106⋅k

julia> eop_IAU2000A = get_iers_eop(Val(:IAU2000A));

julia> r_ecef_to_eci(ITRF(), GCRF(), date_to_jd(1986, 06, 19, 21, 35, 0), eop_IAU2000A)
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.619267      0.78518     -0.00132979
 -0.78518      -0.619267     3.33516e-5
 -0.000797311   0.00106478   0.999999

julia> r_ecef_to_eci(TIRS(), GCRF(), date_to_jd(1986, 06, 19, 21, 35, 0))
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 -0.619271      0.785176    -0.00133066
 -0.785177     -0.619272     3.45884e-5
 -0.000796885   0.00106623   0.999999

julia> r_ecef_to_eci(Quaternion, ITRF(), GCRF(), date_to_jd(1986, 06, 19, 21, 35, 0), eop_IAU2000A)
Quaternion{Float64}:
  + 0.43631 - 0.000590997⋅i + 0.000305106⋅j + 0.000305106⋅k
```
"""
@inline function r_ecef_to_eci(
    T_ECEF::T_ECEFs,
    T_ECI::T_ECIs,
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    r_ecef_to_eci(DCM, T_ECEF, T_ECI, jd_utc, eop_data)
end

@inline function r_ecef_to_eci(
    T_ECEF::T_ECEFs_IAU_2006,
    T_ECI::T_ECIs_IAU_2006,
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    r_ecef_to_eci(DCM, T_ECEF, T_ECI, jd_utc, eop_data)
end

# Specializations for those cases that EOP Data is not needed.
@inline function r_ecef_to_eci(
    T_ECEF::Val{:PEF},
    T_ECI::Union{Val{:J2000}, Val{:MOD}, Val{:TOD}, Val{:TEME}},
    jd_utc::Number
)
    r_ecef_to_eci(DCM, T_ECEF, T_ECI, jd_utc)
end

@inline function r_ecef_to_eci(
    T_ECEF::Val{:TIRS},
    T_ECI::T_ECIs_IAU_2006,
    jd_utc::Number
)
    r_ecef_to_eci(DCM, T_ECEF, T_ECI, jd_utc)
end

################################################################################
#                                  IAU-76/FK5
################################################################################

#                                 ITRF => GCRF
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:GCRF},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    arcsec_to_rad = π / 648000
    milliarcsec_to_rad = arcsec_to_rad / 1000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p      = eop_data.x(jd_utc) * arcsec_to_rad
    y_p      = eop_data.y(jd_utc) * arcsec_to_rad
    δΔϵ_1980 = eop_data.dEps(jd_utc) * milliarcsec_to_rad
    δΔψ_1980 = eop_data.dPsi(jd_utc) * milliarcsec_to_rad

    # Compute the rotation.
    return r_itrf_to_gcrf_fk5(T, jd_ut1, jd_tt, x_p, y_p, δΔϵ_1980, δΔψ_1980)
end

#                                ITRF => J2000
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:J2000},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    arcsec_to_rad = π / 648000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p = eop_data.x(jd_utc) * arcsec_to_rad
    y_p = eop_data.y(jd_utc) * arcsec_to_rad

    # Compute the rotation.
    return r_itrf_to_gcrf_fk5(T, jd_ut1, jd_tt, x_p, y_p, 0, 0)
end

#                                 ITRF => MOD
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:MOD},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    arcsec_to_rad = π / 648000
    milliarcsec_to_rad = arcsec_to_rad / 1000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p      = eop_data.x(jd_utc) * arcsec_to_rad
    y_p      = eop_data.y(jd_utc) * arcsec_to_rad
    δΔϵ_1980 = eop_data.dEps(jd_utc) * milliarcsec_to_rad
    δΔψ_1980 = eop_data.dPsi(jd_utc) * milliarcsec_to_rad

    # Compute the rotation.
    r_pef_itrf = r_itrf_to_pef_fk5(T, x_p, y_p)
    r_mod_pef  = r_pef_to_mod_fk5(T, jd_ut1, jd_tt, δΔϵ_1980, δΔψ_1980)

    return compose_rotation(r_pef_itrf, r_mod_pef)
end

#                                 ITRF => TOD
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:TOD},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    arcsec_to_rad = π / 648000
    milliarcsec_to_rad = arcsec_to_rad / 1000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p      = eop_data.x(jd_utc) * arcsec_to_rad
    y_p      = eop_data.y(jd_utc) * arcsec_to_rad
    δΔψ_1980 = eop_data.dPsi(jd_utc) * milliarcsec_to_rad

    # Compute the rotation.
    r_pef_itrf = r_itrf_to_pef_fk5(T, x_p, y_p)
    r_tod_pef  = r_pef_to_tod_fk5(T, jd_ut1, jd_tt, δΔψ_1980)

    return compose_rotation(r_pef_itrf, r_tod_pef)
end

#                                 ITRF => TEME
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:TEME},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    arcsec_to_rad = π / 648000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p = eop_data.x(jd_utc) * arcsec_to_rad
    y_p = eop_data.y(jd_utc) * arcsec_to_rad

    # Compute the rotation.
    r_pef_itrf = r_itrf_to_pef_fk5(T, x_p, y_p)
    r_teme_pef = r_pef_to_teme(T, jd_ut1)

    return compose_rotation(r_pef_itrf, r_teme_pef)
end

#                                 PEF => GCRF
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:PEF},
    ::Val{:GCRF},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    arcsec_to_rad = π / 648000
    milliarcsec_to_rad = arcsec_to_rad / 1000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p      = eop_data.x(jd_utc) * arcsec_to_rad
    y_p      = eop_data.y(jd_utc) * arcsec_to_rad
    δΔϵ_1980 = eop_data.dEps(jd_utc) * milliarcsec_to_rad
    δΔψ_1980 = eop_data.dPsi(jd_utc) * milliarcsec_to_rad

    # Compute the rotation.
    r_mod_pef  = r_pef_to_mod_fk5(T, jd_ut1, jd_tt, δΔϵ_1980, δΔψ_1980)
    r_gcrf_mod = r_mod_to_gcrf_fk5(T, jd_tt)

    return compose_rotation(r_mod_pef, r_gcrf_mod)
end

#                                PEF => J2000
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:PEF},
    ::Val{:J2000},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Compute the rotation.
    r_mod_pef  = r_pef_to_mod_fk5(T, jd_ut1, jd_tt, 0, 0)
    r_gcrf_mod = r_mod_to_gcrf_fk5(T, jd_tt)

    return compose_rotation(r_mod_pef, r_gcrf_mod)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:PEF}, ::Val{:J2000}, jd_utc::Number)
    # Since we do not have EOP Data, assume that jd_utc is equal to jd_ut1.
    jd_ut1 = jd_utc
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Compute the rotation.
    r_mod_pef  = r_pef_to_mod_fk5(T, jd_ut1, jd_tt, 0, 0)
    r_gcrf_mod = r_mod_to_gcrf_fk5(T, jd_tt)

    return compose_rotation(r_mod_pef, r_gcrf_mod)
end

#                                 PEF => MOD
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:PEF},
    ::Val{:MOD},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    arcsec_to_rad = π / 648000
    milliarcsec_to_rad = arcsec_to_rad / 1000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p      = eop_data.x(jd_utc) * arcsec_to_rad
    y_p      = eop_data.y(jd_utc) * arcsec_to_rad
    δΔϵ_1980 = eop_data.dEps(jd_utc) * milliarcsec_to_rad
    δΔψ_1980 = eop_data.dPsi(jd_utc) * milliarcsec_to_rad

    # Compute the rotation.
    return r_pef_to_mod_fk5(T, jd_ut1, jd_tt, δΔϵ_1980, δΔψ_1980)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:PEF}, ::Val{:MOD}, jd_utc::Number)
    # Since we do not have EOP Data, assume that jd_utc is equal to jd_ut1.
    jd_ut1 = jd_utc
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Compute the rotation.
    return r_pef_to_mod_fk5(T, jd_ut1, jd_tt, 0, 0)
end

#                                 PEF => TOD
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:PEF},
    ::Val{:TOD},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    arcsec_to_rad = π / 648000
    milliarcsec_to_rad = arcsec_to_rad / 1000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p      = eop_data.x(jd_utc) * arcsec_to_rad
    y_p      = eop_data.y(jd_utc) * arcsec_to_rad
    δΔψ_1980 = eop_data.dPsi(jd_utc) * milliarcsec_to_rad

    # Compute the rotation.
    return r_pef_to_tod_fk5(T, jd_ut1, jd_tt, δΔψ_1980)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:PEF}, ::Val{:TOD}, jd_utc::Number)
    # Since we do not have EOP Data, assume that jd_utc is equal to jd_ut1.
    jd_ut1 = jd_utc
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Compute the rotation.
    return r_pef_to_tod_fk5(T, jd_ut1, jd_tt, 0)
end

#                                 PEF => TEME
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:PEF},
    ::Val{:TEME},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    # Get the time in UT1.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)

    # Compute the rotation.
    return r_pef_to_teme(T, jd_ut1)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:PEF}, ::Val{:TEME}, jd_utc::Number)
    # Since we do not have EOP Data, assume that jd_utc is equal to jd_ut1.
    jd_ut1 = jd_utc

    # Compute the rotation.
    return r_pef_to_teme(T, jd_ut1)
end

################################################################################
#                           IAU-2006/2010 CIO-based
################################################################################

#                                 ITRF => CIRS
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:CIRS},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec_to_rad = π / 648000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p = eop_data.x(jd_utc) * arcsec_to_rad
    y_p = eop_data.y(jd_utc) * arcsec_to_rad

    # Compute the rotation.
    r_tirs_itrf = r_itrf_to_tirs_iau2006(T, jd_tt, x_p, y_p)
    r_cirs_tirs = r_tirs_to_cirs_iau2006(T, jd_ut1)

    return compose_rotation(r_tirs_itrf, r_cirs_tirs)
end

#                                 ITRF => GCRF
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:GCRF},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec_to_rad = π / 648000
    milliarcsec_to_rad = arcsec_to_rad / 1000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p = eop_data.x(jd_utc) * arcsec_to_rad
    y_p = eop_data.y(jd_utc) * arcsec_to_rad
    dx  = eop_data.dX(jd_utc) * milliarcsec_to_rad
    dy  = eop_data.dY(jd_utc) * milliarcsec_to_rad

    # Compute the rotation.
    r_tirs_itrf = r_itrf_to_tirs_iau2006(T, jd_tt, x_p, y_p)
    r_cirs_tirs = r_tirs_to_cirs_iau2006(T, jd_ut1)
    r_gcrf_cirs = r_cirs_to_gcrf_iau2006(T, jd_tt, dx, dy)

    return compose_rotation(r_tirs_itrf, r_cirs_tirs, r_gcrf_cirs)
end

#                                 TIRS => CIRS
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:TIRS},
    ::Val{:CIRS},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)

    # Compute the rotation.
    return r_tirs_to_cirs_iau2006(T, jd_ut1)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:TIRS}, ::Val{:CIRS}, jd_utc::Number)
    # Since we do not have EOP Data, assume that jd_utc is equal to jd_ut1.
    jd_ut1 = jd_utc

    # Compute the rotation.
    return r_tirs_to_cirs_iau2006(T, jd_ut1)
end

#                                 TIRS => GCRF
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:TIRS},
    ::Val{:GCRF},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    milliarcsec_to_rad = π / 648000000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    dx = eop_data.dX(jd_utc) * milliarcsec_to_rad
    dy = eop_data.dY(jd_utc) * milliarcsec_to_rad

    # Compute the rotation.
    r_cirs_tirs = r_tirs_to_cirs_iau2006(T, jd_ut1)
    r_gcrf_cirs = r_cirs_to_gcrf_iau2006(T, jd_tt, dx, dy)

    return compose_rotation(r_cirs_tirs, r_gcrf_cirs)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:TIRS}, ::Val{:GCRF}, jd_utc::Number)
    # Since we do not have EOP Data, assume that jd_utc is equal to jd_ut1.
    jd_ut1 = jd_utc

    # Get the time in TT.
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Compute the rotation.
    r_cirs_tirs = r_tirs_to_cirs_iau2006(T, jd_ut1)
    r_gcrf_cirs = r_cirs_to_gcrf_iau2006(T, jd_tt)

    return compose_rotation(r_cirs_tirs, r_gcrf_cirs)
end

################################################################################
#                         IAU-2006/2010 equinox-based
################################################################################

#                                 ITRF => ERS
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:ERS},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec_to_rad = π / 648000
    milliarcsec_to_rad = arcsec_to_rad / 1000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p = eop_data.x(jd_utc) * arcsec_to_rad
    y_p = eop_data.y(jd_utc) * arcsec_to_rad

    # Obtain the correction of the nutation in longitude.
    ~, δΔΨ_2000 = deps_dpsi(eop_data, jd_utc)
    δΔΨ_2000 *= milliarcsec_to_rad

    # Compute the rotation.
    r_tirs_itrf = r_itrf_to_tirs_iau2006(T, jd_tt, x_p, y_p)
    r_ERS_TIRS = r_tirs_to_ers_iau2006(T, jd_ut1, jd_tt, δΔΨ_2000)

    return compose_rotation(r_tirs_itrf, r_ERS_TIRS)
end

#                                 ITRF => MOD
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:MOD06},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec_to_rad = π / 648000
    milliarcsec_to_rad = arcsec_to_rad / 1000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p = eop_data.x(jd_utc) * arcsec_to_rad
    y_p = eop_data.y(jd_utc) * arcsec_to_rad

    # Obtain the correction of the nutation in obliquity and longitude.
    δΔϵ_2000, δΔΨ_2000 = deps_dpsi(eop_data, jd_utc)
    δΔϵ_2000 *= milliarcsec_to_rad
    δΔΨ_2000 *= milliarcsec_to_rad

    # Compute the rotation.
    r_tirs_itrf = r_itrf_to_tirs_iau2006(T, jd_tt, x_p, y_p)
    r_mod_tirs = r_tirs_to_mod_iau2006(T, jd_ut1, jd_tt, δΔϵ_2000, δΔΨ_2000)

    return compose_rotation(r_tirs_itrf, r_mod_tirs)
end

#                                ITRF => MJ2000
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:MJ2000},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec_to_rad = π / 648000
    milliarcsec_to_rad = arcsec_to_rad / 1000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    x_p = eop_data.x(jd_utc) * arcsec_to_rad
    y_p = eop_data.y(jd_utc) * arcsec_to_rad

    # Obtain the correction of the nutation in obliquity and longitude.
    δΔϵ_2000, δΔΨ_2000 = deps_dpsi(eop_data, jd_utc)
    δΔϵ_2000 *= milliarcsec_to_rad
    δΔΨ_2000 *= milliarcsec_to_rad

    # Compute the rotation.
    r_tirs_itrf  = r_itrf_to_tirs_iau2006(T, jd_tt, x_p, y_p)
    r_mod_tirs   = r_tirs_to_mod_iau2006(T, jd_ut1, jd_tt, δΔϵ_2000, δΔΨ_2000)
    r_mj2000_mod = r_mod_to_mj2000_iau2006(T, jd_tt)

    return compose_rotation(r_tirs_itrf, r_mod_tirs, r_mj2000_mod)
end

# NOTE: We do not implement the conversion ITRF => GCRF using equinox-based
# IAU-2006/2010 theory because the CIO-based approach is faster and more
# precise.

#                                  TIRS => ERS
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:TIRS},
    ::Val{:ERS},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    milliarcsec_to_rad = π / 648000000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Obtain the correction of the nutation in longitude.
    ~, δΔΨ_2000 = deps_dpsi(eop_data, jd_utc)
    δΔΨ_2000 *= milliarcsec_to_rad

    # Compute the rotation.
    return r_tirs_to_ers_iau2006(T, jd_ut1, jd_tt, δΔΨ_2000)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:TIRS}, ::Val{:ERS}, jd_utc::Number)
    # Since we do not have EOP Data, assume that jd_utc is equal to jd_ut1.
    jd_ut1 = jd_utc
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Compute the rotation.
    return r_tirs_to_ers_iau2006(T, jd_ut1, jd_tt)
end

#                                 TIRS => MOD
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:TIRS},
    ::Val{:MOD06},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    milliarcsec_to_rad = π / 648000000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Obtain the correction of the nutation in obliquity and longitude.
    δΔϵ_2000, δΔΨ_2000 = deps_dpsi(eop_data, jd_utc)
    δΔϵ_2000 *= milliarcsec_to_rad
    δΔΨ_2000 *= milliarcsec_to_rad

    # Compute the rotation.
    return r_tirs_to_mod_iau2006(T, jd_ut1, jd_tt, δΔϵ_2000, δΔΨ_2000)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:TIRS}, ::Val{:MOD06}, jd_utc::Number)
    # Since we do not have EOP Data, assume that jd_utc is equal to jd_ut1.
    jd_ut1 = jd_utc
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Compute the rotation.
    return r_tirs_to_mod_iau2006(T, jd_ut1, jd_tt)
end

#                                TIRS => MJ2000
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:TIRS},
    ::Val{:MJ2000},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    milliarcsec_to_rad = π / 648000000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Obtain the correction of the nutation in obliquity and longitude.
    δΔϵ_2000, δΔΨ_2000 = deps_dpsi(eop_data, jd_utc)
    δΔϵ_2000 *= milliarcsec_to_rad
    δΔΨ_2000 *= milliarcsec_to_rad

    # Compute the rotation.
    r_mod_tirs   = r_tirs_to_mod_iau2006(T, jd_ut1, jd_tt, δΔϵ_2000, δΔΨ_2000)
    r_mj2000_mod = r_mod_to_mj2000_iau2006(T, jd_tt)

    return compose_rotation(r_mod_tirs, r_mj2000_mod)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:TIRS}, ::Val{:MJ2000}, jd_utc::Number)
    # Since we do not have EOP Data, assume that jd_utc is equal to jd_ut1.
    jd_ut1 = jd_utc
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Compute the rotation.
    r_mod_tirs   = r_tirs_to_mod_iau2006(T, jd_ut1, jd_tt)
    r_mj2000_mod = r_mod_to_mj2000_iau2006(T, jd_tt)

    return compose_rotation(r_mod_tirs, r_mj2000_mod)
end

# NOTE: We do not implement the conversion TIRS => GCRF using equinox-based
# IAU-2006/2010 theory because the CIO-based approach is faster and more
# precise.
