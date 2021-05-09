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
    r_ecef_to_eci([T,] ECEF, ECI, JD_UTC::Number [, eop_data])

Compute the rotation from an Earth-Centered, Earth-Fixed (`ECEF`) reference
frame to an Earth-Centered Inertial (`ECI`) reference frame at the Julian Day
[UTC] `JD_UTC`. The rotation description that will be used is given by `T`,
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

* `DCM`: The rotation will be described by a Direction Cosine Matrix.
* `Quaternion`: The rotation will be described by a Quaternion.

If no value is specified, then it falls back to `DCM`.

# Conversion model

The model that will be used to compute the rotation is automatically inferred
given the selection of the origin and destination frames. **Notice that mixing
IAU-76/FK5 and IAU-2006/2010 frames is not supported yet.**

# ECEF Frame

The ECEF frame is selected by the parameter `ECEF`. The possible values are:

* `ITRF()`: ECEF will be selected as the International Terrestrial Reference
            Frame (ITRF).
* `PEF()`: ECEF will be selected as the Pseudo-Earth Fixed (PEF) reference
           frame.
* `TIRS()`: ECEF will be selected as the Terrestrial Intermediate Reference
            System (TIRS).

# ECI Frame

The ECI frame is selected by the parameter `ECI`. The possible values are:

* `TEME()`: ECI will be selected as the True Equator Mean Equinox (TEME)
            reference frame.
* `TOD()`: ECI will be selected as the True of Date (TOD).
* `MOD()`: ECI will be selected as the Mean of Date (MOD).
* `J2000()`: ECI will be selected as the J2000 reference frame.
* `GCRF()`: ECI will be selected as the Geocentric Celestial Reference Frame
            (GCRF).
* `CIRS()`: ECI will be selected as the Celestial Intermediate Reference System
            (CIRS).
* `ERS()`: ECI will be selected as the Earth Reference System (ERS).
* `MOD06()`: ECI will be selected as the Mean of Date (MOD) according to the
             definition in IAU-2006/2010 theory.
* `MJ2000()`: ECI will be selected as the J2000 mean equatorial frame (MJ2000).

!!! note

    The frames `MOD()` and `MOD06()` are virtually the same. However, we
    selected different names to make clear which theory are being used since
    mixing transformation between frames from IAU-76/FK5 and IAU-2006/2010 must
    be performed with caution.

# EOP Data

The conversion between the frames depends on EOP Data (see `get_iers_eop` and
`read_iers_eop`). If IAU-76/FK5 model is used, then the type of `eop_data` must
be `EOPData_IAU1980`. Otherwise, if IAU-2006/2010 model is used, then the type
of `eop_data` must be `EOPData_IAU2000A`. The following table shows the
requirements for EOP data given the selected frames.

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

`¹`: In this case, the Julian Time UTC will be assumed equal to Julian Time UT1
to compute the Greenwich Mean Sidereal Time. This is an approximation, but
should be sufficiently accurate for some applications. Notice that, if EOP Data
is provided, the Julian Day UT1 will be accurately computed.

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
julia> eop_IAU1980 = get_iers_eop(:IAU1980);

julia> r_ecef_to_eci(DCM, ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267      0.78518     -0.00132979
 -0.78518      -0.619267     3.33492e-5
 -0.000797313   0.00106478   0.999999

julia> r_ecef_to_eci(ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267      0.78518     -0.00132979
 -0.78518      -0.619267     3.33492e-5
 -0.000797313   0.00106478   0.999999

julia> r_ecef_to_eci(PEF(), J2000(), DatetoJD(1986, 06, 19, 21, 35, 0))
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619271      0.785176    -0.00133066
 -0.785177     -0.619272     3.45854e-5
 -0.000796885   0.00106622   0.999999

julia> r_ecef_to_eci(PEF(), J2000(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267      0.78518     -0.00133066
 -0.78518      -0.619267     3.45854e-5
 -0.000796879   0.00106623   0.999999

julia> r_ecef_to_eci(Quaternion, ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
Quaternion{Float64}:
  + 0.4363098936462618 - 0.0005909969666939257.i + 0.00030510511316206974.j + 0.8997962182293519.k

julia> eop_IAU2000A = get_iers_eop(:IAU2000A);

julia> r_ecef_to_eci(ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU2000A)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267      0.78518     -0.00132979
 -0.78518      -0.619267     3.33502e-5
 -0.000797312   0.00106478   0.999999

julia> r_ecef_to_eci(TIRS(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0))
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619271      0.785176    -0.00133066
 -0.785177     -0.619272     3.45884e-5
 -0.000796885   0.00106623   0.999999

julia> r_ecef_to_eci(Quaternion, ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU2000A)
Quaternion{Float64}:
  + 0.4363098936309669 - 0.000590996988144556.i + 0.0003051056555230158.j + 0.8997962182365703.k
```
"""
@inline function r_ecef_to_eci(
    T_ECEF::T_ECEFs,
    T_ECI::T_ECIs,
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    r_ecef_to_eci(DCM, T_ECEF, T_ECI, JD_UTC, eop_data)
end

@inline function r_ecef_to_eci(
    T_ECEF::T_ECEFs_IAU_2006,
    T_ECI::T_ECIs_IAU_2006,
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    r_ecef_to_eci(DCM, T_ECEF, T_ECI, JD_UTC, eop_data)
end

# Specializations for those cases that EOP Data is not needed.
@inline function r_ecef_to_eci(
    T_ECEF::Val{:PEF},
    T_ECI::Union{Val{:J2000}, Val{:MOD}, Val{:TOD}, Val{:TEME}},
    JD_UTC::Number
)
    r_ecef_to_eci(DCM, T_ECEF, T_ECI, JD_UTC)
end

@inline function r_ecef_to_eci(
        T_ECEF::Val{:TIRS},
        T_ECI::T_ECIs_IAU_2006,
        JD_UTC::Number
)
    r_ecef_to_eci(DCM, T_ECEF, T_ECI, JD_UTC)
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
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x(JD_UTC)*arcsec2rad
    y_p      = eop_data.y(JD_UTC)*arcsec2rad
    δΔϵ_1980 = eop_data.dEps(JD_UTC)*arcsec2rad
    δΔψ_1980 = eop_data.dPsi(JD_UTC)*arcsec2rad

    # Compute the rotation.
    return r_itrf_to_gcrf_fk5(T, JD_UT1, JD_TT, x_p, y_p, δΔϵ_1980, δΔψ_1980)
end

#                                ITRF => J2000
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:J2000},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x(JD_UTC)*arcsec2rad
    y_p      = eop_data.y(JD_UTC)*arcsec2rad

    # Compute the rotation.
    return r_itrf_to_gcrf_fk5(T, JD_UT1, JD_TT, x_p, y_p, 0, 0)
end

#                                 ITRF => MOD
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:MOD},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x(JD_UTC)*arcsec2rad
    y_p      = eop_data.y(JD_UTC)*arcsec2rad
    δΔϵ_1980 = eop_data.dEps(JD_UTC)*arcsec2rad
    δΔψ_1980 = eop_data.dPsi(JD_UTC)*arcsec2rad

    # Compute the rotation.
    r_PEF_ITRF = r_itrf_to_pef_fk5(T, x_p, y_p)
    r_MOD_PEF  = r_pef_to_mod_fk5(T, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)

    return compose_rotation(r_PEF_ITRF, r_MOD_PEF)
end

#                                 ITRF => TOD
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:TOD},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x(JD_UTC)*arcsec2rad
    y_p      = eop_data.y(JD_UTC)*arcsec2rad
    δΔψ_1980 = eop_data.dPsi(JD_UTC)*arcsec2rad

    # Compute the rotation.
    r_PEF_ITRF = r_itrf_to_pef_fk5(T, x_p, y_p)
    r_TOD_PEF  = r_pef_to_tod_fk5(T, JD_UT1, JD_TT, δΔψ_1980)

    return compose_rotation(r_PEF_ITRF, r_TOD_PEF)
end

#                                 ITRF => TEME
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:TEME},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x(JD_UTC)*arcsec2rad
    y_p      = eop_data.y(JD_UTC)*arcsec2rad

    # Compute the rotation.
    r_PEF_ITRF = r_itrf_to_pef_fk5(T, x_p, y_p)
    r_TEME_PEF = rPEFtoTEME(T, JD_UT1)

    return compose_rotation(r_PEF_ITRF, r_TEME_PEF)
end

#                                 PEF => GCRF
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:PEF},
    ::Val{:GCRF},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x(JD_UTC)*arcsec2rad
    y_p      = eop_data.y(JD_UTC)*arcsec2rad
    δΔϵ_1980 = eop_data.dEps(JD_UTC)*arcsec2rad
    δΔψ_1980 = eop_data.dPsi(JD_UTC)*arcsec2rad

    # Compute the rotation.
    r_MOD_PEF  = r_pef_to_mod_fk5(T, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)
    r_GCRF_MOD = r_mod_to_gcrf_fk5(T, JD_TT)

    return compose_rotation(r_MOD_PEF, r_GCRF_MOD)
end

#                                PEF => J2000
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:PEF},
    ::Val{:J2000},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Compute the rotation.
    r_MOD_PEF  = r_pef_to_mod_fk5(T, JD_UT1, JD_TT, 0, 0)
    r_GCRF_MOD = r_mod_to_gcrf_fk5(T, JD_TT)

    return compose_rotation(r_MOD_PEF, r_GCRF_MOD)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:PEF}, ::Val{:J2000}, JD_UTC::Number)
    # Since we do not have EOP Data, assume that JD_UTC is equal to JD_UT1.
    JD_UT1 = JD_UTC
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Compute the rotation.
    r_MOD_PEF  = r_pef_to_mod_fk5(T, JD_UT1, JD_TT, 0, 0)
    r_GCRF_MOD = r_mod_to_gcrf_fk5(T, JD_TT)

    return compose_rotation(r_MOD_PEF, r_GCRF_MOD)
end

#                                 PEF => MOD
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:PEF},
    ::Val{:MOD},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x(JD_UTC)*arcsec2rad
    y_p      = eop_data.y(JD_UTC)*arcsec2rad
    δΔϵ_1980 = eop_data.dEps(JD_UTC)*arcsec2rad
    δΔψ_1980 = eop_data.dPsi(JD_UTC)*arcsec2rad

    # Compute the rotation.
    return r_pef_to_mod_fk5(T, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:PEF}, ::Val{:MOD}, JD_UTC::Number)
    # Since we do not have EOP Data, assume that JD_UTC is equal to JD_UT1.
    JD_UT1 = JD_UTC
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Compute the rotation.
    return r_pef_to_mod_fk5(T, JD_UT1, JD_TT, 0, 0)
end

#                                 PEF => TOD
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:PEF},
    ::Val{:TOD},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x(JD_UTC)*arcsec2rad
    y_p      = eop_data.y(JD_UTC)*arcsec2rad
    δΔψ_1980 = eop_data.dPsi(JD_UTC)*arcsec2rad

    # Compute the rotation.
    return r_pef_to_tod_fk5(T, JD_UT1, JD_TT, δΔψ_1980)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:PEF}, ::Val{:TOD}, JD_UTC::Number)

    # Since we do not have EOP Data, assume that JD_UTC is equal to JD_UT1.
    JD_UT1 = JD_UTC
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Compute the rotation.
    return r_pef_to_tod_fk5(T, JD_UT1, JD_TT, 0)
end

#                                 PEF => TEME
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:PEF},
    ::Val{:TEME},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    # Get the time in UT1.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)

    # Compute the rotation.
    return rPEFtoTEME(T, JD_UT1)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:PEF}, ::Val{:TEME}, JD_UTC::Number)
    # Since we do not have EOP Data, assume that JD_UTC is equal to JD_UT1.
    JD_UT1 = JD_UTC

    # Compute the rotation.
    return rPEFtoTEME(T, JD_UT1)
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
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    x_p = eop_data.x(JD_UTC)*arcsec2rad
    y_p = eop_data.y(JD_UTC)*arcsec2rad

    # Compute the rotation.
    r_TIRS_ITRF = rITRFtoTIRS_iau2006(T, JD_TT, x_p, y_p)
    r_CIRS_TIRS = rTIRStoCIRS_iau2006(T, JD_UT1)

    return compose_rotation(r_TIRS_ITRF, r_CIRS_TIRS)
end

#                                 ITRF => GCRF
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:GCRF},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    x_p = eop_data.x(JD_UTC)*arcsec2rad
    y_p = eop_data.y(JD_UTC)*arcsec2rad
    dX  = eop_data.dX(JD_UTC)*arcsec2rad
    dY  = eop_data.dY(JD_UTC)*arcsec2rad

    # Compute the rotation.
    r_TIRS_ITRF = rITRFtoTIRS_iau2006(T, JD_TT, x_p, y_p)
    r_CIRS_TIRS = rTIRStoCIRS_iau2006(T, JD_UT1)
    r_GCRF_CIRS = rCIRStoGCRF_iau2006(T, JD_TT, dX, dY)

    return compose_rotation(r_TIRS_ITRF, r_CIRS_TIRS, r_GCRF_CIRS)
end

#                                 TIRS => CIRS
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:TIRS},
    ::Val{:CIRS},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)

    # Compute the rotation.
    return rTIRStoCIRS_iau2006(T, JD_UT1)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:TIRS}, ::Val{:CIRS}, JD_UTC::Number)
    # Since we do not have EOP Data, assume that JD_UTC is equal to JD_UT1.
    JD_UT1 = JD_UTC

    # Compute the rotation.
    return rTIRStoCIRS_iau2006(T, JD_UT1)
end

#                                 TIRS => GCRF
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:TIRS},
    ::Val{:GCRF},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    dX  = eop_data.dX(JD_UTC)*arcsec2rad
    dY  = eop_data.dY(JD_UTC)*arcsec2rad

    # Compute the rotation.
    r_CIRS_TIRS = rTIRStoCIRS_iau2006(T, JD_UT1)
    r_GCRF_CIRS = rCIRStoGCRF_iau2006(T, JD_TT, dX, dY)

    return compose_rotation(r_CIRS_TIRS, r_GCRF_CIRS)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:TIRS}, ::Val{:GCRF}, JD_UTC::Number)
    # Since we do not have EOP Data, assume that JD_UTC is equal to JD_UT1.
    JD_UT1 = JD_UTC

    # Get the time in TT.
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Compute the rotation.
    r_CIRS_TIRS = rTIRStoCIRS_iau2006(T, JD_UT1)
    r_GCRF_CIRS = rCIRStoGCRF_iau2006(T, JD_TT)

    return compose_rotation(r_CIRS_TIRS, r_GCRF_CIRS)
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
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    x_p = eop_data.x(JD_UTC)*arcsec2rad
    y_p = eop_data.y(JD_UTC)*arcsec2rad

    # Obtain the correction of the nutation in longitude.
    ~, δΔΨ_2000 = dEps_dPsi(eop_data, JD_UTC)
    δΔΨ_2000 *= arcsec2rad

    # Compute the rotation.
    r_TIRS_ITRF = rITRFtoTIRS_iau2006(T, JD_TT, x_p, y_p)
    r_ERS_TIRS = rTIRStoERS_iau2006(T, JD_UT1, JD_TT, δΔΨ_2000)

    return compose_rotation(r_TIRS_ITRF, r_ERS_TIRS)
end

#                                 ITRF => MOD
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:MOD06},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    x_p = eop_data.x(JD_UTC)*arcsec2rad
    y_p = eop_data.y(JD_UTC)*arcsec2rad

    # Obtain the correction of the nutation in obliquity and longitude.
    δΔϵ_2000, δΔΨ_2000 = dEps_dPsi(eop_data, JD_UTC)
    δΔϵ_2000 *= arcsec2rad
    δΔΨ_2000 *= arcsec2rad

    # Compute the rotation.
    r_TIRS_ITRF = rITRFtoTIRS_iau2006(T, JD_TT, x_p, y_p)
    r_MOD_TIRS = rTIRStoMOD_iau2006(T, JD_UT1, JD_TT, δΔϵ_2000, δΔΨ_2000)

    return compose_rotation(r_TIRS_ITRF, r_MOD_TIRS)
end

#                                ITRF => MJ2000
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:ITRF},
    ::Val{:MJ2000},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    x_p = eop_data.x(JD_UTC)*arcsec2rad
    y_p = eop_data.y(JD_UTC)*arcsec2rad

    # Obtain the correction of the nutation in obliquity and longitude.
    δΔϵ_2000, δΔΨ_2000 = dEps_dPsi(eop_data, JD_UTC)
    δΔϵ_2000 *= arcsec2rad
    δΔΨ_2000 *= arcsec2rad

    # Compute the rotation.
    r_TIRS_ITRF  = rITRFtoTIRS_iau2006(T, JD_TT, x_p, y_p)
    r_MOD_TIRS   = rTIRStoMOD_iau2006(T, JD_UT1, JD_TT, δΔϵ_2000, δΔΨ_2000)
    r_MJ2000_MOD = rMODtoMJ2000_iau2006(T, JD_TT)

    return compose_rotation(r_TIRS_ITRF, r_MOD_TIRS, r_MJ2000_MOD)
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
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Obtain the correction of the nutation in longitude.
    ~, δΔΨ_2000 = dEps_dPsi(eop_data, JD_UTC)
    δΔΨ_2000 *= arcsec2rad

    # Compute the rotation.
    return rTIRStoERS_iau2006(T, JD_UT1, JD_TT, δΔΨ_2000)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:TIRS}, ::Val{:ERS}, JD_UTC::Number)
    # Since we do not have EOP Data, assume that JD_UTC is equal to JD_UT1.
    JD_UT1 = JD_UTC
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Compute the rotation.
    return rTIRStoERS_iau2006(T, JD_UT1, JD_TT)
end

#                                 TIRS => MOD
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:TIRS},
    ::Val{:MOD06},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Obtain the correction of the nutation in obliquity and longitude.
    δΔϵ_2000, δΔΨ_2000 = dEps_dPsi(eop_data, JD_UTC)
    δΔϵ_2000 *= arcsec2rad
    δΔΨ_2000 *= arcsec2rad

    # Compute the rotation.
    return rTIRStoMOD_iau2006(T, JD_UT1, JD_TT, δΔϵ_2000, δΔΨ_2000)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:TIRS}, ::Val{:MOD06}, JD_UTC::Number)
    # Since we do not have EOP Data, assume that JD_UTC is equal to JD_UT1.
    JD_UT1 = JD_UTC
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Compute the rotation.
    return rTIRStoMOD_iau2006(T, JD_UT1, JD_TT)
end

#                                TIRS => MJ2000
# ==============================================================================

function r_ecef_to_eci(
    T::T_ROT,
    ::Val{:TIRS},
    ::Val{:MJ2000},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Obtain the correction of the nutation in obliquity and longitude.
    δΔϵ_2000, δΔΨ_2000 = dEps_dPsi(eop_data, JD_UTC)
    δΔϵ_2000 *= arcsec2rad
    δΔΨ_2000 *= arcsec2rad

    # Compute the rotation.
    r_MOD_TIRS   = rTIRStoMOD_iau2006(T, JD_UT1, JD_TT, δΔϵ_2000, δΔΨ_2000)
    r_MJ2000_MOD = rMODtoMJ2000_iau2006(T, JD_TT)

    return compose_rotation(r_MOD_TIRS, r_MJ2000_MOD)
end

function r_ecef_to_eci(T::T_ROT, ::Val{:TIRS}, ::Val{:MJ2000}, JD_UTC::Number)
    # Since we do not have EOP Data, assume that JD_UTC is equal to JD_UT1.
    JD_UT1 = JD_UTC
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Compute the rotation.
    r_MOD_TIRS   = rTIRStoMOD_iau2006(T, JD_UT1, JD_TT)
    r_MJ2000_MOD = rMODtoMJ2000_iau2006(T, JD_TT)

    return compose_rotation(r_MOD_TIRS, r_MJ2000_MOD)
end

# NOTE: We do not implement the conversion TIRS => GCRF using equinox-based
# IAU-2006/2010 theory because the CIO-based approach is faster and more
# precise.
