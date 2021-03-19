# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions related with the equinox-based model IAU-2006 with IAU-2010
#   conventions.
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

export rTIRStoERS_iau2006, rERStoTIRS_iau2006

################################################################################
#                      IAU-2006 equinox-based reductions
################################################################################
#
# The conversion between the Geocentric Celestial Reference Frame (GCRF) to the
# International Terrestrial Reference Frame (ITRF) is done by means of:
#
#              GCRF <=> MJ2000 <=> MOD <=> ERS <=> TIRS <=> ITRF
#
# in which:
#   - TIRS: Terrestrial Intermediate Reference System.
#   - ERS: Earth Reference System.
#   - MOD: Mean of Date.
#   - MJ2000: Mean of J2000 (Mean Equator, Mean Equinox dynamical system).
#   - GCRF: Geocentric Celestial Reference Frame.
#
# Every rotation will be coded as a function using the IAU-2006 theory with the
# IAU-2010 conventions. The API is:
#
#   function r<Origin Frame>to<Destination Frame>_iau2006
#
# The arguments vary depending on the origin and destination frame and should be
# verified using the function documentation.
#
################################################################################
#
################################################################################
#                               Single Rotations
################################################################################

#                                 TIRS <=> ERS
# ==============================================================================

"""
    rTIRStoERS_iau2006([T::Type,] JD_UT1::Number, JD_TT::Number)

Compute the rotation that aligns the Terrestrial Intermediate Reference System
(TIRS) with the Earth Reference System (ERS) at the Julian Day
`JD_UT1` [UT1] and `JD_TT` [Terrestrial Time]. This algorithm uses the IAU-2006
theory. Notice that one can provide corrections for the nutation longitude
(`δΔψ_2000`) \\[rad] that are usually obtained from IERS EOP
Data (see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TIRS frame with the ERS frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The reference frames TIRS and ERS are separated by a rotation about the Z-axis
of the Greenwhich apparent sidereal angle (GAST). This angle is computed using
the IAU-2006 theory, which consist of obtaining the Earth Rotation Angle (ERA)
and subtracting the result of the Equation of Origins (EO).

"""
rTIRStoERS_iau2006(JD_UT1::Number, JD_TT::Number, δΔΨ_2000::Number = 0) =
    rTIRStoERS_iau2006(DCM, JD_UT1, JD_TT)

function rTIRStoERS_iau2006(T::Type, JD_UT1::Number, JD_TT::Number,
                            δΔΨ_2000::Number = 0)
    # In this theory, the rotation of Earth is taken into account by the Earth
    # Rotation Angle, which is the angle between the Conventional International
    # Origin (CIO) and the Terrestrial Intermediate Origin (TIO) [1]. The latter
    # is a reference meridian on Earth that is located about 100m away from
    # Greenwich meridian along the equator of the Celestial Intermediate Pole
    # (CIP) [1].
    θ_era = 2π*(0.7790572732640 + 1.00273781191135448*(JD_UT1 - JD_J2000))
    θ_era = mod(θ_era, 2π)

    # Compute the Equation of the Origins (EO).
    ~, ~, ~, EO = nutation_eo_iau2006(JD_TT, 0, δΔΨ_2000)

    # Compute the Greenwich apparent sidereal angle (GAST).
    θ_gast2000 = θ_era - EO

    # Compute the rotation between the TIRS and ERS.
    return angle_to_rot(T, -θ_gast2000, 0, 0, :ZYX)
end

"""
    rERStoTIRS_iau2006(JD_UT1::Number, JD_TT::Number, δΔΨ_2000::Number = 0)

Compute the rotation that aligns the Earth Reference System (ERS) with the
Terrestrial Intermediate Reference System (TIRS) at the Julian Day `JD_UT1`
[UT1] and `JD_TT` [Terrestrial Time]. This algorithm uses the IAU-2006 theory.
Notice that one can provide corrections for the nutation longitude (`δΔψ_2000`)
\\[rad] that are usually obtained from IERS EOP Data (see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the ERS frame with the TIRS frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The reference frames TIRS and ERS are separated by a rotation about the Z-axis
of the Greenwhich apparent sidereal angle (GAST). This angle is computed using
the IAU-2006 theory, which consist of obtaining the Earth Rotation Angle (ERA)
and subtracting the result of the Equation of Origins (EO).

"""
rERStoTIRS_iau2006(JD_UT1::Number, JD_TT::Number, δΔΨ_2000::Number = 0) =
    rERStoTIRS_iau2006(DCM, JD_UT1, JD_TT)

function rERStoTIRS_iau2006(T::Type, JD_UT1::Number, JD_TT::Number,
                            δΔΨ_2000::Number = 0)
    # In this theory, the rotation of Earth is taken into account by the Earth
    # Rotation Angle, which is the angle between the Conventional International
    # Origin (CIO) and the Terrestrial Intermediate Origin (TIO) [1]. The latter
    # is a reference meridian on Earth that is located about 100m away from
    # Greenwich meridian along the equator of the Celestial Intermediate Pole
    # (CIP) [1].
    θ_era = 2π*(0.7790572732640 + 1.00273781191135448*(JD_UT1 - JD_J2000))
    θ_era = mod(θ_era, 2π)

    # Compute the Equation of the Origins (EO).
    ~, ~, ~, EO = nutation_eo_iau2006(JD_TT)

    # Compute the Greenwich apparent sidereal angle (GAST).
    θ_gast2000 = θ_era - EO

    # Compute the rotation between the TIRS and ERS.
    return angle_to_rot(T, +θ_gast2000, 0, 0, :ZYX)
end
