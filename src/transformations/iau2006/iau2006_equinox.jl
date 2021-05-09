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
#   [2] Capitaine, N., Wallace, P. T (2006). High precision methods for locating
#       the celestial intermediate pole and origin. Astronomy & Astrophysics.
#
#   [3] IERS (2010). Transformation between the International Terrestrial
#       Reference System and the Geocentric Celestial Reference System. IERS
#       Technical Note No. 36, Chapter 5.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export r_tirs_to_ers_iau2006, r_ers_to_tirs_iau2006,
       r_ers_to_mod_iau2006, r_mod_to_ers_iau2006,
       r_mod_to_mj2000_iau2006, r_mj2000_to_mod_iau2006,
       r_mj2000_to_gcrf_iau2006, r_gcrf_to_mj2000_iau2006,
       r_tirs_to_mod_iau2006, r_mod_to_tirs_iau2006

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
#   - MOD: Mean of Date reference frame.
#   - MJ2000: J2000 mean equatorial frame (mean Equator, mean equinox dynamical
#             system).
#   - GCRF: Geocentric Celestial Reference Frame.
#
# Every rotation will be coded as a function using the IAU-2006 theory with the
# IAU-2010 conventions. The API is:
#
#   function r_<Origin Frame>_to_<Destination Frame>_iau2006
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
    r_tirs_to_ers_iau2006([T::Type,] JD_UT1::Number, JD_TT::Number, δΔΨ_2000::Number = 0)

Compute the rotation that aligns the Terrestrial Intermediate Reference System
(TIRS) with the Earth Reference System (ERS) at the Julian Day
`JD_UT1` [UT1] and `JD_TT` [Terrestrial Time]. This algorithm uses the IAU-2006
theory.

Notice that one can provide corrections for the nutation in longitude
(`δΔψ_2000`) \\[rad] that are usually obtained from IERS EOP Data (see
[`get_iers_eop`](@ref) and [`dEps_dPsi`](@ref)). This corrections are related to
Free Core Nutation (FCN) that models the effect of a liquid Earth core.

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
function r_tirs_to_ers_iau2006(
    JD_UT1::Number,
    JD_TT::Number,
    δΔΨ_2000::Number = 0
)
    return r_tirs_to_ers_iau2006(DCM, JD_UT1, JD_TT)
end

function r_tirs_to_ers_iau2006(
    T::Type,
    JD_UT1::Number,
    JD_TT::Number,
    δΔΨ_2000::Number = 0
)
    # In this theory, the rotation of Earth is taken into account by the Earth
    # Rotation Angle, which is the angle between the Conventional International
    # Origin (CIO) and the Terrestrial Intermediate Origin (TIO) [1]. The latter
    # is a reference meridian on Earth that is located about 100m away from
    # Greenwich meridian along the equator of the Celestial Intermediate Pole
    # (CIP) [1].
    θ_era = 2π * (0.7790572732640 + 1.00273781191135448 * (JD_UT1 - JD_J2000))
    θ_era = mod(θ_era, 2π)

    # Compute the Equation of the Origins (EO).
    ~, ~, ~, EO = nutation_eo_iau2006(JD_TT, 0, δΔΨ_2000)

    # Compute the Greenwich apparent sidereal angle (GAST).
    θ_gast2000 = θ_era - EO

    # Compute the rotation between the TIRS and ERS.
    return angle_to_rot(T, -θ_gast2000, 0, 0, :ZYX)
end

"""
    r_ers_to_tirs_iau2006(JD_UT1::Number, JD_TT::Number, δΔΨ_2000::Number = 0)

Compute the rotation that aligns the Earth Reference System (ERS) with the
Terrestrial Intermediate Reference System (TIRS) at the Julian Day `JD_UT1`
[UT1] and `JD_TT` [Terrestrial Time]. This algorithm uses the IAU-2006 theory.

Notice that one can provide corrections for the nutation in longitude
(`δΔψ_2000`) \\[rad] that are usually obtained from IERS EOP Data (see
[`get_iers_eop`](@ref) and [`dEps_dPsi`](@ref)). This corrections are related to
Free Core Nutation (FCN) that models the effect of a liquid Earth core.

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
function r_ers_to_tirs_iau2006(
    JD_UT1::Number,
    JD_TT::Number,
    δΔΨ_2000::Number = 0
)
    return r_ers_to_tirs_iau2006(DCM, JD_UT1, JD_TT)
end

function r_ers_to_tirs_iau2006(
    T::Type,
    JD_UT1::Number,
    JD_TT::Number,
    δΔΨ_2000::Number = 0
)
    return inv_rotation(r_tirs_to_ers_iau2006(T, JD_UT1, JD_TT, δΔΨ_2000))
end

#                                ERS <=> MOD
# ==============================================================================

"""
    r_ers_to_mod_iau2006([T::Type,] JD_TT::Number, δΔϵ_2000::Number = 0, δΔΨ_2000::Number = 0)

Compute the rotation that aligns the Earth Reference System (ERS) with the
Mean of Date (MOD) reference frame at Julian day `JD_TT` [Terrestrial Time].
This algorithm uses the IAU-2006 theory.

Notice that one can provide corrections for the nutation in obliquity
(`δΔϵ_2000`) and in longitude (`δΔψ_2000`) \\[rad] that are usually obtained
from IERS EOP Data (see [`get_iers_eop`](@ref) and [`dEps_dPsi`](@ref)). This
corrections are related to Free Core Nutation (FCN) that models the effect of a
liquid Earth core.

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the ERS frame with the MOD frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The reference systems ERS and MOD are separated by the nutation of the pole.

"""
function r_ers_to_mod_iau2006(
    JD_TT::Number,
    δΔϵ_2000::Number = 0,
    δΔΨ_2000::Number = 0
)
    return r_ers_to_mod_iau2006(DCM, JD_TT, δΔϵ_2000, δΔΨ_2000)
end

function r_ers_to_mod_iau2006(
    T::Type,
    JD_TT::Number,
    δΔϵ_2000::Number = 0,
    δΔΨ_2000::Number = 0
)
    # Compute the angles used to compute the nutation.
    mϵ_2000, Δϵ_2000, ΔΨ_2000, ~ = nutation_eo_iau2006(JD_TT, δΔϵ_2000, δΔΨ_2000)
    return angle_to_rot(T, mϵ_2000 + Δϵ_2000, ΔΨ_2000, -mϵ_2000, :XZX)
end

"""
    r_mod_to_ers_iau2006([T::Type,] JD_TT::Number, δΔϵ_2000::Number = 0, δΔΨ_2000::Number = 0)

Compute the rotation that aligns the Mean of Date (MOD) reference frame with the
Earth Reference System (ERS) at Julian day `JD_TT` [Terrestrial Time]. This
algorithm uses the IAU-2006 theory.

Notice that one can provide corrections for the nutation in obliquity
(`δΔϵ_2000`) and in longitude (`δΔψ_2000`) \\[rad] that are usually obtained
from IERS EOP Data (see [`get_iers_eop`](@ref) and [`dEps_dPsi`](@ref)). This
corrections are related to Free Core Nutation (FCN) that models the effect of a
liquid Earth core.

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the MOD frame with the ERS frame. The rotation
representation is selected by the optional parameter `T`.

"""
function r_mod_to_ers_iau2006(
    JD_TT::Number,
    δΔϵ_2000::Number = 0,
    δΔΨ_2000::Number = 0
)
    return r_mod_to_ers_iau2006(DCM, JD_TT, δΔϵ_2000, δΔΨ_2000)
end

function r_mod_to_ers_iau2006(
    T::Type,
    JD_TT::Number,
    δΔϵ_2000::Number = 0,
    δΔΨ_2000::Number = 0
)
    return inv_rotation(r_ers_to_mod_iau2006(T, JD_TT, δΔϵ_2000, δΔΨ_2000))
end

#                                MOD <=> MJ2000
# ==============================================================================

"""
    r_mod_to_mj2000_iau2006([T::Type,] JD_TT::Number)

Compute the rotation that aligns the Mean of Date (MOD) reference frame with the
J2000 mean equatorial frame at Julian day `JD_TT` [Terrestrial Time]. This
algorithm uses the IAU-2006 theory.

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the MOD frame with the MJ2000 frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The J2000 reference frame here is not equal to the previous definition in FK5
theory. It is the reason why it is internally called `MJ2000`. According to [3]:

> The mean equinox of J2000.0 to be considered is not the “rotational dynamical
> mean equinox of J2000.0” as used in the past, but the “inertial dynamical mean
> equinox of J2000.0” to which the recent numerical or analytical solutions
> refer.  The latter is associated with the ecliptic in the inertial sense,
> which is the plane perpendicular to the angular momentum vector of the orbital
> motion of the Earth-Moon barycenter as computed from the velocity of the
> barycenter relative to an inertial system. The rotational equinox is
> associated with the ecliptic in the rotational sense, which is perpendicular
> to the angular momentum vector computed from the velocity referred to the
> rotating orbital plane of the Earth-Moon barycenter. (The difference between
> the two angular momenta is the angular momentum associated with the rotation
> of the orbital plane.)

"""
r_mod_to_mj2000_iau2006(JD_TT::Number) = r_mod_to_mj2000_iau2006(DCM, JD_TT)

function r_mod_to_mj2000_iau2006(T::Type, JD_TT::Number)
    # Compute the angles used in the precession model.
    Ψ_a, ω_a, χ_a = precession_iau2006(JD_TT)
    ϵ_0 = 84381.406 * π / 648_000

    # NOTE: According to [2], the matrix as written in [1, p. 218] rotates the
    # MJ2000 to the MOD. Hence, we need the inverse matrix. Furthermore, the
    # equation in [1, p. 218, eq. 3-73] uses mϵ_2000 instead of ϵ_0 as in
    # [2, eq. 12].
    return compose_rotation(
        angle_to_rot(T, -χ_a, ω_a, Ψ_a, :ZXZ),
        angle_to_rot(T, -ϵ_0, 0, 0, :XYZ)
    )
end

"""
    r_mj2000_to_mod_iau2006([T::Type,] JD_TT::Number)

Compute the rotation that aligns the J2000 mean equatorial frame with the Mean
of Date (MOD) reference frame with the at Julian day `JD_TT` [Terrestrial Time].
This algorithm uses the IAU-2006 theory.

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the MJ2000 frame with the MOD frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The J2000 reference frame here is not equal to the previous definition in FK5
theory. It is the reason why it is internally called `MJ2000`. According to [3]:

> The mean equinox of J2000.0 to be considered is not the “rotational dynamical
> mean equinox of J2000.0” as used in the past, but the “inertial dynamical mean
> equinox of J2000.0” to which the recent numerical or analytical solutions
> refer.  The latter is associated with the ecliptic in the inertial sense,
> which is the plane perpendicular to the angular momentum vector of the orbital
> motion of the Earth-Moon barycenter as computed from the velocity of the
> barycenter relative to an inertial system. The rotational equinox is
> associated with the ecliptic in the rotational sense, which is perpendicular
> to the angular momentum vector computed from the velocity referred to the
> rotating orbital plane of the Earth-Moon barycenter. (The difference between
> the two angular momenta is the angular momentum associated with the rotation
> of the orbital plane.)

"""
r_mj2000_to_mod_iau2006(JD_TT::Number) = r_mj2000_to_mod_iau2006(DCM, JD_TT)

function r_mj2000_to_mod_iau2006(T::Type, JD_TT::Number)
    return inv_rotation(r_mod_to_mj2000_iau2006(T, JD_TT))
end

################################################################################
#                               MJ2000 <=> GCRF
################################################################################

"""
    r_mj2000_to_gcrf_iau2006([T::Type,] JD_TT::Number = 0)

Compute the rotation that aligns the J2000 mean equatorial frame with the
Geocentric Celestial Reference Frame (GCRF). This algorithm uses the IAU-2006
theory. Notice that this rotation is just a bias matrix that does not depend on
the date. However, this function receives the argument `JD_TT` just to keep the
API compatibility.

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the MJ2000 frame with the MOD frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

According to [1], the frame bias that converts MJ2000 <=> GCRF is not a precise
transformation for all the times.

"""
r_mj2000_to_gcrf_iau2006(JD_TT::Number = 0) = r_mj2000_to_gcrf_iau2006(DCM, JD_TT)

function r_mj2000_to_gcrf_iau2006(T::Type, JD_TT::Number = 0)
    # Auxiliary variables.
    d2r = π / 180
    a2d = 1 / 3600
    a2r = a2d * d2r

    δα₀ = -0.0146 * a2r
    ξ₀  = -0.041775 * sin(84381.448 * a2r) * a2r
    η₀  = -0.0068192 * a2r

    return angle_to_rot(T, η₀, -ξ₀, -δα₀, :XYZ)
end

"""
    r_gcrf_to_mj2000_iau2006([T::Type,] JD_TT::Number = 0)

Compute the rotation that aligns the Geocentric Celestial Reference Frame (GCRF)
with the J2000 mean equatorial frame. This algorithm uses the IAU-2006 theory.
Notice that this rotation is just a bias matrix that does not depend on the
date. However, this function receives the argument `JD_TT` just to keep the API
compatibility.

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the MJ2000 frame with the MOD frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

According to [1], the frame bias that converts MJ2000 <=> GCRF is not a precise
transformation for all the times.

"""
r_gcrf_to_mj2000_iau2006(JD_TT::Number = 0) = r_gcrf_to_mj2000_iau2006(DCM, JD_TT)

function r_gcrf_to_mj2000_iau2006(T::Type, JD_TT::Number = 0)
    return inv_rotation(r_mj2000_to_gcrf_iau2006(T, JD_TT))
end

################################################################################
#                              Multiple Rotations
################################################################################

# The functions with multiple rotations must be added here only when the it will
# decrease the computational burden compared to calling the functions with the
# single rotations.

"""
    r_tirs_to_mod_iau2006([T::Type,] JD_UT1::Number, JD_TT::Number, δΔϵ_2000::Number = 0, δΔΨ_2000::Number = 0)

Compute the rotation that aligns the Terrestrial Intermediate Reference System
(TIRS) with the Mean of Date (MOD) reference frame at the Julian Day
`JD_UT1` [UT1] and `JD_TT` [Terrestrial Time]. This algorithm uses the IAU-2006
theory.

Notice that one can provide corrections for the nutation in obliquity
(`δΔϵ_2000`) and in longitude (`δΔψ_2000`) \\[rad] that are usually obtained
from IERS EOP Data (see [`get_iers_eop`](@ref) and [`dEps_dPsi`](@ref)). This
corrections are related to Free Core Nutation (FCN) that models the effect of a
liquid Earth core.

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TIRS frame with the ERS frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

This composed rotation TIRS <=> ERS <=> MOD is implemented as a new function
because the single rotations TIRS <=> ERS and ERS <=> MOD call the function
`nutation_eo`, which has a high computational burden. In this case, the composed
algorithm is about 2x faster than calling those function separately.

"""
function r_tirs_to_mod_iau2006(
    JD_UT1::Number,
    JD_TT::Number,
    δΔϵ_2000::Number = 0,
    δΔΨ_2000::Number = 0
)
    return r_tirs_to_mod_iau2006(DCM, JD_UT1, JD_TT, δΔϵ_2000, δΔΨ_2000)
end

function r_tirs_to_mod_iau2006(
    T::Type,
    JD_UT1::Number,
    JD_TT::Number,
    δΔϵ_2000::Number = 0,
    δΔΨ_2000::Number = 0
)
    # In this theory, the rotation of Earth is taken into account by the Earth
    # Rotation Angle, which is the angle between the Conventional International
    # Origin (CIO) and the Terrestrial Intermediate Origin (TIO) [1]. The latter
    # is a reference meridian on Earth that is located about 100m away from
    # Greenwich meridian along the equator of the Celestial Intermediate Pole
    # (CIP) [1].
    θ_era = 2π * (0.7790572732640 + 1.00273781191135448 * (JD_UT1 - JD_J2000))
    θ_era = mod(θ_era, 2π)

    # Compute the Equation of the Origins (EO).
    mϵ_2000, Δϵ_2000, ΔΨ_2000, EO = nutation_eo_iau2006(JD_TT, δΔϵ_2000, δΔΨ_2000)

    # Compute the Greenwich apparent sidereal angle (GAST).
    θ_gast2000 = θ_era - EO

    # Compute the rotation between the TIRS and ERS.
    r_ERS_TIRS = angle_to_rot(T, -θ_gast2000, 0, 0, :ZYX)

    # Compute the rotation between ERS and MOD.
    r_MOD_ERS = angle_to_rot(T, mϵ_2000 + Δϵ_2000, ΔΨ_2000, -mϵ_2000, :XZX)

    # Return the composed rotation.
    return compose_rotation(r_ERS_TIRS, r_MOD_ERS)
end

"""
    r_mod_to_tirs_iau2006([T::Type,] JD_UT1::Number, JD_TT::Number, δΔϵ_2000::Number = 0, δΔΨ_2000::Number = 0)

Compute the rotation that aligns the Mean of Date (MOD) reference frame with the
Terrestrial Intermediate Reference System (TIRS) at the Julian Day `JD_UT1`
[UT1] and `JD_TT` [Terrestrial Time]. This algorithm uses the IAU-2006
theory.

Notice that one can provide corrections for the nutation in obliquity
(`δΔϵ_2000`) and in longitude (`δΔψ_2000`) \\[rad] that are usually obtained
from IERS EOP Data (see [`get_iers_eop`](@ref) and [`dEps_dPsi`](@ref)). This
corrections are related to Free Core Nutation (FCN) that models the effect of a
liquid Earth core.

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TIRS frame with the ERS frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

This composed rotation TIRS <=> ERS <=> MOD is implemented as a new function
because the single rotations TIRS <=> ERS and ERS <=> MOD call the function
`nutation_eo`, which has a high computational burden. In this case, the composed
algorithm is about 2x faster than calling those function separately.

"""
function r_mod_to_tirs_iau2006(
    JD_UT1::Number,
    JD_TT::Number,
    δΔϵ_2000::Number = 0,
    δΔΨ_2000::Number = 0
)
    return r_mod_to_tirs_iau2006(DCM, JD_UT1, JD_TT, δΔϵ_2000, δΔΨ_2000)
end

function r_mod_to_tirs_iau2006(
    T::Type,
    JD_UT1::Number,
    JD_TT::Number,
    δΔϵ_2000::Number = 0,
    δΔΨ_2000::Number = 0
)
    return inv_rotation(r_tirs_to_mod_iau2006(T, JD_UT1, JD_TT, δΔϵ_2000, δΔΨ_2000))
end
