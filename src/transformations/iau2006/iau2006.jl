# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions related with the CIO-based model IAU-2006 with IAU-2010
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

export rITRFtoTIRS_iau2006, rTIRStoITRF_iau2006
export rTIRStoCIRS_iau2006, rCIRStoTIRS_iau2006
export rCIRStoGCRF_iau2006, rGCRFtoCIRS_iau2006

################################################################################
#                             IAU-2006 Reductions
################################################################################
#
# The conversion between the Geocentric Celestial Reference Frame (GCRF) to the
# International Terrestrial Reference Frame (ITRF) is done by means of:
#
#                       GCRF <=> CIRS <=> TIRS <=> ITRF
#
# in which:
#   - TIRS: Terrestrial Intermediate Reference System.
#   - CIRS: Celestial Intermediate Reference System.
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

#                                ITRF <=> TIRS
# ==============================================================================

"""
    rITRFtoTIRS_iau2006([T::Type,] JD_TT::Number, x_p::Number, y_p::Number)

Compute the rotation that aligns the International Terrestrial Reference Frame
(ITRF) with the Terrestrial Intermediate Reference System (TIRS) considering the
polar motion represented by the angles `x_p` [rad] and `y_p` [rad] that are
obtained from IERS EOP Data (see `get_iers_eop`).

`x_p` is the polar motion displacement about X-axis, which is the IERS Reference
Meridian direction (positive south along the 0˚ longitude meridian). `y_p` is
the polar motion displacement about Y-axis (90˚W or 270˚E meridian).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the ITRF frame with the TIRS frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The ITRF is defined based on the International Reference Pole (IRP), which is
the location of the terrestrial pole agreed by international committees [1]. The
Terrestrial Intermediate Reference Frame (TIRS), on the other hand, is defined
based on the Earth axis of rotation, or the Celestial Intermediate Pole (CIP).
Hence, TIRS XY-plane contains the True Equator. Furthermore, since the recovered
latitude and longitude are sensitive to the CIP, then it should be computed
considering the TIRS frame.

The TIRS and PEF (IAU-76/FK5) are virtually the same reference frame, but
according to [1] it is convenient to separate the names as the exact formulae
differ.

"""
rITRFtoTIRS_iau2006(JD_TT::Number, x_p::Number, y_p::Number) =
    rITRFtoTIRS_iau2006(DCM, JD_TT, x_p, y_p)

function rITRFtoTIRS_iau2006(T::Type, JD_TT::Number, x_p::Number, y_p::Number)
    # Convert Julian days to Julian centuries.
    T_TT = (JD_TT - JD_J2000)/36525

    # Notice that this rotation has an additional one, called `sl`, from the
    # IAU-76/FK5 theory that accounts for the instantaneous prime meridian
    # called TIO locator [1, p. 212].
    sl  = (-0.000047*π/648000)*T_TT # [rad]

    return angle_to_rot(T, y_p, x_p, -sl, :XYZ)
end

"""
    rTIRStoITRF_iau2006([T::Type,] JD_TT::Number, x_p::Number, y_p::Number)

Compute the rotation that aligns the Terrestrial Intermediate Reference System
(TIRS) with the International Terrestrial Reference Frame (ITRF) considering the
polar motion represented by the angles `x_p` [rad] and `y_p` [rad] that are
obtained from IERS EOP Data (see `get_iers_eop`).

`x_p` is the polar motion displacement about X-axis, which is the IERS Reference
Meridian direction (positive south along the 0˚ longitude meridian). `y_p` is
the polar motion displacement about Y-axis (90˚W or 270˚E meridian).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TIRS frame with the ITRF frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The ITRF is defined based on the International Reference Pole (IRP), which is
the location of the terrestrial pole agreed by international committees [1]. The
Terrestrial Intermediate Reference Frame (TIRS), on the other hand, is defined
based on the Earth axis of rotation, or the Celestial Intermediate Pole (CIP).
Hence, TIRS XY-plane contains the True Equator. Furthermore, since the recovered
latitude and longitude are sensitive to the CIP, then it should be computed
considering the TIRS frame.

The TIRS and PEF (IAU-76/FK5) are virtually the same reference frame, but
according to [1] it is convenient to separate the names as the exact formulae
differ.

"""
rTIRStoITRF_iau2006(JD_TT::Number, x_p::Number, y_p::Number) =
    rTIRStoITRF_iau2006(DCM, JD_TT, x_p, y_p)

function rTIRStoITRF_iau2006(T::Type, JD_TT::Number, x_p::Number, y_p::Number)
    # Convert Julian days to Julian centuries.
    T_TT = (JD_TT - JD_J2000)/36525

    # Notice that this rotation has an additional one, called `sl`, from the
    # IAU-76/FK5 theory that accounts for the instantaneous prime meridian
    # called TIO locator [1, p. 212].
    sl  = (-0.000047*π/648000)*T_TT # [rad]

    return angle_to_rot(T, sl, -x_p, -y_p, :ZYX)
end

#                                TIRS <=> CIRS
# ==============================================================================

"""
    rTIRStoCIRS_iau2006([T::Type,] JD_UT1::Number)

Compute the rotation that aligns the Terrestrial Intermediate Reference System
(TIRS) with the Celestial Intermediate Reference System (CIRS) at the Julian Day
`JD_UT1` [UT1]. This algorithm uses the IAU-2006 theory.

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TIRS frame with the CIRS frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The reference frames TIRS and CIRS are separated by a rotation about the Z-axis
of the Earth Rotation Angle, which is the angle between the Conventional
International Origin (CIO) and the Terrestrial Intermediate Origin (TIO) [1].
The latter is a reference meridian on Earth that is located about 100m away from
Greenwich meridian along the equator of the Celestial Intermediate Pole (CIP)
[1].

"""
rTIRStoCIRS_iau2006(JD_UT1::Number) = rTIRStoCIRS_iau2006(DCM, JD_UT1)

function rTIRStoCIRS_iau2006(T::Type, JD_UT1::Number)
    # In this theory, the rotation of Earth is taken into account by the Earth
    # Rotation Angle, which is the angle between the Conventional International
    # Origin (CIO) and the Terrestrial Intermediate Origin (TIO) [1]. The latter
    # is a reference meridian on Earth that is located about 100m away from
    # Greenwich meridian along the equator of the Celestial Intermediate Pole
    # (CIP) [1].
    θ_era = 2π*(0.7790572732640 + 1.00273781191135448*(JD_UT1 - JD_J2000))

    return angle_to_rot(T, -θ_era, 0, 0, :ZXY)
end

"""
    rCIRStoTIRS_iau2006([T::Type,] JD_UT1::Number)

Compute the rotation that aligns the Celestial Intermediate Reference System
(CIRS) with the Terrestrial Intermediate Reference System (TIRS) at the Julian
Day `JD_UT1` [UT1]. This algorithm uses the IAU-2006 theory.

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the CIRS frame with the TIRS frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The reference frames TIRS and CIRS are separated by a rotation about the Z-axis
of the Earth Rotation Angle, which is the angle between the Conventional
International Origin (CIO) and the Terrestrial Intermediate Origin (TIO) [1].
The latter is a reference meridian on Earth that is located about 100m away from
Greenwich meridian along the equator of the Celestial Intermediate Pole (CIP)
[1].

"""
rCIRStoTIRS_iau2006(JD_UT1::Number) = rCIRStoTIRS_iau2006(DCM, JD_UT1)

function rCIRStoTIRS_iau2006(T::Type, JD_UT1::Number)
    # In this theory, the rotation of Earth is taken into account by the Earth
    # Rotation Angle, which is the angle between the Conventional International
    # Origin (CIO) and the Terrestrial Intermediate Origin (TIO) [1]. The latter
    # is a reference meridian on Earth that is located about 100m away from
    # Greenwich meridian along the equator of the Celestial Intermediate Pole
    # (CIP) [1].
    θ_era = 2π*(0.7790572732640 + 1.00273781191135448*(JD_UT1 - JD_J2000))

    return angle_to_rot(T, θ_era, 0, 0, :ZXY)
end

#                                CIRS <=> GCRF
# ==============================================================================

"""
    rCIRStoGCRF_iau2006([T::Type,] JD_TT::Number, dX::Number = 0, dY::Number = 0)

Compute the rotation that aligns the Celestial Intermediate Reference System
(CIRS) with the Geocentric Celestial Reference Frame (GCRF) at the Julian Day
`JD_TT` [TT] and considering the IERS EOP Data `dX` [rad] and `dY` [rad] \\(see
`get_iers_eop`). This algorithm uses the IAU-2006 theory.

The IERS EOP Data `dX` and `dY` accounts for the free-core nutation and time
dependent effects of the Celestial Intermediate Pole (CIP) position with respect
to the GCRF.

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the CIRS frame with the GCRF frame. The rotation
representation is selected by the optional parameter `T`.

"""
rCIRStoGCRF_iau2006(JD_TT::Number, dX::Number = 0, dY::Number = 0) =
    rCIRStoGCRF_iau2006(DCM, JD_TT, dX, dY)

function rCIRStoGCRF_iau2006(::Type{DCM}, JD_TT::Number, dX::Number = 0,
                             dY::Number = 0)

    # Compute the CIP position w.r.t. GCRS.
    X, Y, s = precession_nutation_iau2006(JD_TT)

    # Add the corrections.
    X += dX
    Y += dY

    # Auxiliary variables.
    X² = X^2
    Y² = Y^2
    XY = X*Y

    # Compute the rotation matrix
    # ==========================================================================

    # This is the approximate value for:
    #
    #   a = 1/(1 + cos(d)), d = atan( sqrt( ( X^2 + Y^2 )/( 1 - X^2 - Y^2 ) ) )
    a = 1/2 + 1/8*(X² + Y²)

    D = DCM(1 - a*X²,    -a*XY, X,
               -a*XY, 1 - a*Y², Y,
                 -X ,      -Y , 1 - a*(X² + Y²))'

    return D*create_rotation_matrix(s, :Z)
end

rCIRStoGCRF_iau2006(::Type{Quaternion}, JD_TT::Number, dX::Number = 0,
                    dY::Number = 0) =
    dcm_to_quat(rCIRStoGCRF_iau2006(DCM, JD_TT, dX, dY))

"""
    rGCRFtoCIRS_iau2006([T::Type,] JD_TT::Number, dX::Number = 0, dY::Number = 0)

Compute the rotation that aligns the Geocentric Celestial Reference Frame (GCRF)
with the Celestial Intermediate Reference System (CIRS) at the Julian Day
`JD_TT` [TT] and considering the IERS EOP Data `dX` [rad] and `dY` [rad] \\(see
`get_iers_eop`). This algorithm uses the IAU-2006 theory.

The IERS EOP Data `dX` and `dY` accounts for the free-core nutation and time
dependent effects of the Celestial Intermediate Pole (CIP) position with respect
to the GCRF.

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the GCRF frame with the CIRS frame. The rotation
representation is selected by the optional parameter `T`.

"""
rGCRFtoCIRS_iau2006(JD_TT::Number, dX::Number = 0, dY::Number = 0) =
    rGCRFtoCIRS_iau2006(DCM, JD_TT, dX, dY)

function rGCRFtoCIRS_iau2006(::Type{DCM}, JD_TT::Number, dX::Number = 0,
                             dY::Number = 0)

    # Compute the CIP position w.r.t. GCRS.
    X, Y, s = precession_nutation_iau2006(JD_TT)

    # Add the corrections.
    X += dX
    Y += dY

    # Auxiliary variables.
    X² = X^2
    Y² = Y^2
    XY = X*Y

    # Compute the rotation matrix
    # ==========================================================================

    # This is the approximate value for:
    #
    #   a = 1/(1 + cos(d)), d = atan( sqrt( ( X^2 + Y^2 )/( 1 - X^2 - Y^2 ) ) )
    a = 1/2 + 1/8*(X² + Y²)

    D = DCM(1 - a*X²,    -a*XY, X,
               -a*XY, 1 - a*Y², Y,
                 -X ,      -Y , 1 - a*(X² + Y²))

    return create_rotation_matrix(-s, :Z)*D
end

rGCRFtoCIRS_iau2006(::Type{Quaternion}, JD_TT::Number, dX::Number = 0,
                    dY::Number = 0) =
    dcm_to_quat(rGCRFtoCIRS_iau2006(DCM, JD_TT, dX, dY))
