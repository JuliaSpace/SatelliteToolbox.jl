# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions related with True Equator Mean Equinox (TEME) reference frame.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Remarks
#
# As mentioned in [1, p. 233], there is not an official definition of the TEME
# frame. Hence, in this package, it is considered the definition presented in
# [1, p. 233] in which the complete form of the Equation of Equinoxes is used.
# This seems to be the case when comparing the values shown in Table 3-6 [1, p.
# 232].
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export r_teme_to_tod,  r_tod_to_teme
export r_teme_to_mod,  r_mod_to_teme
export r_teme_to_gcrf, r_gcrf_to_teme
export r_teme_to_pef,  r_pef_to_teme

#                                 TEME <=> TOD
# ==============================================================================

"""
    r_teme_to_tod([T,] JD_TT::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the True Equator Mean Equinox (TEME) frame with
the True of Date (TOD) frame at the Julian Day `JD_TT` [Terrestrial Time]. This
algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p. 233]. Notice
that one can provide corrections for the nutation in obliquity (`δΔϵ_1980`)
\\[rad] and in longitude (`δΔψ_1980`) \\[rad] that are usually obtained from
IERS EOP Data (see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TEME frame with the TOD frame. The rotation
representation is selected by the optional parameter `T`.

"""
function r_teme_to_tod(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0)
    return r_teme_to_tod(DCM, JD_TT, δΔϵ_1980, δΔψ_1980)
end

function r_teme_to_tod(
    T::Type,
    JD_TT::Number,
    δΔϵ_1980::Number = 0,
    δΔψ_1980::Number = 0
)
    # Compute the nutation in the Julian Day (Terrestrial Time) `JD_TT`.
    mϵ_1980, Δϵ_1980, Δψ_1980 = nutation_fk5(JD_TT)

    # Add the corrections to the nutation in obliquity and longitude.
    Δϵ_1980 += δΔϵ_1980
    Δψ_1980 += δΔψ_1980

    # Compute the obliquity.
    ϵ_1980 = mϵ_1980 + Δϵ_1980

    # Evaluate the Delaunay parameters associated with the Moon in the interval
    # [0,2π]°.
    #
    # The parameters here were updated as stated in the errata [2].
    T_TT = (JD_TT - JD_J2000)/36525
    r    = 360
    Ω_m  = @evalpoly(
        T_TT,
        125.04452222,
        -5r - 134.1362608,
        0.0020708,
        2.2e-6
    )
    Ω_m = mod(Ω_m, 360) * π / 180

    # Compute the equation of Equinoxes.
    #
    # According to [2], the constant unit before `sin(2Ω_m)` is also in [rad].
    Eq_equinox1982 = Δψ_1980*cos(mϵ_1980) +
        (0.002640sin(1Ω_m) + 0.000063sin(2Ω_m)) * π / 648000

    # Compute the rotation.
    return angle_to_rot(T, -Eq_equinox1982, 0, 0, :ZYX)
end

"""
    r_tod_to_teme([T,] JD_TT::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the True of Date (TOD) frame with the True
Equator Mean Equinox (TEME) frame at the Julian Day `JD_TT` [Terrestrial Time].
This algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p. 233].
Notice that one can provide corrections for the nutation in obliquity
(`δΔϵ_1980`) \\[rad] and in longitude (`δΔψ_1980`) \\[rad] that are usually
obtained from IERS EOP Data (see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TOD frame with the TEME frame. The rotation
representation is selected by the optional parameter `T`.

"""
function r_tod_to_teme(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0)
    return r_tod_to_teme(DCM, JD_TT, δΔϵ_1980, δΔψ_1980)
end

function r_tod_to_teme(
    T::Type,
    JD_TT::Number,
    δΔϵ_1980::Number = 0,
    δΔψ_1980::Number = 0
)
    return inv_rotation(r_teme_to_tod(T, JD_TT, δΔϵ_1980, δΔψ_1980))
end

#                                 TEME <=> MOD
# ==============================================================================

"""
    r_teme_to_mod([T,] JD_TT::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the True Equator Mean Equinox (TEME) frame with
the Mean of Date (MOD) frame at the Julian Day `JD_TT` [Terrestrial Time]. This
algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p. 233]. Notice
that one can provide corrections for the nutation in obliquity (`δΔϵ_1980`)
\\[rad] and in longitude (`δΔψ_1980`) \\[rad] that are usually obtained from
IERS EOP Data (see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TEME frame with the MOD frame. The rotation
representation is selected by the optional parameter `T`.

"""
function r_teme_to_mod(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0)
    return r_teme_to_mod(DCM, JD_TT, δΔϵ_1980, δΔψ_1980)
end

function r_teme_to_mod(
    T::Type,
    JD_TT::Number,
    δΔϵ_1980::Number = 0,
    δΔψ_1980::Number = 0
)
    # Notice that, in this case, we will not use `r_teme_to_tod`, and `rTODtoMOD`
    # because this would call the function `nutation` twice, leading to a huge
    # performance drop. Hence, the code of those two functions is almost
    # entirely rewritten here.

    # Compute the nutation in the Julian Day (Terrestrial Time) `JD_TT`.
    mϵ_1980, Δϵ_1980, Δψ_1980 = nutation_fk5(JD_TT)

    # Add the corrections to the nutation in obliquity and longitude.
    Δϵ_1980 += δΔϵ_1980
    Δψ_1980 += δΔψ_1980

    # Compute the obliquity.
    ϵ_1980 = mϵ_1980 + Δϵ_1980

    # Evaluate the Delaunay parameters associated with the Moon in the interval
    # [0,2π]°.
    #
    # The parameters here were updated as stated in the errata [2].
    T_TT = (JD_TT - JD_J2000) / 36525
    r    = 360
    Ω_m  = @evalpoly(
        T_TT,
        125.04452222,
        -5r - 134.1362608,
        0.0020708,
        2.2e-6
    )
    Ω_m = mod(Ω_m, 360) * π / 180

    # Compute the equation of Equinoxes.
    #
    # According to [2], the constant unit before `sin(2Ω_m)` is also in [rad].
    Eq_equinox1982 = Δψ_1980*cos(mϵ_1980) +
        (0.002640sin(1Ω_m) + 0.000063sin(2Ω_m)) * π / 648000

    # Compute the rotation TEME => TOD.
    r_TOD_TEME = angle_to_rot(T, -Eq_equinox1982, 0, 0, :ZYX)

    # Compute the rotation TOD => MOD.
    r_MOD_TOD = angle_to_rot(T, ϵ_1980, Δψ_1980, -mϵ_1980, :XZX)

    # Return the full rotation.
    return compose_rotation(r_TOD_TEME, r_MOD_TOD)
end

"""
    r_mod_to_teme([T,] JD_TT::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the Mean of Date (MOD) frame with the True
Equator Mean Equinox (TEME) frame at the Julian Day `JD_TT` [Terrestrial Time].
This algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p. 233].
Notice that one can provide corrections for the nutation in obliquity
(`δΔϵ_1980`) \\[rad] and in longitude (`δΔψ_1980`) \\[rad] that are usually
obtained from IERS EOP Data (see `get_iers_eop`).  .

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the MOD frame with the TEME frame. The rotation
representation is selected by the optional parameter `T`.

"""
function r_mod_to_teme(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0)
    return r_mod_to_teme(DCM, JD_TT, δΔϵ_1980, δΔψ_1980)
end

function r_mod_to_teme(
    T::Type,
    JD_TT::Number,
    δΔϵ_1980::Number = 0,
    δΔψ_1980::Number = 0
)
    return inv_rotation(r_teme_to_mod(T, JD_TT, δΔϵ_1980, δΔψ_1980))
end

#                                TEME <=> GCRF
# ==============================================================================

"""
    r_teme_to_gcrf([T,] JD_TT::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the True Equator Mean Equinox (TEME) frame with
the Geocentric Celestial Reference Frame (GCRF) at the Julian Day `JD_TT`
[Terrestrial Time]. This algorithm uses the IAU-76/FK5 theory and TEME
definition in [1, p. 233]. Notice that one can provide corrections for the
nutation in obliquity (`δΔϵ_1980`) \\[rad] and in longitude (`δΔψ_1980`) \\[rad]
that are usually obtained from IERS EOP Data (see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TEME frame with the GCRF frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The EOP data related to the nutation of the obliquity (`δΔϵ_1980`) and the
nutation of the longitude (`δΔψ_1980`) can be omitted. In this case, the GCRF
frame is what is usually called J2000 reference frame.

"""
function r_teme_to_gcrf(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0)
    return r_teme_to_gcrf(DCM, JD_TT, δΔϵ_1980, δΔψ_1980)
end

function r_teme_to_gcrf(
    T::Type,
    JD_TT::Number,
    δΔϵ_1980::Number = 0,
    δΔψ_1980::Number = 0
)
    # Compute the rotation TEME => MOD.
    r_MOD_TEME  = r_teme_to_mod(T, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Compute the rotation MOD => GCRF.
    r_GCRF_MOD = r_mod_to_gcrf_fk5(T, JD_TT)

    # Return the full rotation.
    return compose_rotation(r_MOD_TEME, r_GCRF_MOD)
end

"""
    r_gcrf_to_teme([T,] JD_TT::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the GCRF frame with the True Equator Mean
Equinox (TEME) frame at the Julian Day `JD_TT` [Terrestrial Time]. This
algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p. 233]. Notice
that one can provide corrections for the nutation in obliquity (`δΔϵ_1980`)
\\[rad] and in longitude (`δΔψ_1980`) \\[rad] that are usually obtained from
IERS EOP Data (see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the GCRF frame with the TEME frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The EOP data related to the nutation of the obliquity (`δΔϵ_1980`) and the
nutation of the longitude (`δΔψ_1980`) can be omitted. In this case, the GCRF
frame is what is usually called J2000 reference frame.

"""
function r_gcrf_to_teme(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0)
    return r_gcrf_to_teme(DCM, JD_TT, δΔϵ_1980, δΔψ_1980)
end

function r_gcrf_to_teme(
    T::Type,
    JD_TT::Number,
    δΔϵ_1980::Number = 0,
    δΔψ_1980::Number = 0
)
    return inv_rotation(r_teme_to_gcrf(T, JD_TT, δΔϵ_1980, δΔψ_1980))
end

#                                 TEME <=> PEF
# ==============================================================================

"""
    r_teme_to_pef([T,] JD_TT::Number)

Compute the rotation that aligns the True Equator Mean Equinox (TEME) frame with
the Pseudo-Earth Fixed (PEF) frame at the Julian Day `JD_TT` [Terrestrial Time].
This algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p. 233].

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TEME frame with the PEF frame. The rotation
representation is selected by the optional parameter `T`.

"""
r_teme_to_pef(JD_UT1::Number) = r_teme_to_pef(DCM, JD_UT1)

function r_teme_to_pef(T::Type, JD_UT1::Number)
    # Compute the Greenwich Mean Sidereal Time.
    θ_gmst = jd_to_gmst(JD_UT1)

    # Compute the rotation.
    return angle_to_rot(T, θ_gmst, 0, 0, :ZYX)
end

"""
    r_pef_to_teme([T,] JD_TT::Number)

Compute the rotation that aligns the Pseudo-Earth Fixed (PEF) frame with the
True Equator Mean Equinox (TEME) frame at the Julian Day `JD_TT` [Terrestrial
Time]. This algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p.
233].

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the PEF frame with the TEME frame. The rotation
representation is selected by the optional parameter `T`.

"""
r_pef_to_teme(JD_UT1::Number) = r_pef_to_teme(DCM, JD_UT1)
r_pef_to_teme(T::Type, JD_UT1::Number) = inv_rotation(r_teme_to_pef(T, JD_UT1))
