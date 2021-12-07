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
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Remarks
# ==============================================================================
#
# As mentioned in [1](p. 233), there is not an official definition of the TEME
# frame. Hence, in this package, it is considered the definition presented in
# [1](p. 233) in which the complete form of the Equation of Equinoxes is used.
# This seems to be the case when comparing the values shown in Table 3-6 [1](p.
# 232).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export r_teme_to_tod,  r_tod_to_teme
export r_teme_to_mod,  r_mod_to_teme
export r_teme_to_gcrf, r_gcrf_to_teme
export r_teme_to_pef,  r_pef_to_teme

#                                 TEME <=> TOD
# ==============================================================================

"""
    r_teme_to_tod([T,] jd_tt::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the True Equator Mean Equinox (TEME) frame with
the True of Date (TOD) frame at the Julian Day `jd_tt` [Terrestrial Time]. This
algorithm uses the IAU-76/FK5 theory and TEME definition in **[1]**(p. 233).
Notice that one can provide corrections for the nutation in obliquity
(`δΔϵ_1980`) [rad] and in longitude (`δΔψ_1980`) [rad] that are usually obtained
from IERS EOP Data (see [`get_iers_eop`](@ref)).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TEME frame with the TOD frame. The rotation
representation is selected by the optional parameter `T`.

# References

- **[1]**: Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
function r_teme_to_tod(jd_tt::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0)
    return r_teme_to_tod(DCM, jd_tt, δΔϵ_1980, δΔψ_1980)
end

function r_teme_to_tod(
    T::Type,
    jd_tt::Number,
    δΔϵ_1980::Number = 0,
    δΔψ_1980::Number = 0
)
    # Compute the nutation in the Julian Day (Terrestrial Time) `jd_tt`.
    mϵ_1980, Δϵ_1980, Δψ_1980 = nutation_fk5(jd_tt)

    # Add the corrections to the nutation in obliquity and longitude.
    Δϵ_1980 += δΔϵ_1980
    Δψ_1980 += δΔψ_1980

    # Evaluate the Delaunay parameters associated with the Moon in the interval
    # [0,2π]°.
    #
    # The parameters here were updated as stated in the errata [2].
    t_tt = (jd_tt - JD_J2000) / 36525
    r    = 360
    Ω_m  = @evalpoly(
        t_tt,
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
    r_tod_to_teme([T,] jd_tt::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the True of Date (TOD) frame with the True
Equator Mean Equinox (TEME) frame at the Julian Day `jd_tt` [Terrestrial Time].
This algorithm uses the IAU-76/FK5 theory and TEME definition in **[1]**(p.
233). Notice that one can provide corrections for the nutation in obliquity
(`δΔϵ_1980`) [rad] and in longitude (`δΔψ_1980`) [rad] that are usually obtained
from IERS EOP Data (see [`get_iers_eop`](@ref)).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TOD frame with the TEME frame. The rotation
representation is selected by the optional parameter `T`.

# References

- **[1]**: Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
function r_tod_to_teme(jd_tt::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0)
    return r_tod_to_teme(DCM, jd_tt, δΔϵ_1980, δΔψ_1980)
end

function r_tod_to_teme(
    T::Type,
    jd_tt::Number,
    δΔϵ_1980::Number = 0,
    δΔψ_1980::Number = 0
)
    return inv_rotation(r_teme_to_tod(T, jd_tt, δΔϵ_1980, δΔψ_1980))
end

#                                 TEME <=> MOD
# ==============================================================================

"""
    r_teme_to_mod([T,] jd_tt::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the True Equator Mean Equinox (TEME) frame with
the Mean of Date (MOD) frame at the Julian Day `jd_tt` [Terrestrial Time]. This
algorithm uses the IAU-76/FK5 theory and TEME definition in **[1]**(p. 233).
Notice that one can provide corrections for the nutation in obliquity
(`δΔϵ_1980`) [rad] and in longitude (`δΔψ_1980`) [rad] that are usually obtained
from IERS EOP Data (see [`get_iers_eop`](@ref)).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TEME frame with the MOD frame. The rotation
representation is selected by the optional parameter `T`.

# References

- **[1]**: Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
function r_teme_to_mod(jd_tt::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0)
    return r_teme_to_mod(DCM, jd_tt, δΔϵ_1980, δΔψ_1980)
end

function r_teme_to_mod(
    T::Type,
    jd_tt::Number,
    δΔϵ_1980::Number = 0,
    δΔψ_1980::Number = 0
)
    # Notice that, in this case, we will not use `r_teme_to_tod`, and `rTODtoMOD`
    # because this would call the function `nutation` twice, leading to a huge
    # performance drop. Hence, the code of those two functions is almost
    # entirely rewritten here.

    # Compute the nutation in the Julian Day (Terrestrial Time) `jd_tt`.
    mϵ_1980, Δϵ_1980, Δψ_1980 = nutation_fk5(jd_tt)

    # Add the corrections to the nutation in obliquity and longitude.
    Δϵ_1980 += δΔϵ_1980
    Δψ_1980 += δΔψ_1980

    # Compute the obliquity.
    ϵ_1980 = mϵ_1980 + Δϵ_1980

    # Evaluate the Delaunay parameters associated with the Moon in the interval
    # [0,2π]°.
    #
    # The parameters here were updated as stated in the errata [2].
    t_tt = (jd_tt - JD_J2000) / 36525
    r    = 360
    Ω_m  = @evalpoly(
        t_tt,
        125.04452222,
        -5r - 134.1362608,
        0.0020708,
        2.2e-6
    )
    Ω_m = mod(Ω_m, 360) * π / 180

    # Compute the equation of Equinoxes.
    #
    # According to [2], the constant unit before `sin(2Ω_m)` is also in [rad].
    Eq_equinox1982 = Δψ_1980 * cos(mϵ_1980) +
        (0.002640sin(1Ω_m) + 0.000063sin(2Ω_m)) * π / 648000

    # Compute the rotation TEME => TOD.
    r_tod_teme = angle_to_rot(T, -Eq_equinox1982, 0, 0, :ZYX)

    # Compute the rotation TOD => MOD.
    r_mod_tod = angle_to_rot(T, ϵ_1980, Δψ_1980, -mϵ_1980, :XZX)

    # Return the full rotation.
    return compose_rotation(r_tod_teme, r_mod_tod)
end

"""
    r_mod_to_teme([T,] jd_tt::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the Mean of Date (MOD) frame with the True
Equator Mean Equinox (TEME) frame at the Julian Day `jd_tt` [Terrestrial Time].
This algorithm uses the IAU-76/FK5 theory and TEME definition in **[1]**(p.
233). Notice that one can provide corrections for the nutation in obliquity
(`δΔϵ_1980`) [rad] and in longitude (`δΔψ_1980`) [rad] that are usually
obtained from IERS EOP Data (see [`get_iers_eop`](@ref)).  .

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the MOD frame with the TEME frame. The rotation
representation is selected by the optional parameter `T`.

# References

- **[1]**: Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
function r_mod_to_teme(jd_tt::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0)
    return r_mod_to_teme(DCM, jd_tt, δΔϵ_1980, δΔψ_1980)
end

function r_mod_to_teme(
    T::Type,
    jd_tt::Number,
    δΔϵ_1980::Number = 0,
    δΔψ_1980::Number = 0
)
    return inv_rotation(r_teme_to_mod(T, jd_tt, δΔϵ_1980, δΔψ_1980))
end

#                                TEME <=> GCRF
# ==============================================================================

"""
    r_teme_to_gcrf([T,] jd_tt::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the True Equator Mean Equinox (TEME) frame with
the Geocentric Celestial Reference Frame (GCRF) at the Julian Day `jd_tt`
[Terrestrial Time]. This algorithm uses the IAU-76/FK5 theory and TEME
definition in **[1]**(p. 233). Notice that one can provide corrections for the
nutation in obliquity (`δΔϵ_1980`) [rad] and in longitude (`δΔψ_1980`) [rad]
that are usually obtained from IERS EOP Data (see [`get_iers_eop`](@ref)).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

!!! info
    The EOP data related to the nutation of the obliquity (`δΔϵ_1980`) and the
    nutation of the longitude (`δΔψ_1980`) can be omitted. In this case, the
    GCRF frame is what is usually called J2000 reference frame.

# Returns

The rotation that aligns the TEME frame with the GCRF frame. The rotation
representation is selected by the optional parameter `T`.

# References

- **[1]**: Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
function r_teme_to_gcrf(jd_tt::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0)
    return r_teme_to_gcrf(DCM, jd_tt, δΔϵ_1980, δΔψ_1980)
end

function r_teme_to_gcrf(
    T::Type,
    jd_tt::Number,
    δΔϵ_1980::Number = 0,
    δΔψ_1980::Number = 0
)
    # Compute the rotation TEME => MOD.
    r_mod_teme  = r_teme_to_mod(T, jd_tt, δΔϵ_1980, δΔψ_1980)

    # Compute the rotation MOD => GCRF.
    r_gcrf_mod = r_mod_to_gcrf_fk5(T, jd_tt)

    # Return the full rotation.
    return compose_rotation(r_mod_teme, r_gcrf_mod)
end

"""
    r_gcrf_to_teme([T,] jd_tt::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the GCRF frame with the True Equator Mean
Equinox (TEME) frame at the Julian Day `jd_tt` [Terrestrial Time]. This
algorithm uses the IAU-76/FK5 theory and TEME definition in **[1]**(p. 233).
Notice that one can provide corrections for the nutation in obliquity
(`δΔϵ_1980`) [rad] and in longitude (`δΔψ_1980`) [rad] that are usually obtained
from IERS EOP Data (see [`get_iers_eop`](@ref)).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

!!! info
    The EOP data related to the nutation of the obliquity (`δΔϵ_1980`) and the
    nutation of the longitude (`δΔψ_1980`) can be omitted. In this case, the
    GCRF frame is what is usually called J2000 reference frame.

# Returns

The rotation that aligns the GCRF frame with the TEME frame. The rotation
representation is selected by the optional parameter `T`.

# References

- **[1]**: Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
function r_gcrf_to_teme(jd_tt::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0)
    return r_gcrf_to_teme(DCM, jd_tt, δΔϵ_1980, δΔψ_1980)
end

function r_gcrf_to_teme(
    T::Type,
    jd_tt::Number,
    δΔϵ_1980::Number = 0,
    δΔψ_1980::Number = 0
)
    return inv_rotation(r_teme_to_gcrf(T, jd_tt, δΔϵ_1980, δΔψ_1980))
end

#                                 TEME <=> PEF
# ==============================================================================

"""
    r_teme_to_pef([T,] jd_tt::Number)

Compute the rotation that aligns the True Equator Mean Equinox (TEME) frame with
the Pseudo-Earth Fixed (PEF) frame at the Julian Day `jd_tt` [Terrestrial Time].
This algorithm uses the IAU-76/FK5 theory and TEME definition in **[1]**(p.
233).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the TEME frame with the PEF frame. The rotation
representation is selected by the optional parameter `T`.

# References

- **[1]**: Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
r_teme_to_pef(jd_ut1::Number) = r_teme_to_pef(DCM, jd_ut1)

function r_teme_to_pef(T::Type, jd_ut1::Number)
    # Compute the Greenwich Mean Sidereal Time.
    θ_gmst = jd_to_gmst(jd_ut1)

    # Compute the rotation.
    return angle_to_rot(T, θ_gmst, 0, 0, :ZYX)
end

"""
    r_pef_to_teme([T,] jd_tt::Number)

Compute the rotation that aligns the Pseudo-Earth Fixed (PEF) frame with the
True Equator Mean Equinox (TEME) frame at the Julian Day `jd_tt` [Terrestrial
Time]. This algorithm uses the IAU-76/FK5 theory and TEME definition in
**[1]**(p. 233).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Returns

The rotation that aligns the PEF frame with the TEME frame. The rotation
representation is selected by the optional parameter `T`.

# References

- **[1]**: Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
r_pef_to_teme(jd_ut1::Number) = r_pef_to_teme(DCM, jd_ut1)
r_pef_to_teme(T::Type, jd_ut1::Number) = inv_rotation(r_teme_to_pef(T, jd_ut1))
