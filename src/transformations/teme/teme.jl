#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
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
#
# Changelog
#
# 2018-05-24: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export rTEMEtoTOD,  rTODtoTEME
export rTEMEtoMOD,  rMODtoTEME
export rTEMEtoGCRF, rGCRFtoTEME
export rTEMEtoPEF,  rPEFtoTEME

#                                 TEME <=> TOD
# ==============================================================================

"""
    function rTEMEtoTOD([T,] JD_TT::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the True Equator Mean Equinox (TEME) frame with
the True of Date (TOD) frame at the Julian Day (Terrestrial Time) `JD_TT`. This
algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p. 233]. Notice
that one can provide corrections for the nutation in obliquity (`δΔϵ`) and in
longitude (`δΔψ`) that are usually obtained from IERS EOP Data (see
`get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Args

* `T`: (OPTIONAL) Type of the rotation representation (**Default**: `DCM`).
* `JD_TT`: Julian Day [Terrestrial Time].
* `δΔϵ_1980`: (OPTIONAL) Correction in the nutation of the obliquity [rad]
              (**Default** = 0).
* `δΔψ_1980`: (OPTIONAL) Correction in the nutation of the longitude [rad]
              (**Default** = 0).

# Returns

The rotation that aligns the TEME frame with the TOD frame. The rotation
representation is selected by the optional parameter `T`.

"""
rTEMEtoTOD(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
    rTEMEtoTOD(DCM, JD_TT, δΔϵ_1980, δΔψ_1980)

function rTEMEtoTOD(T::Type,
                    JD_TT::Number,
                    δΔϵ_1980::Number = 0,
                    δΔψ_1980::Number = 0)

    # Compute the nutation in the Julian Day (Terrestrial Time) `JD_TT`.
    (mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(JD_TT)

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
    Ω_m  = 125.04452222 - (5r + 134.1362608)*T_TT +
                          0.0020708*T_TT^2 +
                          2.2e-6*T_TT^3
    Ω_m  = mod(Ω_m, 360)*pi/180

    # Compute the equation of Equinoxes.
    #
    # According to [2], the constant unit before `sin(2Ω_m)` is also in [rad].
    Eq_equinox1982 = Δψ_1980*cos(mϵ_1980) +
        ( 0.002640*sin(1Ω_m) + 0.000063*sin(2Ω_m) )*pi/648000

    # Compute the rotation.
    angle2rot(T, -Eq_equinox1982, 0, 0, :ZYX)
end

"""
    function rTODtoTEME([T,] JD_TT::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the True of Date (TOD) frame with the True
Equator Mean Equinox (TEME) frame at the Julian Day (Terrestrial Time) `JD_TT`.
This algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p. 233].
Notice that one can provide corrections for the nutation in obliquity (`δΔϵ`)
and in longitude (`δΔψ`) that are usually obtained from IERS EOP Data (see
`get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Args

* `T`: (OPTIONAL) Type of the rotation representation (**Default**: `DCM`).
* `JD_TT`: Julian Day [Terrestrial Time].
* `δΔϵ_1980`: (OPTIONAL) Correction in the nutation of the obliquity [rad]
              (**Default** = 0).
* `δΔψ_1980`: (OPTIONAL) Correction in the nutation of the longitude [rad]
              (**Default** = 0).

# Returns

The rotation that aligns the TOD frame with the TEME frame. The rotation
representation is selected by the optional parameter `T`.

"""
rTODtoTEME(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
    rTODtoTEME(DCM, JD_TT, δΔϵ_1980, δΔψ_1980)

rTODtoTEME(T::Type, JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
    inv_rotation(rTEMEtoTOD(T, JD_TT, δΔϵ_1980, δΔψ_1980))

#                                 TEME <=> MOD
# ==============================================================================

"""
    function rTEMEtoMOD([T,] JD_TT::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the True Equator Mean Equinox (TEME) frame with
the Mean of Date (MOD) frame at the Julian Day (Terrestrial Time) `JD_TT`. This
algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p. 233]. Notice
that one can provide corrections for the nutation in obliquity (`δΔϵ`) and in
longitude (`δΔψ`) that are usually obtained from IERS EOP Data (see
`get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Args

* `T`: (OPTIONAL) Type of the rotation representation (**Default**: `DCM`).
* `JD_TT`: Julian Day [Terrestrial Time].
* `δΔϵ_1980`: (OPTIONAL) Correction in the nutation of the obliquity [rad]
              (**Default** = 0).
* `δΔψ_1980`: (OPTIONAL) Correction in the nutation of the longitude [rad]
              (**Default** = 0).

# Returns

The rotation that aligns the TEME frame with the MOD frame. The rotation
representation is selected by the optional parameter `T`.

"""
rTEMEtoMOD(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
    rTEMEtoMOD(DCM, JD_TT, δΔϵ_1980, δΔψ_1980)

function rTEMEtoMOD(T::Type,
                    JD_TT::Number,
                    δΔϵ_1980::Number = 0,
                    δΔψ_1980::Number = 0)

    # Notice that, in this case, we will not use `rTEMEtoTOD`, and `rTODtoMOD`
    # because this would call the function `nutation` twice, leading to a huge
    # performance drop. Hence, the code of those two functions is almost
    # entirely rewritten here.

    # Compute the nutation in the Julian Day (Terrestrial Time) `JD_TT`.
    (mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(JD_TT)

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
    Ω_m  = 125.04452222 - (5r + 134.1362608)*T_TT +
                          0.0020708*T_TT^2 +
                          2.2e-6*T_TT^3
    Ω_m  = mod(Ω_m, 360)*pi/180

    # Compute the equation of Equinoxes.
    #
    # According to [2], the constant unit before `sin(2Ω_m)` is also in [rad].
    Eq_equinox1982 = Δψ_1980*cos(mϵ_1980) +
        ( 0.002640*sin(1Ω_m) + 0.000063*sin(2Ω_m) )*pi/648000

    # Compute the rotation TEME => TOD.
    r_TOD_TEME = angle2rot(T, -Eq_equinox1982, 0, 0, :ZYX)

    # Compute the rotation TOD => MOD.
    r_MOD_TOD = angle2rot(T, ϵ_1980, Δψ_1980, -mϵ_1980, :XZX)

    # Return the full rotation.
    compose_rotation(r_TOD_TEME, r_MOD_TOD)
end

"""
    function rMODtoTEME([T,] JD_TT::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the Mean of Date (MOD) frame with the True
Equator Mean Equinox (TEME) frame at the Julian Day (Terrestrial Time) `JD_TT`.
This algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p. 233].
Notice that one can provide corrections for the nutation in obliquity (`δΔϵ`)
and in longitude (`δΔψ`) that are usually obtained from IERS EOP Data (see
`get_iers_eop`).  .

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Args

* `T`: (OPTIONAL) Type of the rotation representation (**Default**: `DCM`).
* `JD_TT`: Julian Day [Terrestrial Time].
* `δΔϵ_1980`: (OPTIONAL) Correction in the nutation of the obliquity [rad]
              (**Default** = 0).
* `δΔψ_1980`: (OPTIONAL) Correction in the nutation of the longitude [rad]
              (**Default** = 0).

# Returns

The rotation that aligns the MOD frame with the TEME frame. The rotation
representation is selected by the optional parameter `T`.

"""
rMODtoTEME(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
    rMODtoTEME(DCM, JD_TT, δΔϵ_1980, δΔψ_1980)

rMODtoTEME(T::Type, JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
    inv_rotation(rTEMEtoMOD(T, JD_TT, δΔϵ_1980, δΔψ_1980))

#                                TEME <=> GCRF
# ==============================================================================

"""
    function rTEMEtoGCRF([T,] JD_TT::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the True Equator Mean Equinox (TEME) frame with
the Geocentric Celestial Reference Frame (GCRF) at the Julian Day (Terrestrial
Time) `JD_TT`. This algorithm uses the IAU-76/FK5 theory and TEME definition in
[1, p. 233]. Notice that one can provide corrections for the nutation in
obliquity (`δΔϵ`) and in longitude (`δΔψ`) that are usually obtained from IERS
EOP Data (see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Args

* `T`: (OPTIONAL) Type of the rotation representation (**Default**: `DCM`).
* `JD_TT`: Julian Day [Terrestrial Time].
* `δΔϵ_1980`: (OPTIONAL) Correction in the nutation of the obliquity [rad]
              (**Default** = 0).
* `δΔψ_1980`: (OPTIONAL) Correction in the nutation of the longitude [rad]
              (**Default** = 0).

# Returns

The rotation that aligns the TEME frame with the GCRF frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The EOP data related to the nutation of the obliquity (`δΔϵ_1980`) and the
nutation of the longitude (`δΔψ_1980`) can be omitted. In this case, the GCRF
frame is what is usually called J2000 reference frame.

"""
rTEMEtoGCRF(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
    rTEMEtoGCRF(DCM, JD_TT, δΔϵ_1980, δΔψ_1980)

function rTEMEtoGCRF(T::Type,
                     JD_TT::Number,
                     δΔϵ_1980::Number = 0,
                     δΔψ_1980::Number = 0)

    # Compute the rotation TEME => MOD.
    r_MOD_TEME  = rTEMEtoMOD(T, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Compute the rotation MOD => GCRF.
    r_GCRF_MOD = rMODtoGCRF_fk5(T, JD_TT)

    # Return the full rotation.
    compose_rotation(r_MOD_TEME, r_GCRF_MOD)
end

"""
    function rGCRFtoTEME([T,] JD_TT::Number [, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0])

Compute the rotation that aligns the GCRF frame with the True Equator Mean
Equinox (TEME) frame at the Julian Day (Terrestrial Time) `JD_TT`. This
algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p. 233]. Notice
that one can provide corrections for the nutation in obliquity (`δΔϵ`) and in
longitude (`δΔψ`) that are usually obtained from IERS EOP Data (see
`get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Args

* `T`: (OPTIONAL) Type of the rotation representation (**Default**: `DCM`).
* `JD_TT`: Julian Day [Terrestrial Time].

# Returns

The rotation that aligns the GCRF frame with the TEME frame. The rotation
representation is selected by the optional parameter `T`.

# Remarks

The EOP data related to the nutation of the obliquity (`δΔϵ_1980`) and the
nutation of the longitude (`δΔψ_1980`) can be omitted. In this case, the GCRF
frame is what is usually called J2000 reference frame.

"""
rGCRFtoTEME(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
    rGCRFtoTEME(DCM, JD_TT, δΔϵ_1980, δΔψ_1980)

rGCRFtoTEME(T::Type, JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
   inv_rotation(rTEMEtoGCRF(T, JD_TT, δΔϵ_1980, δΔψ_1980))

#                                 TEME <=> PEF
# ==============================================================================

"""
    function rTEMEtoPEF([T,] JD_TT::Number)

Compute the rotation that aligns the True Equator Mean Equinox (TEME) frame with
the Pseudo-Earth Fixed (PEF) frame at the Julian Day (Terrestrial Time) `JD_TT`.
This algorithm uses the IAU-76/FK5 theory and TEME definition in [1, p. 233].

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Args

* `T`: (OPTIONAL) Type of the rotation representation (**Default**: `DCM`).
* `JD_TT`: Julian Day [Terrestrial Time].

# Returns

The rotation that aligns the TEME frame with the PEF frame. The rotation
representation is selected by the optional parameter `T`.

"""
rTEMEtoPEF(JD_UT1::Number) = rTEMEtoPEF(DCM, JD_UT1)

function rTEMEtoPEF(T::Type, JD_UT1::Number)
    # Compute the Greenwich Mean Sidereal Time.
    θ_gmst = JDtoGMST(JD_UT1)

    # Compute the rotation.
    angle2rot(T, θ_gmst, 0, 0, :ZYX)
end

"""
    function rPEFtoTEME([T,] JD_TT::Number)

Compute the rotation that aligns the Pseudo-Earth Fixed (PEF) frame with the
True Equator Mean Equinox (TEME) frame at the Julian Day (Terrestrial Time)
`JD_TT`.  This algorithm uses the IAU-76/FK5 theory and TEME definition in [1,
p. 233].

The rotation type is described by the optional variable `T`. If it is `DCM`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`DCM`.

# Args

* `T`: (OPTIONAL) Type of the rotation representation (**Default**: `DCM`).
* `JD_TT`: Julian Day [Terrestrial Time].

# Returns

The rotation that aligns the PEF frame with the TEME frame. The rotation
representation is selected by the optional parameter `T`.

"""
rPEFtoTEME(JD_UT1::Number)          = rPEFtoTEME(DCM, JD_UT1)
rPEFtoTEME(T::Type, JD_UT1::Number) = inv_rotation(rTEMEtoPEF(T, JD_UT1))
