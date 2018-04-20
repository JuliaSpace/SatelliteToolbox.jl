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
#   Functions related with the model IAU-76/FK5.
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
# Changelog
#
# 2018-04-09: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

using Rotations

export rTODtoMOD_fk5,  rMODtoTOD_fk5
export rMODtoGCRF_fk5, rGCRFtoMOD_fk5

export rTODtoGCRF_fk5, rGCRFtoTOD_fk5

################################################################################
#                            IAU-76/FK5 Reductions
################################################################################
#
# The conversion between the Geocentric Celestial Reference Frame (GCRF) to the
# International Terrestrial Reference Frame (ITRF) is done by means of:
#
#                    GCRF <=> MOD <=> TOD <=> PEF <=> ITRF
#
# in which:
#   - MOD: Mean of Date frame.
#   - TOD: True of Date frame.
#   - PEF: Pseudo-Earth fixed frame.
#
# Every rotation will be coded as a function using the IAU-76/FK5 theory.
# Additionally, composed rotations will also available. In general, the API is:
#
#   function r<Origin Frame>to<Destination Frame>_fk5
#
# The arguments vary depending on the origin and destination frame and should be
# verified using the function documentation.
#
################################################################################

################################################################################
#                               Single Rotations
################################################################################

#                                 TOD <=> MOD
# ==============================================================================

"""
### function rTODtoMOD_fk5([T,] JD_TT::Number [, δΔϵ_1980::Number, δΔψ_1980::Number])

Compute the rotation that aligns the True of Date (TOD) frame with the Mean of
Date (MOD) frame at the Julian Day (Terrestrial Time) `JD_TT`. This algorithm
uses the IAU-76/FK5 theory. Notice that one can provide corrections for the
nutation in obliquity (`δΔϵ`) and in longitude (`δΔψ`) that are usually obtained
from IERS EOP Data (see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT** = `Matrix`).
* JD_TT: Julian Day [Terrestrial Time].
* δΔϵ_1980: (OPTIONAL) Correction in the nutation of the obliquity [rad]
            (**DEFAULT** = 0).
* δΔψ_1980: (OPTIONAL) Correction in the nutation of the longitude [rad]
            (**DEFAULT** = 0).

##### Returns

The rotation that aligns the TOD frame with the MOD frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

The True of Date (TOD) frame is rotated into the Mean of Date (MOD) frame
considering the 1980 IAU Theory of Nutation. The IERS EOP corrections must be
added if one wants to make the rotation consistent with the Geocentric Celestial
Reference Systems (GCRS).

"""
rTODtoMOD_fk5(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
    rTODtoMOD_fk5(Matrix, JD_TT, δΔϵ_1980, δΔψ_1980)

function rTODtoMOD_fk5(::Type{Matrix},
                       JD_TT::Number,
                       δΔϵ_1980::Number = 0,
                       δΔψ_1980::Number = 0)

    # Compute the nutation in the Juliay Day (Terrestrial Time) `JD_TT`.
    (mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(JD_TT)

    # Add the corrections to the nutation in obliquity and longitude.
    Δϵ_1980 += δΔϵ_1980
    Δψ_1980 += δΔψ_1980

    # Compute the obliquity.
    ϵ_1980 = mϵ_1980 + Δϵ_1980

    # Compute and return the Direction Cosine Matrix.
    angle2dcm(ϵ_1980, Δψ_1980, -mϵ_1980, "XZX")
end

function rTODtoMOD_fk5(::Type{Quaternion},
                       JD_TT::Number,
                       δΔϵ_1980::Number = 0,
                       δΔψ_1980::Number = 0)
    # Compute the nutation in the Juliay Day (Terrestrial Time) `JD_TT`.
    (mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(JD_TT)

    # Add the corrections to the nutation in obliquity and longitude.
    Δϵ_1980 += δΔϵ_1980
    Δψ_1980 += δΔψ_1980

    # Compute the obliquity.
    ϵ_1980 = mϵ_1980 + Δϵ_1980

    # Compute the three rotations using Quaternions.
    q1_x = Quaternion([ cos(ϵ_1980/2);   sin(ϵ_1980/2); 0;        0      ])
    q2_z = Quaternion([cos(Δψ_1980/2);        0;        0; sin(Δψ_1980/2)])
    q3_x = Quaternion([cos(mϵ_1980/2); -sin(mϵ_1980/2); 0;        0      ])

    # Return the complete quaternion.
    q1_x*q2_z*q3_x
end

"""
### function rMODtoTOD_fk5([T,] JD_TT::Number [, δΔϵ_1980::Number, δΔψ_1980::Number])

Compute the rotation that aligns the Mean of Date (MOD) frame with the True of
Date (TOD) frame at the Julian Day (Terrestrial Time) `JD_TT`. This algorithm
uses the IAU-76/FK5 theory. Notice that one can provide corrections for the
nutation in obliquity (`δΔϵ`) and in longitude (`δΔψ`) that are usually obtained
from IERS EOP Data (see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT** = `Matrix`).
* JD_TT: Julian Day [Terrestrial Time].
* δΔϵ_1980: (OPTIONAL) Correction in the nutation of the obliquity [rad]
            (**DEFAULT** = 0).
* δΔψ_1980: (OPTIONAL) Correction in the nutation of the longitude [rad]
            (**DEFAULT** = 0).

##### Returns

The rotation that aligns the MOD frame with the TOD frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

The Mean of Date (MOD) frame is rotated into the True of Date (TOD) frame
considering the 1980 IAU Theory of Nutation. The IERS EOP corrections must be
added if one wants to make the rotation consistent with the Geocentric Celestial
Reference Systems (GCRS).

"""
rMODtoTOD_fk5(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
    rTODtoMOD_fk5(Matrix, JD_TT, δΔϵ_1980, δΔψ_1980)'

function rMODtoTOD_fk5(::Type{Matrix},
                       JD_TT::Number,
                       δΔϵ_1980::Number = 0,
                       δΔψ_1980::Number = 0)
    rTODtoMOD_fk5(Matrix, JD_TT, δΔϵ_1980, δΔψ_1980)'
end

function rMODtoTOD_fk5(::Type{Quaternion},
                       JD_TT::Number,
                       δΔϵ_1980::Number = 0,
                       δΔψ_1980::Number = 0)
    conj(rTODtoMOD_fk5(Quaternion, JD_TT, δΔϵ_1980, δΔψ_1980))
end

#                                 MOD <=> GCRF
# ==============================================================================

"""
### function rMODtoGCRF_fk5([T,] JD_TT::Number)

Compute the rotation that aligns the Mean of Date (MOD) frame with the
Geocentric Celestial Reference Frame (GCRF) at the Julian Day (Terrestrial Time)
`JD_TT`. This algorithm uses the IAU-76/FK5 theory.

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT**: `Matrix`).
* JD_TT: Julian Day [Terrestrial Time].

##### Returns

The rotation that aligns the MOD frame with the GCRF frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

The Mean of Date (MOD) frame is rotated into the Geocentric Celestial Reference
Frame (GCRF) considering the IAU 1976 Precession model.

Notice that if the conversion `TOD => MOD` is performed **without** considering
the EOP corrections, then the GCRF obtained by this rotation is what is usually
called the J2000 reference frame.

"""
rMODtoGCRF_fk5(JD_TT::Number) = rMODtoGCRF_fk5(Matrix,JD_TT)

function rMODtoGCRF_fk5(::Type{Matrix},JD_TT::Number)
    (ζ,Θ,z) = precession_fk5(JD_TT)
    angle2dcm(z,-Θ,ζ,"ZYZ")
end

function rMODtoGCRF_fk5(::Type{Quaternion},JD_TT::Number)
    (ζ,Θ,z) = precession_fk5(JD_TT)

    # Compute the three rotations using Quaternions.
    q1_z = Quaternion([cos(z/2); 0;     0    ; sin(z/2)])
    q2_y = Quaternion([cos(Θ/2); 0; -sin(Θ/2);    0    ])
    q3_z = Quaternion([cos(ζ/2); 0;     0    ; sin(ζ/2)])

    # Return the complete quaternion.
    q1_z*q2_y*q3_z
end

"""
### function rGCRFtoMOD_fk5([T,] JD_TT::Number)

Compute the rotation that aligns the Geocentric Celestial Reference Frame (GCRF)
with the Mean of Date (MOD) frame at the Julian Day (Terrestrial Time) `JD_TT`.
This algorithm uses the IAU-76/FK5 theory.

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT**: `Matrix`).
* JD_TT: Julian Day [Terrestrial Time].

##### Returns

The rotation that aligns the GCRF frame with the MOD frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

The Geocentric Celestial Reference Frame (GCRF) is rotated into the Mean of Date
(MOD) frame considering the IAU 1976 Precession model.

Notice that if the conversion `MOD => TOD` is performed **without** considering
the EOP corrections, then the GCRF in this rotation is what is usually called
the J2000 reference frame.

"""

rGCRFtoMOD_fk5(JD_TT::Number) = rMODtoGCRF_fk5(Matrix,JD_TT)'

function rGCRFtoMOD_fk5(::Type{Matrix},JD_TT::Number)
    rMODtoGCRF_fk5(Matrix, JD_TT)'
end

function rGCRFtoMOD_fk5(::Type{Quaternion},JD_TT::Number)
    conj(rMODtoGCRF_fk5(Quaternion, JD_TT))
end

################################################################################
#                              Multiple Rotations
################################################################################

#                                 TOD <=> GCRF
# ==============================================================================

"""
### function rTODtoGCRF_fk5([T,] JD_TT::Number [, δΔϵ_1980::Number, δΔψ_1980::Number])

Compute the rotation that aligns the True of Date (TOD) frame with the
Geocentric Celestial Reference Frame (GCRF) at the Julian Day (Terrestrial Time)
`JD_TT`. This algorithm uses the IAU-76/FK5 theory. Notice that one can provide
corrections for the nutation in obliquity (`δΔϵ`) and in longitude (`δΔψ`) that
are usually obtained from IERS EOP Data (see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT** = `Matrix`).
* JD_TT: Julian Day [Terrestrial Time].
* δΔϵ_1980: (OPTIONAL) Correction in the nutation of the obliquity [rad]
            (**DEFAULT** = 0).
* δΔψ_1980: (OPTIONAL) Correction in the nutation of the longitude [rad]
            (**DEFAULT** = 0).

##### Returns

The rotation that aligns the TOD frame with the GCRF frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

If the EOP correction data is not used, then the GCRF is what is usually called
the J2000 reference frame.

"""
rTODtoGCRF_fk5(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
    rTODtoGCRF_fk5(Matrix, JD_TT, δΔϵ_1980, δΔψ_1980)

function rTODtoGCRF_fk5(::Type{Matrix},
                       JD_TT::Number,
                       δΔϵ_1980::Number = 0,
                       δΔψ_1980::Number = 0)
    # Compute the rotation TOD => MOD.
    D_MOD_TOD  = rTODtoMOD_fk5(Matrix, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Compute the rotation MOD => GCRF.
    D_GCRF_MOD = rMODtoGCRF_fk5(Matrix, JD_TT)

    # Return the full rotation.
    D_GCRF_MOD*D_MOD_TOD
end

function rTODtoGCRF_fk5(::Type{Quaternion},
                       JD_TT::Number,
                       δΔϵ_1980::Number = 0,
                       δΔψ_1980::Number = 0)
    # Compute the rotation TOD => MOD.
    q_MOD_TOD  = rTODtoMOD_fk5(Quaternion, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Compute the rotation MOD => GCRF.
    q_GCRF_MOD = rMODtoGCRF_fk5(Quaternion, JD_TT)

    # Return the full rotation.
    q_MOD_TOD*q_GCRF_MOD
end

"""
### function rGCRFtoTOD_fk5([T,] JD_TT::Number [, δΔϵ_1980::Number, δΔψ_1980::Number])

Compute the rotation that aligns the Geocentric Celestial Reference Frame (GCRF)
with the True of Date (TOD) frame at the Julian Day (Terrestrial Time) `JD_TT`.
This algorithm uses the IAU-76/FK5 theory. Notice that one can provide
corrections for the nutation in obliquity (`δΔϵ`) and in longitude (`δΔψ`) that
are usually obtained from IERS EOP Data (see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT** = `Matrix`).
* JD_TT: Julian Day [Terrestrial Time].
* δΔϵ_1980: (OPTIONAL) Correction in the nutation of the obliquity [rad]
            (**DEFAULT** = 0).
* δΔψ_1980: (OPTIONAL) Correction in the nutation of the longitude [rad]
            (**DEFAULT** = 0).

##### Returns

The rotation that aligns the GCRF frame with the TOD frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

If the EOP correction data is not used, then the GCRF is what is usually called
the J2000 reference frame.

"""

rGCRFtoTOD_fk5(JD_TT::Number, δΔϵ_1980::Number = 0, δΔψ_1980::Number = 0) =
    rTODtoGCRF_fk5(Matrix, JD_TT, δΔϵ_1980, δΔψ_1980)'

function rGCRFtoTOD_fk5(::Type{Matrix},
                        JD_TT::Number,
                        δΔϵ_1980::Number = 0,
                        δΔψ_1980::Number = 0)
    rTODtoGCRF_fk5(Matrix, JD_TT, δΔϵ_1980, δΔψ_1980)'
end

function rGCRFtoTOD_fk5(::Type{Quaternion},
                        JD_TT::Number,
                        δΔϵ_1980::Number = 0,
                        δΔψ_1980::Number = 0)
    conj(rTODtoGCRF_fk5(Quaternion, JD_TT, δΔϵ_1980, δΔψ_1980))
end
