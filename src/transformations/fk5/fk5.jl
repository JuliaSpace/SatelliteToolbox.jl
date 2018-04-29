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
#   [2] Gontier, A. M., Capitaine, N (1991). High-Accuracy Equation of Equinoxes
#       and VLBI Astrometric Modelling. Radio Interferometry: Theory, Techniques
#       and Applications, IAU Coll. 131, ASP Conference Series, Vol. 19.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-04-09: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export rITRFtoPEF_fk5, rPEFtoITRF_fk5
export rPEFtoTOD_fk5,  rTODtoPEF_fk5
export rTODtoMOD_fk5,  rMODtoTOD_fk5
export rMODtoGCRF_fk5, rGCRFtoMOD_fk5

export rITRFtoGCRF_fk5, rGCRFtoITRF_fk5
export rPEFtoGCRF_fk5, rGCRFtoPEF_fk5
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

#                                 ITRF <=> PEF
# ==============================================================================

"""
### function rITRFtoPEF_fk5([T,] x_p::Number, y_p::Number)

Compute the rotation that aligns the International Terrestrial Reference Frame
(ITRF) with the Pseudo-Earth Fixed (PEF) frame considering the polar motion
represented by the angles `x_p` and `y_p` that is obtained from IERS EOP Data
(see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT** = `Matrix`).
* x_p: Polar motion displacement about X-axis, which is the IERS Reference
       Meridian direction (positive south along the 0˚ longitude meridian)
       [rad].
* y_p: Polar motion displacement about Y-axis (90˚W or 270˚E meridian) [rad].

##### Returns

The rotation that aligns the ITRF frame with the PEF frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

The ITRF is defined based on the International Reference Pole (IRP), which is
the location of the terrestrial pole agreed by international committees [1]. The
Pseudo-Earth Fixed, on the other hand, is defined based on the Earth axis of
rotation, or the Celestial Intermediate Pole (CIP). Hence, PEF XY-plane contains
the True Equator. Furthermore, since the recovered latitude and longitude are
sensitive to the CIP, then it should be computed considering the PEF frame.

"""

rITRFtoPEF_fk5(x_p::Number, y_p::Number) = rITRFtoPEF_fk5(Matrix, x_p, y_p)

function rITRFtoPEF_fk5(::Type{Matrix}, x_p::Number, y_p::Number)
    # Notice that `x_p` and `y_p` are displacements in X and Y directions and
    # **not** rotation angles. Hence, a displacement in X is a rotation in Y and
    # a displacement in Y is a rotation in X.
    [  1    0  -x_p;
       0    1  +y_p;
     +x_p -y_p   1 ;]
end

function rITRFtoPEF_fk5(::Type{Quaternion}, x_p::Number, y_p::Number)
    # Notice that `x_p` and `y_p` are displacements in X and Y directions and
    # **not** rotation angles. Hence, a displacement in X is a rotation in Y and
    # a displacement in Y is a rotation in X.
    q = Quaternion([1; +y_p/2; +x_p/2; 0])
    q/norm(q)
end

"""
### function rPEFtoITRF_fk5([T,] x_p::Number, y_p::Number)

Compute the rotation that aligns the Pseudo-Earth Fixed (PEF) with the
International Terrestrial Reference Frame (ITRF) considering the polar motion
represented by the angles `x_p` and `y_p` that is obtained from IERS EOP Data
(see `get_iers_eop`).

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT** = `Matrix`).
* x_p: Polar motion rotation about X-axis, which is the IERS Reference Meridian
       direction (positive south along the 0˚ longitude meridian). [rad].
* y_p: Polar motion rotation about Y-axis (90˚ W or 270˚ E meridian) [rad].

##### Returns

The rotation that aligns the PEF frame with the ITRF. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

The ITRF is defined based on the International Reference Pole (IRP), which is
the location of the terrestrial pole agreed by international committees [1]. The
Pseudo-Earth Fixed, on the other hand, is defined based on the Earth axis of
rotation, or the Celestial Intermediate Pole (CIP). Hence, PEF XY-plane contains
the True Equator. Furthermore, since the recovered latitude and longitude are
sensitive to the CIP, then it should be computed considering the PEF frame.

"""

rPEFtoITRF_fk5(x_p::Number, y_p::Number) = rPEFtoITRF_fk5(Matrix, x_p, y_p)

function rPEFtoITRF_fk5(::Type{Matrix}, x_p::Number, y_p::Number)
    # Notice that `x_p` and `y_p` are displacements in X and Y directions and
    # **not** rotation angles. Hence, a displacement in X is a rotation in Y and
    # a displacement in Y is a rotation in X.
    [  1    0  +x_p;
       0    1  -y_p;
     -x_p +y_p   1 ;]
end

function rPEFtoITRF_fk5(::Type{Quaternion}, x_p::Number, y_p::Number)
    # Notice that `x_p` and `y_p` are displacements in X and Y directions and
    # **not** rotation angles. Hence, a displacement in X is a rotation in Y and
    # a displacement in Y is a rotation in X.
    q = Quaternion([1; -y_p/2; -x_p/2; 0])
    q/norm(q)
end

#                                 PEF <=> TOD
# ==============================================================================

"""
### function rPEFtoTOD_fk5([T,] JD_UT1::Number, JD_TT::Number [, δΔψ_1980::Number])

Compute the rotation that aligns the Pseudo-Earth Fixed (PEF) frame with the
True of Date (TOD) frame at the Julian Day `JD_UT1` (UT1) and `JD_TT`
(Terrestrial Time). This algorithm uses the IAU-76/FK5 theory. Notice that one
can provide correction for the nutation in longitude (`δΔψ`) that are usually
obtained from IERS EOP Data (see `get_iers_eop`).

The Julian Day in UT1 is used to compute the Greenwich Mean Sidereal Time (GMST)
(see `JDtoGMST`), whereas the Julian Day in Terrestrial Time is used to compute
the nutation in the longitude. Notice that the Julian Day in UT1 and in
Terrestrial Time must be equivalent, i.e. must be related to the same instant.
This function **does not** check this.

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT** = `Matrix`).
* JD_UT1: Julian Day [UT1].
* JD_TT: Julian Day [Terrestrial Time].
* δΔψ_1980: (OPTIONAL) Correction in the nutation of the longitude [rad]
            (**DEFAULT** = 0).

##### Returns

The rotation that aligns the PEF frame with the TOD frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

The Pseudo-Earth Fixed (PEF) frame is rotated into the True of Date (TOD) frame
considering the 1980 IAU Theory of Nutation. The IERS EOP corrections must be
added if one wants to make the rotation consistent with the Geocentric Celestial
Reference Systems (GCRS).

"""

rPEFtoTOD_fk5(JD_UT1::Number, JD_TT::Number, δΔψ_1980::Number = 0) =
    rPEFtoTOD_fk5(Matrix, JD_UT1, JD_TT, δΔψ_1980)

function rPEFtoTOD_fk5(::Type{Matrix},
                       JD_UT1::Number,
                       JD_TT::Number,
                       δΔψ_1980::Number = 0)
    # Compute the nutation in the Julian Day (Terrestrial Time) `JD_TT`.
    (mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(JD_TT)

    # Add the corrections to the nutation in obliquity and longitude.
    Δψ_1980 += δΔψ_1980

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

    # Compute the Mean Greenwich Sidereal Time.
    θ_gmst = JDtoGMST(JD_UT1)

    # Compute the Greenwich Apparent Sidereal Time (GAST).
    #
    # TODO: Should GAST be moved to a new function as the GMST?
    θ_gast = θ_gmst + Eq_equinox1982

    # Compute the rotation matrix.
    create_rotation_matrix(-θ_gast, 'Z')
end

function rPEFtoTOD_fk5(::Type{Quaternion},
                       JD_UT1::Number,
                       JD_TT::Number,
                       δΔψ_1980::Number = 0)
    # Compute the nutation in the Juliay Day (Terrestrial Time) `JD_TT`.
    (mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(JD_TT)

    # Add the corrections to the nutation in obliquity and longitude.
    Δψ_1980 += δΔψ_1980

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

    # Compute the Mean Greenwich Sidereal Time.
    θ_gmst = JDtoGMST(JD_UT1)

    # Compute the Greenwich Apparent Sidereal Time (GAST).
    #
    # TODO: Should GAST be moved to a new function as the GMST?
    θ_gast = θ_gmst + Eq_equinox1982

    # Compute the quaternion.
    Quaternion([cos(θ_gast/2); 0; 0; -sin(θ_gast/2)])
end

"""
### function rTODtoPEF_fk5([T,] JD_UT1::Number, JD_TT::Number [, δΔψ_1980::Number])

Compute the rotation that aligns the True of Date (TOD) frame with the
Pseudo-Earth Fixed (PEF) frame at the Julian Day `JD_UT1` (UT1) and `JD_TT`
(Terrestrial Time). This algorithm uses the IAU-76/FK5 theory. Notice that one
can provide correction for the nutation in longitude (`δΔψ`) that are usually
obtained from IERS EOP Data (see `get_iers_eop`).

The Julian Day in UT1 is used to compute the Greenwich Mean Sidereal Time (GMST)
(see `JDtoGMST`), whereas the Julian Day in Terrestrial Time is used to compute
the nutation in the longitude. Notice that the Julian Day in UT1 and in
Terrestrial Time must be equivalent, i.e. must be related to the same instant.
This function **does not** check this.

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT** = `Matrix`).
* JD_UT1: Julian Day [UT1].
* JD_TT: Julian Day [Terrestrial Time].
* δΔψ_1980: (OPTIONAL) Correction in the nutation of the longitude [rad]
            (**DEFAULT** = 0).

##### Returns

The rotation that aligns the TOD frame with the PEF frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

The True of Date (TOD) frame is rotated into the Pseudo-Earth Fixed (PEF) frame
considering the 1980 IAU Theory of Nutation. The IERS EOP corrections must be
added if one wants to make the rotation consistent with the Geocentric Celestial
Reference Systems (GCRS).

"""

rTODtoPEF_fk5(JD_UT1::Number, JD_TT::Number, δΔψ_1980::Number = 0) =
    rPEFtoTOD_fk5(Matrix, JD_UT1, JD_TT, δΔψ_1980)'

function rTODtoPEF_fk5(::Type{Matrix},
                       JD_UT1::Number,
                       JD_TT::Number,
                       δΔψ_1980::Number = 0)
    rPEFtoTOD_fk5(Matrix, JD_UT1, JD_TT, δΔψ_1980)'
end

function rTODtoPEF_fk5(::Type{Quaternion},
                       JD_UT1::Number,
                       JD_TT::Number,
                       δΔψ_1980::Number = 0)
    conj(rPEFtoTOD_fk5(Quaternion, JD_UT1, JD_TT, δΔψ_1980))
end

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

    # Compute the nutation in the Julian Day (Terrestrial Time) `JD_TT`.
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
    # Compute the nutation in the Julian Day (Terrestrial Time) `JD_TT`.
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

#                                ITRF <=> GCRF
# ==============================================================================

"""
### function rITRFtoGCRF_fk5([T,] JD_UT1::Number, JD_TT::Number, x_p::Number, y_p::Number [, δΔϵ_1980::Number, δΔψ_1980::Number])

Compute the rotation that aligns the International Terrestrial Reference Frame
(ITRF) with the Pseudo-Earth Fixed (PEF) frame at the Julian Day `JD_UT1` (UT1)
and `JD_TT` (Terrestrial Time), and considering the IERS EOP Data `x_p`, `y_p`,
`δΔϵ_1980`, and `δΔψ_1980` (see `get_iers_eop`). This algorithm uses the
IAU-76/FK5 theory.

The Julian Day in UT1 is used to compute the Greenwich Mean Sidereal Time (GMST)
(see `JDtoGMST`), whereas the Julian Day in Terrestrial Time is used to compute
the nutation in the longitude. Notice that the Julian Day in UT1 and in
Terrestrial Time must be equivalent, i.e. must be related to the same instant.
This function **does not** check this.

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT** = `Matrix`).
* JD_UT1: Julian Day [UT1].
* JD_TT: Julian Day [Terrestrial Time].
* x_p: Polar motion displacement about X-axis, which is the IERS Reference
       Meridian direction (positive south along the 0˚ longitude meridian)
       [rad].
* y_p: Polar motion displacement about Y-axis (90˚W or 270˚E meridian) [rad].
* δΔϵ_1980: (OPTIONAL) Correction in the nutation of the obliquity [rad]
            (**DEFAULT** = 0).
* δΔψ_1980: (OPTIONAL) Correction in the nutation of the longitude [rad]
            (**DEFAULT** = 0).

##### Returns

The rotation that aligns the ITRF frame with the GCRF frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

The EOP data related to the polar motion (`x_p` and `y_p`) is required, since
this is the only way available to compute the conversion ITRF <=> PEF (the
models are highly imprecise since the motion is still not very well understood
[1]). However, the EOP data related to the nutation of the obliquity
(`δΔϵ_1980`) and the nutation of the longitude (`δΔψ_1980`) can be omitted. In
this case, the GCRF frame is what is usually called J2000 reference frame.

"""

rITRFtoGCRF_fk5(JD_UT1::Number,
                JD_TT::Number,
                x_p::Number,
                y_p::Number,
                δΔϵ_1980::Number = 0,
                δΔψ_1980::Number = 0) =
    rITRFtoGCRF_fk5(Matrix, JD_UT1, JD_TT, x_p, y_p, δΔϵ_1980, δΔψ_1980)

function rITRFtoGCRF_fk5(::Type{Matrix},
                         JD_UT1::Number,
                         JD_TT::Number,
                         x_p::Number,
                         y_p::Number,
                         δΔϵ_1980::Number = 0,
                         δΔψ_1980::Number = 0)
    # Compute the rotation ITRF => PEF.
    D_PEF_ITRF = rITRFtoPEF_fk5(Matrix, x_p, y_p)

    # Compute the rotation PEF => GCRF.
    D_GCRF_PEF = rPEFtoGCRF_fk5(Matrix, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Return the full rotation.
    D_GCRF_PEF*D_PEF_ITRF
end

function rITRFtoGCRF_fk5(::Type{Quaternion},
                         JD_UT1::Number,
                         JD_TT::Number,
                         x_p::Number,
                         y_p::Number,
                         δΔϵ_1980::Number = 0,
                         δΔψ_1980::Number = 0)
    # Compute the rotation ITRF => PEF.
    q_PEF_ITRF = rITRFtoPEF_fk5(Quaternion, x_p, y_p)

    # Compute the rotation PEF => GCRF.
    q_GCRF_PEF = rPEFtoGCRF_fk5(Quaternion, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Return the full rotation.
    q_PEF_ITRF*q_GCRF_PEF
end

"""
### function rGCRFtoITRF_fk5([T,] JD_UT1::Number, JD_TT::Number, x_p::Number, y_p::Number [, δΔϵ_1980::Number, δΔψ_1980::Number])

Compute the rotation that aligns the Pseudo-Earth Fixed (PEF) frame with the
International Terrestrial Reference Frame (ITRF) at the Julian Day `JD_UT1`
(UT1) and `JD_TT` (Terrestrial Time), and considering the IERS EOP Data `x_p`,
`y_p`, `δΔϵ_1980`, and `δΔψ_1980` (see `get_iers_eop`). This algorithm uses the
IAU-76/FK5 theory.

The Julian Day in UT1 is used to compute the Greenwich Mean Sidereal Time (GMST)
(see `JDtoGMST`), whereas the Julian Day in Terrestrial Time is used to compute
the nutation in the longitude. Notice that the Julian Day in UT1 and in
Terrestrial Time must be equivalent, i.e. must be related to the same instant.
This function **does not** check this.

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT** = `Matrix`).
* JD_UT1: Julian Day [UT1].
* JD_TT: Julian Day [Terrestrial Time].
* x_p: Polar motion displacement about X-axis, which is the IERS Reference
       Meridian direction (positive south along the 0˚ longitude meridian)
       [rad].
* y_p: Polar motion displacement about Y-axis (90˚W or 270˚E meridian) [rad].
* δΔϵ_1980: (OPTIONAL) Correction in the nutation of the obliquity [rad]
            (**DEFAULT** = 0).
* δΔψ_1980: (OPTIONAL) Correction in the nutation of the longitude [rad]
            (**DEFAULT** = 0).

##### Returns

The rotation that aligns the GCRF frame with the ITRF frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

The EOP data related to the polar motion (`x_p` and `y_p`) is required, since
this is the only way available to compute the conversion ITRF <=> PEF (the
models are highly imprecise since the motion is still not very well understood
[1]). However, the EOP data related to the nutation of the obliquity
(`δΔϵ_1980`) and the nutation of the longitude (`δΔψ_1980`) can be omitted. In
this case, the GCRF frame is what is usually called J2000 reference frame.

"""

rGCRFtoITRF_fk5(JD_UT1::Number,
                JD_TT::Number,
                x_p::Number,
                y_p::Number,
                δΔϵ_1980::Number = 0,
                δΔψ_1980::Number = 0) =
    rITRFtoGCRF_fk5(Matrix, JD_UT1, JD_TT, x_p, y_p, δΔϵ_1980, δΔψ_1980)'

function rGCRFtoITRF_fk5(::Type{Matrix},
                         JD_UT1::Number,
                         JD_TT::Number,
                         x_p::Number,
                         y_p::Number,
                         δΔϵ_1980::Number = 0,
                         δΔψ_1980::Number = 0)
    rITRFtoGCRF_fk5(Matrix, JD_UT1, JD_TT, x_p, y_p, δΔϵ_1980, δΔψ_1980)'
end

function rGCRFtoITRF_fk5(::Type{Quaternion},
                         JD_UT1::Number,
                         JD_TT::Number,
                         x_p::Number,
                         y_p::Number,
                         δΔϵ_1980::Number = 0,
                         δΔψ_1980::Number = 0)
    conj(rITRFtoGCRF_fk5(Quaternion,
                         JD_UT1,
                         JD_TT,
                         x_p,
                         y_p,
                         δΔϵ_1980,
                         δΔψ_1980))
end

#                                 PEF <=> GCRF
# ==============================================================================

"""
### function rPEFtoGCRF_fk5([T,] JD_UT1::Number, JD_TT::Number [, δΔϵ_1980::Number, δΔψ_1980::Number])

Compute the rotation that aligns the Pseudo-Earth Fixed (PEF) frame with the
Geocentric Celestial Reference Frame (GCRF) at the Julian Day `JD_UT1` (UT1) and
`JD_TT` (Terrestrial Time). This algorithm uses the IAU-76/FK5 theory. Notice
that one can provide corrections for the nutation in obliquity (`δΔϵ`) and in
longitude (`δΔψ`) that are usually obtained from IERS EOP Data (see
`get_iers_eop`).

The Julian Day in UT1 is used to compute the Greenwich Mean Sidereal Time (GMST)
(see `JDtoGMST`), whereas the Julian Day in Terrestrial Time is used to compute
the nutation in the longitude. Notice that the Julian Day in UT1 and in
Terrestrial Time must be equivalent, i.e. must be related to the same instant.
This function **does not** check this.

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT** = `Matrix`).
* JD_UT1: Julian Day [UT1].
* JD_TT: Julian Day [Terrestrial Time].
* δΔϵ_1980: (OPTIONAL) Correction in the nutation of the obliquity [rad]
            (**DEFAULT** = 0).
* δΔψ_1980: (OPTIONAL) Correction in the nutation of the longitude [rad]
            (**DEFAULT** = 0).

##### Returns

The rotation that aligns the PEF frame with the GCRF frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

If the EOP correction data is not used, then the GCRF is what is usually called
the J2000 reference frame.

"""
rPEFtoGCRF_fk5(JD_UT1::Number,
               JD_TT::Number,
               δΔϵ_1980::Number = 0,
               δΔψ_1980::Number = 0) =
    rPEFtoGCRF_fk5(Matrix, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)

function rPEFtoGCRF_fk5(::Type{Matrix},
                       JD_UT1::Number,
                       JD_TT::Number,
                       δΔϵ_1980::Number = 0,
                       δΔψ_1980::Number = 0)

    # Notice that, in this case, we will not use `rPEFtoTOD` and `rTODtoMOD`
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

    # Compute the Mean Greenwich Sidereal Time.
    θ_gmst = JDtoGMST(JD_UT1)

    # Compute the Greenwich Apparent Sidereal Time (GAST).
    #
    # TODO: Should GAST be moved to a new function as the GMST?
    θ_gast = θ_gmst + Eq_equinox1982

    # Compute the rotation PEF => TOD.
    D_TOD_PEF = create_rotation_matrix(-θ_gast, 'Z')

    # Compute the rotation TOD => MOD.
    D_MOD_TOD = angle2dcm(ϵ_1980, Δψ_1980, -mϵ_1980, "XZX")

    # Compute the rotation MOD => GCRF.
    D_GCRF_MOD = rMODtoGCRF_fk5(Matrix, JD_TT)

    # Return the full rotation.
    D_GCRF_MOD*D_MOD_TOD*D_TOD_PEF
end

function rPEFtoGCRF_fk5(::Type{Quaternion},
                       JD_UT1::Number,
                       JD_TT::Number,
                       δΔϵ_1980::Number = 0,
                       δΔψ_1980::Number = 0)

    # Notice that, in this case, we will not use `rPEFtoTOD` and `rTODtoMOD`
    # because this would call the function `nutation` twice, leading to a huge
    # performance drop. Hence, the code of those two functions is almost
    # entirely rewritten here.

    # Compute the nutation in the Juliay Day (Terrestrial Time) `JD_TT`.
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

    # Compute the Mean Greenwich Sidereal Time.
    θ_gmst = JDtoGMST(JD_UT1)

    # Compute the Greenwich Apparent Sidereal Time (GAST).
    #
    # TODO: Should GAST be moved to a new function as the GMST?
    θ_gast = θ_gmst + Eq_equinox1982

    # Compute the rotation PEF => TOD.
    q_TOD_PEF = Quaternion([cos(θ_gast/2); 0; 0; -sin(θ_gast/2)])

    # Compute the rotation TOD => MOD.
    q1_x      = Quaternion([ cos(ϵ_1980/2);   sin(ϵ_1980/2); 0;        0      ])
    q2_z      = Quaternion([cos(Δψ_1980/2);        0;        0; sin(Δψ_1980/2)])
    q3_x      = Quaternion([cos(mϵ_1980/2); -sin(mϵ_1980/2); 0;        0      ])
    q_MOD_TOD = q1_x*q2_z*q3_x

    # Compute the rotation MOD => GCRF.
    q_GCRF_MOD = rMODtoGCRF_fk5(Quaternion, JD_TT)

    # Return the full rotation.
    q_TOD_PEF*q_MOD_TOD*q_GCRF_MOD
end

"""
### function rGCRFtoPEF_fk5([T,] JD_UT1::Number, JD_TT::Number [, δΔϵ_1980::Number, δΔψ_1980::Number])

Compute the rotation that aligns the Geocentric Celestial Reference Frame (GCRF)
with the Pseudo-Earth Fixed (PEF) frame at the Julian Day `JD_UT1` (UT1) and
`JD_TT` (Terrestrial Time). This algorithm uses the IAU-76/FK5 theory. Notice
that one can provide corrections for the nutation in obliquity (`δΔϵ`) and in
longitude (`δΔψ`) that are usually obtained from IERS EOP Data (see
`get_iers_eop`).

The Julian Day in UT1 is used to compute the Greenwich Mean Sidereal Time (GMST)
(see `JDtoGMST`), whereas the Julian Day in Terrestrial Time is used to compute
the nutation in the longitude. Notice that the Julian Day in UT1 and in
Terrestrial Time must be equivalent, i.e. must be related to the same instant.
This function **does not** check this.

The rotation type is described by the optional variable `T`. If it is `Matrix`,
then a DCM will be returned. Otherwise, if it is `Quaternion`, then a Quaternion
will be returned. In case this parameter is omitted, then it falls back to
`Matrix`.

##### Args

* T: (OPTIONAL) Type of the rotation representation (**DEFAULT** = `Matrix`).
* JD_UT1: Julian Day [UT1].
* JD_TT: Julian Day [Terrestrial Time].
* δΔϵ_1980: (OPTIONAL) Correction in the nutation of the obliquity [rad]
            (**DEFAULT** = 0).
* δΔψ_1980: (OPTIONAL) Correction in the nutation of the longitude [rad]
            (**DEFAULT** = 0).

##### Returns

The rotation that aligns the GCRF frame with the PEF frame. The rotation
representation is selected by the optional parameter `T`.

##### Remarks

If the EOP correction data is not used, then the GCRF is what is usually called
the J2000 reference frame.

"""

rGCRFtoPEF_fk5(JD_UT1::Number,
               JD_TT::Number,
               δΔϵ_1980::Number = 0,
               δΔψ_1980::Number = 0) =
    rPEFtoGCRF_fk5(Matrix, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)'

function rGCRFtoPEF_fk5(::Type{Matrix},
                       JD_UT1::Number,
                       JD_TT::Number,
                       δΔϵ_1980::Number = 0,
                       δΔψ_1980::Number = 0)
    rPEFtoGCRF_fk5(Matrix, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)'
end

function rGCRFtoPEF_fk5(::Type{Quaternion},
                       JD_UT1::Number,
                       JD_TT::Number,
                       δΔϵ_1980::Number = 0,
                       δΔψ_1980::Number = 0)
    conj(rPEFtoGCRF_fk5(Quaternion, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980))
end

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
