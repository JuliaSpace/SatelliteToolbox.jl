#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Compute the satellite position.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-05-13: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   The functions related to Geodetic coordinates and ECEF frames were removed
#   because they must be updated to use the correct reference frames. This will
#   be done in the future.
#
# 2017-08-04: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   All DEPRECATED functions were removed.
#
# 2016-07-21: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   WARNING: satellite_position_latlon was renamed to satellite_position_LLA.
#   The old function is still present, but it is marked as DEPRECATED and must
#   not be used for new projects.
#
# 2015-11-05: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Remove the deprecated structure OrbitalParameters.
#
# 2014-08-12: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#   Add support to the structure OrbitalParameters.
#   WARNING: the order of parameters in function satellite_position_i changed.
#
# 2014-07-28: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export satellite_position_i

"""
### function satellite_position_i(a::Number, e::Number, i::Number, RAAN::Number, w::Number, f::Number)

Compute the satellite position in the Earth-Centered Inertial (ECI) reference
frame given the orbital elements `a`, `e`, `i`, `RAAN`, `w`, and `f`.

Notice that the ECI frame used will be the same as the frame of the orbital
elements.

##### Args

* a: Semi-major axis.
* e: Eccentricity.
* i: Inclination [rad].
* RAAN: Right ascension of the ascending node [rad].
* w: Argument of perigee [rad].
* f: True anomaly [rad].

##### Returns

* The satellite position vector represented in the ECI reference frame.
* The unit vector perpendicular to the satellite position vector that lies on
  the orbit plane represented in the ECI reference frame.

###### Remarks

The satellite position vector will have the same unit of the semi-major axis.

"""
function satellite_position_i(a::Number, e::Number, i::Number, RAAN::Number,
                              w::Number, f::Number)
    # Compute the radius from the focus.
    norm_r = a*(1-e^2)/(1+e*cos(f))

    # Let s be the coordinate system in which:
    #     - The X axis points towards the satellite;
    #     - The Z axis is normal to the orbit plane (right-hand direction);
    #     - The Y axis completes a right-hand coordinate system.
    #
    # Thus, the satellite vector represented in the s coordinate frame is:
    r_s = SVector{3}(1,0,0)*norm_r

    # rt is the versor perpendicular to the r vector that lies on the orbit
    # plane.
    rt_s = SVector{3}(0,1,0)

    # Compute the matrix that rotates from the S coordinate frame to the
    # Inertial coordinate Frame.
    Dis = angle2dcm(RAAN, i, w+f, :ZXZ)'

    # Compute the satellite vector represented in the Inertial coordinate
    # frame.
    r_i = Dis*r_s

    # Compute unit vector `rt` represented in the Inertial coordinate frame.
    rt_i = Dis*rt_s

    # Return.
    return r_i, rt_i
end
