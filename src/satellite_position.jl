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
#    Compute the satellite position.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2015-11-05: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#    Remove the deprecated structure OrbitalParameters.
#
# 2014-08-12: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Add support to the structure OrbitalParameters.
#    WARNING: the order of parameters in function satellite_position_i changed.
#
# 2014-07-28: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

import Rotations: angle2dcm!

export satellite_position_i

"""
### function satellite_position_i(a::Real, e::Real, i::Real, RAAN::Real, w::Real, f::Real)

Compute the satellite position on the Inertial coordinate frame.

##### Args

* a: Semi-major axis.
* e: Eccentricity.
* i: Inclination [rad].
* RAAN: Right ascension of the ascending node [rad].
* w: Argument of perigee [rad].
* f: True anomaly [rad].

##### Returns

* The satellite position vector represented in the Inertial coordinate frame.
* The versor perpendicular to the satellite position vector that lies on the
orbit plane represented in the Inertial coordinate frame.

###### Remarks

The satellite position vector will have the same unit of the semi-major axis.

"""

function satellite_position_i(a::Real, e::Real, i::Real, RAAN::Real,
                              w::Real, f::Real)
    # Compute the radius from the focus.
    norm_r = a*(1-e^2)/(1+e*cos(f))

    # Let s be the coordinate system in which:
    #     - The X axis points towards the satellite;
    #     - The Z axis is normal to the orbit plane (right-hand direction);
    #     - The Y axis completes a right-hand coordinate system.
    #
    # Thus, the satellite vector represented in the s coordinate frame is:
    r_s = [1;0;0]*norm_r

    # rt is the versor perpendicular to the r vector that lies on the orbit
    # plane.
    rt_s = [0;1;0]

    # Compute the matrix that rotates from the S coordinate frame to the
    # Inertial coordinate Frame.
    Dsi = Array(Float64,(3,3))
    angle2dcm!(Dsi, RAAN, i, w+f, "ZXZ")

    # Compute the satellite vector represented in the Inertial coordinate
    # frame.
    r_i = Dsi'*r_s

    # Compute versor rt represented in the Inertial coordinate frame.
    rt_i = Dsi'*rt_s

    # Return.
    return r_i, rt_i
end
