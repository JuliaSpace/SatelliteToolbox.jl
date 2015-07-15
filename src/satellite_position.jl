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
# 2014-08-12: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Add support to the structure OrbitalParameters.
#    WARNING: the order of parameters in function satellite_position_i changed.
#
# 2014-07-28: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

import Rotations: angle2dcm!

export  satellite_position_i, sun_position_i

#==#
# 
# @brief Compute the satellite position on the Inertial coordinate frame.
#
# @param[in] sat Orbital parameters.
#
# @return The satellite vector and the versor that is perpendicular to the
# satellite vector and lies on the orbit plane, both represented in the Inertial
# coordinate frame.
#
# @note The dimension of the satellite vector will be the same of that of the
# semi-major axis.
#
#==#

function satellite_position_i{T}(sat::OrbitalParameters{T})
  return satellite_position_i(sat.a, sat.e, sat.i, sat.RAAN, sat.w, sat.f)
end

#==#
# 
# @brief Compute the satellite position on the Inertial coordinate frame.
#
# @param[in] a Semi-major axis.
# @param[in] e Eccentricity.
# @param[in] i Inclination [rad].
# @param[in] RAAN Right ascension of the ascending node [rad].
# @param[in] w Argument of perigee [rad].
# @param[in] f True anomaly [rad].
#
# @return The satellite vector and the versor that is perpendicular to the
# satellite vector and lies on the orbit plane, both represented in the Inertial
# coordinate frame.
#
# @note The dimension of the satellite vector will be the same of that of the
# semi-major axis.
#
#==#

function satellite_position_i{T}(a::T, e::T, i::T, RAAN::T, w::T, f::T)
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
    Dsi = Array(T,(3,3))
    angle2dcm!(Dsi, RAAN, i, w+f, "ZXZ")

    # Compute the satellite vector represented in the Inertial coordinate
    # frame.
    r_i = Dsi'*r_s

    # Compute versor rt represented in the Inertial coordinate frame.
    rt_i = Dsi'*rt_s

    # Return.
    return r_i, rt_i
end
