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
#    Many functions to perform a preliminary analysis of an optical payload.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2015-07-15: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

#==#
#
# @brief Compute the minimum half FOV of a ground repeating sun-synchronous
# (GRSS) orbit to cover the entire Equator within the revisit interval.
#
# @param[in] h Orbit altitude in the Equator [m].
# @param[in] T Orbit period [s].
# @param[in] i Inclination [rad].
# @param[in] To Orbit cycle [days].
#
# @return The minimum half FOV [rad].
#
#==#

function minimum_half_FOV_grss(h::Real, T::Real, i::Real, To::Integer)
    adjacent_track_angle_grss(h, T, i, To, 0.0)
end

#==#
#
# @brief Compute the minimum half FOV of a ground repeating sun-synchronous
# (GRSS) orbit to cover the entire Equator within the revisit interval.
#
# @param[in] h Orbit altitude in the Equator [m].
# @param[in] a Semi-major axis [m].
# @param[in] e Eccentricity [s].
# @param[in] i Inclination [rad].
# @param[in] To Orbit cycle [days].
#
# @return The minimum half FOV [rad].
#
#==#

function minimum_half_FOV_grss(h::Real, a::Real, e::Real, i::Real, To::Integer)
    adjacent_track_angle_grss(h, a, e, i, To, 0.0)
end
    
#==#
#
# @brief Compute the minimum swath of a ground repeating sun-synchronous (GRSS)
# orbit to cover the entire Equator within the revisit interval.
#
# @param[in] T Orbit period [s].
# @param[in] i Inclination [rad].
# @param[in] To Orbit cycle [days].
#
# @return The minimum swath [m].
#
#==#

function minimum_swath_grss(T::Real, i::Real, To::Integer)
    adjacent_track_distance_grss(T, i, To, 0.0)
end

#==#
#
# @brief Compute the minimum swath of a ground repeating sun-synchronous (GRSS)
# orbit to cover the entire Equator within the revisit interval.
#
# @param[in] a Semi-major axis [m].
# @param[in] e Eccentricity [s].
# @param[in] i Inclination [rad].
# @param[in] To Orbit cycle [days].
#
# @return The minimum swath [m].
#
#==#

function minimum_swath_orbit_grss(a::Real, e::Real, i::Real, To::Integer)
    adjacent_track_distance_grss(a, e, i, To, 0.0)
end

#==#
#
# @brief Compute the swath width given the orbit altitude and the half FOV.
#
# @param[in] h Orbit altitude [m].
# @param[in] HalfFOV Half field of view [rad].
#
#==#

function swath_width(h::Real, HalfFOV::Real)
    gamma = pi - asin((R0+h)/R0 + sin(HalfFOV))
    alpha = pi - gamma - HalfFOV
    S = R0*alpha
end


