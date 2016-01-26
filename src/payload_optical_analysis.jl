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

export minimum_swath_grss, minimum_half_FOV_grss, swath_width

"""
### function minimum_half_FOV_grss(h::Real, T::Real, i::Real, To::Integer)

Compute the minimum half FOV of a ground repeating Sun-synchronous (GRSS) orbit
to cover the entire Equator within the revisit interval.

##### Args

* h: Orbit altitude in the Equator [m].
* T: Orbit period [s].
* i: Inclination [rad].
* To: Orbit cycle [days].

##### Returns

* The minimum half FOV [rad].

"""

function minimum_half_FOV_grss(h::Real, T::Real, i::Real, To::Integer)
    adjacent_track_angle_grss(h, T, i, To, 0.0)
end

"""
### function minimum_half_FOV_grss(h::Real, a::Real, e::Real, i::Real, To::Integer)

Compute the minimum half FOV of a ground repeating Sun-synchronous (GRSS) orbit
to cover the entire Equator within the revisit interval.

##### Args

* h: Orbit altitude in the Equator [m].
* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].
* To: Orbit cycle [days].

##### Returns

* The minimum half FOV [rad].

"""

function minimum_half_FOV_grss(h::Real, a::Real, e::Real, i::Real, To::Integer)
    adjacent_track_angle_grss(h, a, e, i, To, 0.0)
end

"""
### function minimum_swath_grss(T::Real, i::Real, To::Integer)

Compute the minimum swath of a ground repeating Sun-synchronous (GRSS) orbit to
cover the entire Equator within the revisit interval.

##### Args

* T: Orbit period [s].
* i: Inclination [rad].
* To: Orbit cycle [days].

##### Returns

* The minimum swath [m].

"""

function minimum_swath_grss(T::Real, i::Real, To::Integer)
    adjacent_track_distance_grss(T, i, To, 0.0)
end

"""
### function minimum_swath_grss(a::Real, e::Real, i::Real, To::Integer)

Compute the minimum swath of a ground repeating Sun-synchronous (GRSS) orbit to
cover the entire Equator within the revisit interval.

##### Args

* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].
* To: Orbit cycle [days].

##### Returns

* The minimum swath [m].

"""

function minimum_swath_grss(a::Real, e::Real, i::Real, To::Integer)
    adjacent_track_distance_grss(a, e, i, To, 0.0)
end

"""
### function swath_width(h::real, HalfFOV::real)

Compute the swath width given the orbit altitude and the half FOV.

##### Args

* h: Orbit altitude [m].
* HalfFOV: Half field of view [rad].

##### Returns

* The swath width [m].

"""

function swath_width(h::Real, HalfFOV::Real)
    gamma = pi - asin((R0+h)/R0*sin(HalfFOV))
    alpha = pi - gamma - HalfFOV
    S = R0*alpha
end


