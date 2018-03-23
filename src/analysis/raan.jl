# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
#    Functions related to the analysis of the Right Ascension of the Ascending
#    Node (RAAN).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2017-08-04: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export compute_RAAN_lt, sim_RAAN_J2

"""
### function compute_RAAN_lt(t0::Int, asc_node_lt::Real)

Compute the RAAN given a data and a local time.

##### Args

* t0: Launch date [number of days since 01/01/2000].
* asc_node_lt: Desired local time for the ascending node [hour].

##### Returns

* The RAAN in the interval [0, 2π].

##### Remarks

The sun position is computed at noon of the day t0.

"""

function compute_RAAN_lt(t0::Int, asc_node_lt::Real)
    println("WARNING: The function compute_RAAN_lt(t0::Int, asc_node_lt::Real) is deprecated!")
    println("Use the function compute_RAAN_lt(JD::Real, asc_node_lt::Real) instead.\n")

    JD = JD_J2000 + t0

    compute_RAAN_lt(JD, asc_node_lt)
end

"""
### function compute_RAAN_lt(JD::Real, asc_node_lt::Real)

Compute the RAAN given a data and a local time.

##### Args

* JD: Julian day.
* asc_node_lt: Desired local time for the ascending node [hour].

##### Returns

* The RAAN in the interval [0, 2π].

"""

function compute_RAAN_lt(JD::Real, asc_node_lt::Real)
    # Get the sun position at noon (UT) represented in the Inertial ref. frame.
    s_i = sun_position_i(JD)

    # Get the desired angle between the Sun and the ascending node [deg].
    alpha = (asc_node_lt-12.0)*float(pi)/12.0

    # Get the right ascension of the Sun in the Inertial ref. frame. This is the
    # Sun apparent local time.
    SALT = atan2(s_i[2],s_i[1])

    # Get the equation of time to compute the Sun mean local time [rad].
    eot = equation_of_time(JD)

    # Compute the Sun mean local time.
    SMLT = SALT + eot

    # Compute the desired RAAN in the interval 0, 2*pi.
    RAAN = mod(SMLT+alpha, 2*pi)
end

"""
### function sim_RAAN_J2(a::Real, e::Real, i::Real, RAAN_0::Real, numDays::Integer)

Simulate the RAAN of an orbit considering J2 perturbations.

##### Args

* a: Semi-major axis of the satellite orbit [m].
* e: Orbit eccentricity.
* i: Orbit inclination [rad].
* RAAN_0: Initial right ascension of the ascending node [rad].
* numDays: Number of days of the analysis.

##### Returns

* The RAAN computed for each day in radians (0-2\pi).

"""

function sim_RAAN_J2(a::Real,
                     e::Real,
                     i::Real,
                     RAAN_0::Real,
                     numDays::Integer)

    # Initialization of variables.
    days = collect(0:1:numDays-1) # Vector of the days in which the RAAN will be
                                  # simulated.

    # Output vector.
    RAAN = Array{Float64}(numDays,1)

    # RAAN rotation rate [rad/day].
    dOmega = dRAAN(a, e, i, :J2)*24.0*3600.0

    # Simulate the RAAN for each day considering just the J2 perturbations.
    RAAN = mod(RAAN_0 + dOmega.*days,2*pi)

    [days RAAN]
end

