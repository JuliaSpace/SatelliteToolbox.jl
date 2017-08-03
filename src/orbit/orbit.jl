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
#    Many auxiliary functions for orbit computations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2014-12-18: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export compute_RAAN_lt, sim_RAAN_J2
export dRAAN_J2, dw_J2, n_J0, n_J2, t_J0, t_J2

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
    Si = sun_position_i(JD)

    # Get the desired angle between the Sun and the ascending node [deg].
    alpha = (asc_node_lt-12.0)*float(pi)/12.0

    # Get the ascension of the Sun in the Inertial ref. frame.
    S_asc_i = atan2(Si[2],Si[1])

    # Compute the desired RAAN in the interval 0, 2*pi.
    RAAN = mod(S_asc_i+alpha, 2*float(pi))
end

"""
### function dRAAN_J2(a::Real, e::Real, i::Real)

Compute the perturbation of the RAAN using terms up to J2.

##### Args

* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].

##### Returns

* The perturbation of the RAAN (J2) [rad/s].

"""

function dRAAN_J2(a::Real, e::Real, i::Real)
    # Check if the perigee is inside Earth.
    if ( !is_orbit_valid(a,e) )
        throw(OrbitInvalidPerigee(a*(1.0-e)))
    end

    # Semi-lactum rectum.
    p = a*(1.0-e^2)

    # Unperturbed orbit period.
    n0 = n_J0(a)

    # Perturbation of the right ascension of the ascending node.
    -3.0/2.0*R0^2/(p^2)*n0*J2*cos(i)
end

"""
### function dw_J2(a::Real, e::Real, i::Real)

Compute the perturbation of the argument of perigee using terms up to J2.

##### Args

* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].

##### Returns

* The perturbation of the argument of perigee (J2) [rad/s].

"""

function dw_J2(a::Real, e::Real, i::Real)
    # Check if the perigee is inside Earth.
    if ( !is_orbit_valid(a,e) )
        throw(OrbitInvalidPerigee(a*(1.0-e)))
    end

    # Semi-lactum rectum.
    p = a*(1.0-e^2)

    # Unperturbed orbit period.
    n0 = n_J0(a)

    # Perturbation of the argument of perigee.
    3.0*R0^2*J2/(4.0*p^2)*n0*(5.0*cos(i)^2-1.0)
end

"""
### function is_orbit_valid(a::Real, e::Real)

Verify if the orbit is valid.

##### Args

* a: Semi-major axis [m].
* e: Eccentricity.

##### Returns

* **TRUE**: The orbit is valid.
* **FALSE**: The orbit is invalid.

"""

function is_orbit_valid(a::Real, e::Real)
    # Check if the arguments are valid.
    if ( a < 0. )
        throw(ArgumentError("The semi-major axis must be greater than 0."))
    end

    if !( 0. <= e < 1. )
        throw(ArgumentError("The eccentricity must be within the interval 0 <= e < 1."))
    end

    # Check if the perigee is inside Earth.
    (a*(1.-e) > R0)
end

"""
### function n_J0(a::Real)

Return the orbit angular velocity neglecting the perturbations (two-body).

##### Args

* a: Semi-major axis [m].

##### Returns

* The unperturbed orbit angular velocity [rad/s].

"""

function n_J0(a::Real)
    # Check if the arguments are valid.
    if ( a < 0. )
        throw(ArgumentError("The semi-major axis must be greater than 0."))
    end

    sqrt(m0/Float64(a)^3)
end

"""
### function n_J2(a::Real, e::Real, i::Real)

Return the orbit angular velocity considering the perturbations up to J2 terms.

##### Args

* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].

##### Returns

* The perturbed orbit angular velocity [rad/s].

"""

function n_J2(a::Real, e::Real, i::Real)
    # Check if the perigee is inside Earth.
    if ( !is_orbit_valid(a,e) )
        throw(OrbitInvalidPerigee(a*(1-e)))
    end

    # Semi-lactum rectum.
    p = a*(1.0-e^2)

    # Unperturbed orbit period.
    n0 = n_J0(a)

    # Orbit period considering the perturbations (up to J2).
    n0 + 3.0*R0^2*J2/(4.0*p^2)*n0*(sqrt(1.0-e^2)*(3.0*cos(i)^2-1.0) +
                                   (5.0*cos(i)^2-1.0))
end

"""
### function sim_RAAN_J2(JD0::Real, a::Real, e::Real, i::Real, RAAN_0::Real, numDays::Integer)

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
    dOmega = dRAAN_J2(a, e, i)*24.0*3600.0

    # Simulate the RAAN for each day considering just the J2 perturbations.
    RAAN = mod(RAAN_0 + dOmega.*days,2*pi)

    [days RAAN]
end

"""
### function t_J0(a::Real)

Return the orbit period neglecting the perturbations (two-body).

##### Args

* a: Semi-major axis [m].

##### Returns

* The unperturbed orbit period [s].

"""

function t_J0(a::Real)
    2.0*pi/n_J0(a)
end

"""
### function t_J2(a::Real, e::Real, i::Real)

Return the orbit period considering the perturbations up to J2 terms.

##### Args

* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].

##### Returns

* The perturbed orbit period [s].

"""

function t_J2(a::Real, e::Real, i::Real)
    2.0*pi/n_J2(a, e, i)
end
