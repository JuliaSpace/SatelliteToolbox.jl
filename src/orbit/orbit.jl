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

export dRAAN_J2, dw_J2, n_J0, n_J2, t_J0, t_J2

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
