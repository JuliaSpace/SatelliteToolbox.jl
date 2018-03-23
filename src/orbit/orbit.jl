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
# 2018-03-23: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#    Change API. The functions now receive a symbol that specifies which kind of
#    perturbations must be considered to compute the values.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export angvel, dArgPer, dRAAN, period

"""
### function angvel(a::Number, e::Number, i::Number, pert::Symbol = :J2)

Compute the angular velocity of an object in an orbit with semi-major axis `a`,
eccentricity `e`, and inclination `i`, using the perturbation terms specified by
the symbol `pert`.

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.

##### Args

* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].
* pert: (OPTIONAL) Symbol that defines the perturbation (**DEFAULT** = `:J2`).

##### Returns

The angular velocity of an object in the specified orbit [rad/s].

"""

function angvel(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    # Check if the perigee is inside Earth.
    if ( !is_orbit_valid(a,e) )
        throw(OrbitInvalidPerigee(a*(1-e)))
    end

    # Unperturbed orbit period.
    n0 = sqrt(m0/Float64(a)^3)

    # Perturbation computed using a Keplerian orbit.
    if pert == :J0
        return n0
    # Perturbation computed using perturbations terms up to J2.
    elseif pert == :J2
        # Semi-lactum rectum.
        p = a*(1.0-e^2)

        # Orbit period considering the perturbations (up to J2).
        return n0 + 3.0*R0^2*J2/(4.0*p^2)*n0*(sqrt(1.0-e^2)*(3.0*cos(i)^2-1.0) +
                                              (5.0*cos(i)^2-1.0))
    else
        throw(ArgumentError("The perturbation parameter $pert is not defined."))
    end

end

"""
### function dArgPer(a::Number, e::Number, i::Number, pert::Symbol = :J2)

Compute the time-derivative of the argument of perigee of a orbit with
semi-major axis `a`, eccentricity `e`, and inclination `i`, using the
perturbation terms specified by the symbol `pert`.

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.

##### Args

* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].
* pert: (OPTIONAL) Symbol that defines the perturbation (**DEFAULT** = `:J2`).

##### Returns

The perturbation of the argument of perigee [rad/s].

"""

function dArgPer(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    # Check if the perigee is inside Earth.
    if ( !is_orbit_valid(a,e) )
        throw(OrbitInvalidPerigee(a*(1.0-e)))
    end

    # Perturbation computed using a Keplerian orbit.
    if pert == :J0
        return 0.0
    # Perturbation computed using perturbations terms up to J2.
    elseif pert == :J2
        # Semi-lactum rectum.
        p = a*(1.0-e^2)

        # Unperturbed orbit period.
        n0 = period(a, e, i, :J0)

        # Perturbation of the argument of perigee.
        return 3.0*R0^2*J2/(4.0*p^2)*n0*(5.0*cos(i)^2-1.0)
    else
        throw(ArgumentError("The perturbation parameter $pert is not defined."))
    end

end

"""
### function dRAAN(a::Number, e::Number, i::Number, pert::Symbol = :J2)

Compute the time-derivative of the right ascencion of the ascending node of a
orbit with semi-major axis `a`, eccentricity `e`, and inclination `i`, using the
perturbation terms specified by the symbol `pert`.

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.

##### Args

* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].
* pert: Symbol that defines the perturbation (DEFAULT = `:J2`).

##### Returns

The time derivative of the RAAN [rad/s].

"""

function dRAAN(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    # Check if the perigee is inside Earth.
    if ( !is_orbit_valid(a,e) )
        throw(OrbitInvalidPerigee(a*(1.0-e)))
    end

    # Perturbation computed using a Keplerian orbit.
    if pert == :J0
        return 0.0
    # Perturbation computed using perturbations terms up to J2.
    elseif pert == :J2
        # Semi-lactum rectum.
        p = a*(1.0-e^2)

        # Unperturbed orbit period.
        n0 = period(a, e, i, :J0)

        # Perturbation of the right ascension of the ascending node.
        return -3.0/2.0*R0^2/(p^2)*n0*J2*cos(i)
    else
        throw(ArgumentError("The perturbation parameter $pert is not defined."))
    end
end


"""
### function is_orbit_valid(a::Number, e::Number)

Verify if the orbit is valid.

##### Args

* a: Semi-major axis [m].
* e: Eccentricity.

##### Returns

* **TRUE**: The orbit is valid.
* **FALSE**: The orbit is invalid.

"""

function is_orbit_valid(a::Number, e::Number)
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
### function period(a::Number, e::Number, i::Number, pert::Symbol = :J2)

Compute the period of an object in an orbit with semi-major axis `a`,
eccentricity `e`, and inclination `i`, using the perturbation terms specified by
the symbol `pert`.

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.

##### Args

* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].
* pert: (OPTIONAL) Symbol that defines the perturbation (**DEFAULT** = `:J2`).

##### Returns

The orbit period [s].

"""

function period(a::Number, e::Number, i::Number, pert::Symbol)
    n = angvel(a, e, i, pert)
    2.0*pi/n
end
