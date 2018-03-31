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
#   Functions to convert anomalies related to the orbit.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-03-27: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export M_to_E, M_to_f
export E_to_f, E_to_M
export f_to_E, f_to_M

################################################################################
#                              From Mean Anomaly
################################################################################

"""
### function M_to_E(e::Number, M::Number, tol::Number = 1e-10)

Compute the eccentric anomaly given the eccentricity `e` and the mean anomaly
`M`. This function uses the Newton-Raphson algorithm and the tolerance to accept
the solution is `tol`.

##### Args

* e: Eccentricity.
* M: Mean anomaly [rad].
* tol: (OPTIONAL) Tolerance for the Newton-Raphson method (**DEFAULT**: 1e-10).

##### Returns

The eccentric anomaly in the interval [0,2π].

"""

function M_to_E(e::Number, M::Number, tol::Number = 1e-10)
    # Compute the eccentric anomaly using the Newton-Raphson method.
    # ==============================================================

    # Make sure that M is in the interval [0,2π].
    M = mod(M,2*pi)

    # Initial guess.
    #
    # See [1, p. 75].
    E = (M > pi) ? M - e : M + e

    sin_E = sin(E)
    cos_E = cos(E)

    # Newton-Raphson iterations.
    while ( abs(E - e*sin_E - M) > tol )
        E = E - (E - e*sin_E - M)/(1-e*cos_E)

        sin_E = sin(E)
        cos_E = cos(E)
    end

    # Return the eccentric anomaly in the interval [0, 2*π].
    mod(E, 2*pi)
end

"""
### function M_to_f(e::Number, M::Number, tol::Number = 1e-10)

Compute the true anomaly given the eccentricity `e` and the mean anomaly `M`.
This function uses the Newton-Raphson algorithm and the tolerance to accept the
solution is `tol`.

##### Args

* e: Eccentricity.
* M: Mean anomaly [rad].
* tol: (OPTIONAL) Tolerance for the Newton-Raphson method (**DEFAULT**: 1e-10).

##### Returns

The true anomaly in the interval [0,2π].

"""

function M_to_f(e::Number, M::Number, tol::Number = 1e-10)
    # Compute the eccentric anomaly.
    E = M_to_E(e, M, tol)

    # Compute the true anomaly in the interval [0,2π].
    E_to_f(e,E)
end

################################################################################
#                            From Eccentric Anomaly
################################################################################

"""
### function E_to_f(e::Number, E::Number)

Compute the true anomaly given the eccentricity `e` and the eccentric anomaly
`E`.

##### Args

* e: Eccentricity.
* E: Eccentric anomaly [rad].

##### Returns

The true anomaly in the interval [0, 2π].

"""

function E_to_f(e::Number, E::Number)
    # Compute the true anomaly in the interval [0, 2*π].
    mod(2*atan2( sqrt(1+e)*sin(E/2), sqrt(1-e)*cos(E/2) ), 2*pi)
end

"""
### function E_to_M(e::Number, E::Number)

Compute the mean anomaly given the eccentricity `e` and the eccentric anomaly
`E`.

##### Args

* e: Eccentricity.
* E: Eccentric anomaly [rad].

##### Returns

The mean anomaly in the interval [0, 2π].

"""

function E_to_M(e::Number, E::Number)
    mod(E - e*sin(E), 2*pi)
end

################################################################################
#                              From True Anomaly
################################################################################

"""
### function f_to_E(e::Number,f::Number)

Compute the eccentric anomaly given the eccentricity `e` and the true anomaly
`f`.

##### Args

* e: Eccentricity.
* f: True anomaly [rad].

##### Returns

The eccentric anomaly in the interval [0,2π].

"""

function f_to_E(e::Number, f::Number)
    mod(2*atan2( sqrt(1-e)*sin(f/2), sqrt(1+e)*cos(f/2) ), 2*pi)
end

"""
### function f_to_E(orb::Orbit)

Compute the eccentric anomaly given the orbit `orb`.

##### Args

* orb: Orbit (see `Orbit`).

##### Returns

The eccentric anomaly in the interval [0,2π].

"""

function f_to_E(orb::Orbit)
    f_to_E(orb.e, orb.f)
end

"""
### function f_to_M(e::Number, f::Number)

Compute the mean anomaly given the eccentricity `e` and the true anomaly `f`.

##### Args

* e: Eccentricity.
* f: True anomaly [rad].

##### Returns

The mean anomaly in the interval [0,2π].

"""

function f_to_M(e::Number, f::Number)
    # Compute the eccentric anomaly.
    E = f_to_E(e, f)

    # Compute the true anomaly in the interval [0, 2π].
    E_to_M(e,E)
end

"""
### function f_to_M(orb::Orbit)

Compute the mean anomaly given the orbit `orb`.

##### Args

* orb: Orbit (see `Orbit`).

##### Returns

The mean anomaly in the interval [0,2π].

"""

function f_to_M(orb::Orbit)
    f_to_M(orb.e, orb.f)
end
