#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export M_to_E, M_to_f
export E_to_f, E_to_M
export f_to_E, f_to_M

################################################################################
#                              From Mean Anomaly
################################################################################

"""
    function M_to_E(e::Number, M::Number, tol::Number = 1e-10)

Compute the eccentric anomaly (0,2π) \\[rad] given the eccentricity `e` and the
mean anomaly `M` [rad]. This function uses the Newton-Raphson algorithm and the
tolerance to accept the solution is `tol`.

"""
function M_to_E(e::Number, M::Number, tol::Number = 1e-10)
    # Compute the eccentric anomaly using the Newton-Raphson method.
    # ==============================================================

    # Make sure that M is in the interval [0,2π].
    M = mod(M,2π)

    # Initial guess.
    #
    # See [1, p. 75].
    E = (M > π) ? M - e : M + e

    sin_E, cos_E = sincos(E)

    # Newton-Raphson iterations.
    while ( abs(E - e*sin_E - M) > tol )
        E = E - (E - e*sin_E - M)/(1-e*cos_E)

        sin_E, cos_E = sincos(E)
    end

    # Return the eccentric anomaly in the interval [0, 2π].
    mod(E, 2π)
end

"""
    function M_to_f(e::Number, M::Number, tol::Number = 1e-10)

Compute the true anomaly (0,2π) \\[rad] given the eccentricity `e` and the mean
anomaly `M` [rad]. This function uses the Newton-Raphson algorithm and the
tolerance to accept the solution is `tol`.

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
    function E_to_f(e::Number, E::Number)

Compute the true anomaly (0,2π) \\[rad] given the eccentricity `e` and the
eccentric anomaly `E` [rad].

"""
function E_to_f(e::Number, E::Number)
    sin_Eo2, cos_Eo2 = sincos(E/2)

    # Compute the true anomaly in the interval [0, 2*π].
    mod( 2atan(sqrt(1+e)*sin_Eo2, sqrt(1-e)*cos_Eo2) , 2π )
end

"""
    function E_to_M(e::Number, E::Number)

Compute the mean anomaly (0,2π) \\[rad] given the eccentricity `e` and the
eccentric anomaly `E` [rad].

"""
function E_to_M(e::Number, E::Number)
    mod(E - e*sin(E), 2*pi)
end

################################################################################
#                              From True Anomaly
################################################################################

"""
    function f_to_E(e::Number,f::Number)

Compute the eccentric anomaly (0,2π) \\[rad] given the eccentricity `e` and
the true anomaly `f` [rad].

"""
function f_to_E(e::Number, f::Number)
    mod(2*atan( sqrt(1-e)*sin(f/2), sqrt(1+e)*cos(f/2) ), 2*pi)
end

"""
    function f_to_E(orb::Orbit)

Compute the eccentric anomaly (0,2π) \\[rad] given the orbit `orb` (see
`Orbit`).

"""
function f_to_E(orb::Orbit)
    f_to_E(orb.e, orb.f)
end

"""
    function f_to_M(e::Number, f::Number)

Compute the mean anomaly (0,2π) \\[rad] given the eccentricity `e` and the
true anomaly `f` [rad].

"""
function f_to_M(e::Number, f::Number)
    # Compute the eccentric anomaly.
    E = f_to_E(e, f)

    # Compute the true anomaly in the interval [0, 2π].
    E_to_M(e,E)
end

"""
    function f_to_M(orb::Orbit)

Compute the mean anomaly (0,2π) \\[rad] given the orbit `orb` (see `Orbit`).

"""
function f_to_M(orb::Orbit)
    f_to_M(orb.e, orb.f)
end
