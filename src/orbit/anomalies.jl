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
    M_to_E(e::Number, M::Number, tol::Number = 1e-10)

Compute the eccentric anomaly (0,2π) \\[rad] given the eccentricity `e` and the
mean anomaly `M` [rad]. This function uses the Newton-Raphson algorithm and the
tolerance to accept the solution is `tol`.

"""
@inline function M_to_E(e::Number, M::Number, tol::Number = 1e-10)
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
    return mod(E, 2π)
end

"""
    M_to_f(e::Number, M::Number, tol::Number = 1e-10)

Compute the true anomaly (0,2π) \\[rad] given the eccentricity `e` and the mean
anomaly `M` [rad]. This function uses the Newton-Raphson algorithm and the
tolerance to accept the solution is `tol`.

"""
@inline function M_to_f(e::Number, M::Number, tol::Number = 1e-10)
    # Compute the eccentric anomaly.
    E = M_to_E(e, M, tol)

    # Compute the true anomaly in the interval [0,2π].
    return E_to_f(e,E)
end

################################################################################
#                            From Eccentric Anomaly
################################################################################

"""
    E_to_f(e::Number, E::Number)

Compute the true anomaly (0,2π) \\[rad] given the eccentricity `e` and the
eccentric anomaly `E` [rad].

"""
@inline function E_to_f(e::Number, E::Number)
    sin_Eo2, cos_Eo2 = sincos(E/2)

    # Compute the true anomaly in the interval [0, 2*π].
    return mod( 2atan(sqrt(1+e)*sin_Eo2, sqrt(1-e)*cos_Eo2) , 2π )
end

"""
    E_to_M(e::Number, E::Number)

Compute the mean anomaly (0,2π) \\[rad] given the eccentricity `e` and the
eccentric anomaly `E` [rad].

"""
@inline E_to_M(e::Number, E::Number) = mod(E - e*sin(E), 2π)

################################################################################
#                              From True Anomaly
################################################################################

"""
    f_to_E(e::Number,f::Number)

Compute the eccentric anomaly (0,2π) \\[rad] given the eccentricity `e` and
the true anomaly `f` [rad].

"""
@inline function f_to_E(e::Number, f::Number)
    sin_fo2, cos_fo2 = sincos(f/2)

    return mod( 2atan(sqrt(1-e)*sin_fo2, sqrt(1+e)*cos_fo2), 2π )
end

"""
    f_to_M(e::Number, f::Number)

Compute the mean anomaly (0,2π) \\[rad] given the eccentricity `e` and the
true anomaly `f` [rad].

"""
@inline function f_to_M(e::Number, f::Number)
    # Compute the eccentric anomaly.
    E = f_to_E(e, f)

    # Compute the true anomaly in the interval [0, 2π].
    return E_to_M(e,E)
end
