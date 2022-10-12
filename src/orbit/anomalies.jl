# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions to convert anomalies related to the orbit.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export M_to_E, M_to_f
export E_to_f, E_to_M
export f_to_E, f_to_M

################################################################################
#                              From Mean Anomaly
################################################################################

"""
    M_to_E(e::T1, M::T2; kwargs...) where {T1, T2}

Compute the eccentric anomaly (0,2π) \\[rad] given the eccentricity `e` and the
mean anomaly `M` [rad].

This function uses the Newton-Raphson algorithm to solve the Kepler's equation.

# Keywords

- `tol::Union{Nothing, Number}`: Tolerance to accept the solution from
    Newton-Raphson algorithm. If `tol` is `nothing`, then it will be
    `eps(T)`, where `T` is a floating-point type obtained from the promotion of
    `T1` and `T2`. (**Default** = `nothing`)
- `max_iterations::Number`: Maximum number of iterations allowed for the
    Newton-Raphson algorithm. If it is lower than 1, then it is set to 10.
    (**Default** = 10)
"""
@inline function M_to_E(
    e::T1,
    M::T2;
    max_iterations::Integer = 10,
    tol::Union{Nothing, Number} = nothing
) where {T1, T2}
    T = float(promote_type(T1, T2))

    # Compute the eccentric anomaly using the Newton-Raphson method.
    # ==============================================================

    # Make sure that M is in the interval [0,2π].
    M = mod(M, T(2π))

    # Initial guess.
    #
    # See [1, p. 75].
    E = (M > π) ? M - e : M + e

    sin_E, cos_E = sincos(E)

    # Check the tolerance.
    δ = isnothing(tol) ? eps(T) : T(tol)

    # Check the maximum number of iterations.
    if max_iterations < 1
        max_iterations = 10
    end

    # Newton-Raphson iterations.
    for i in 1:max_iterations
        abs(E - e * sin_E - M) ≤ δ && break
        E = E - (E - e * sin_E - M) / (1 - e * cos_E)
        sin_E, cos_E = sincos(E)
    end

    # Return the eccentric anomaly in the interval [0, 2π].
    return mod(E, T(2π))
end

"""
    M_to_f(e::T1, M::T2; kwargs...) where {T1, T2}

Compute the true anomaly (0,2π) \\[rad] given the eccentricity `e` and the mean
anomaly `M` [rad].

This function uses the Newton-Raphson algorithm to solve the Kepler's equation.

# Keywords

- `tol::Union{Nothing, Number}`: Tolerance to accept the solution from
    Newton-Raphson algorithm. If `tol` is `nothing`, then it will be
    `eps(T)`, where `T` is a floating-point type obtained from the promotion of
    `T1` and `T2`. (**Default** = `nothing`)
- `max_iterations::Number`: Maximum number of iterations allowed for the
    Newton-Raphson algorithm. If it is lower than 1, then it is set to 10.
    (**Default** = 10)
"""
@inline function M_to_f(e::T1, M::T2; kwargs...) where {T1, T2}
    # Compute the eccentric anomaly.
    E = M_to_E(e, M; kwargs...)

    # Compute the true anomaly in the interval [0,2π].
    return E_to_f(e, E)
end

################################################################################
#                            From Eccentric Anomaly
################################################################################

"""
    E_to_f(e::T1, E::T2) where {T1, T2}

Compute the true anomaly (0,2π) \\[rad] given the eccentricity `e` and the
eccentric anomaly `E` [rad].
"""
@inline function E_to_f(e::T1, E::T2) where {T1, T2}
    T = float(promote_type(T1, T2))

    sin_Eo2, cos_Eo2 = sincos(E / 2)

    # Compute the true anomaly in the interval [0, 2*π].
    return mod(2atan(sqrt(1 + e) * sin_Eo2, sqrt(1 - e) * cos_Eo2), T(2π))
end

"""
    E_to_M(e::T1, E::T2) where {T1, T2}

Compute the mean anomaly (0,2π) \\[rad] given the eccentricity `e` and the
eccentric anomaly `E` [rad].
"""
@inline function E_to_M(e::T1, E::T2) where {T1, T2}
    T = float(promote_type(T1, T2))
    return mod(E - e * sin(E), T(2π))
end

################################################################################
#                              From True Anomaly
################################################################################

"""
    f_to_E(e::T1, f::T2) where {T1, T2}

Compute the eccentric anomaly (0,2π) \\[rad] given the eccentricity `e` and
the true anomaly `f` [rad].
"""
@inline function f_to_E(e::T1, f::T2) where {T1, T2}
    T = float(promote_type(T1, T2))
    sin_fo2, cos_fo2 = sincos(f / 2)
    return mod(2atan(sqrt(1 - e) * sin_fo2, sqrt(1 + e) * cos_fo2), T(2π) )
end

"""
    f_to_M(e::T1, f::T2) where {T1, T2}

Compute the mean anomaly (0,2π) \\[rad] given the eccentricity `e` and the
true anomaly `f` [rad].
"""
@inline function f_to_M(e::T1, f::T2) where {T1, T2}
    # Compute the eccentric anomaly.
    E = f_to_E(e, f)

    # Compute the true anomaly in the interval [0, 2π].
    return E_to_M(e, E)
end
