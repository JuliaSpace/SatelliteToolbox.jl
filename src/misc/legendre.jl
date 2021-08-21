# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions related to the associated Legendre functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Holmes, S. A. and W. E. Featherstone, 2002. A unified approach to the
#       Clenshaw summation and the recursive computation of very high degree and
#       order normalised associated Legendre functions. Journal of Geodesy,
#       76(5), pp. 279-299.
#
#       For more info.: http://mitgcm.org/~mlosch/geoidcookbook/node11.html
#
#   [2] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [3] Schmidt, A (1917). Erdmagnetismus, Enzykl. Math. Wiss., 6, pp. 265–396.
#
#   [4] Winch, D. E., Ivers, D. J., Turner, J. P. R., Stening R. J (2005).
#       Geomagnetism and Schmidt quasi-normalization. Geophysical Journal
#       International, 160(2), pp. 487-504.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export denormalize_legendre, legendre!, legendre

export legendre_fully_normalized!, legendre_fully_normalized
export legendre_schmidt_quasi_normalized!, legendre_schmidt_quasi_normalized
export legendre_conventional!, legendre_conventional

"""
    legendre!([N,] P::AbstractMatrix, ϕ::Number, ph_term::Bool = false, n_max::Integer = -1, m_max::Integer = -1)

Compute the associated Legendre function `P_n,m[cos(ϕ)]`. The maximum degree and
order that will be computed are given by the parameters `n_max` and `m_max`. If
they are negative, then the dimensions of matrix `P` will be used.

The result will be stored at matrix `P`.

The optional parameter `N` can be used to select the normalization. The
following values are valid:

- `Val(:full)`: Compute the fully normalized associated Legendre function (see
    [`legendre_fully_normalized!`](@ref)).
- `Val(:schmidt)`: Compute the Schmidt quasi-normalized associated Legendre
    function (see [`legendre_schmidt_quasi_normalized!`](@ref)).
- `Val(:conv)`: Compute the conventional associated Legendre function (see
    [`legendre_conventional!`](@ref)).

If `N` is omitted, then the full normalization will be used.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.
"""
function legendre!(
    P::AbstractMatrix,
    ϕ::Number,
    ph_term::Bool = false,
    n_max::Integer = -1,
    m_max::Integer = -1
)
    return legendre_fully_normalized!(P, float(ϕ), ph_term, n_max, m_max)
end

function legendre!(
    ::Val{:full},
    P::AbstractMatrix,
    ϕ::Number,
    ph_term::Bool = false,
    n_max::Integer = -1,
    m_max::Integer = -1
)
    return legendre_fully_normalized!(P, float(ϕ), ph_term, n_max, m_max)
end

function legendre!(
    ::Val{:schmidt},
    P::AbstractMatrix,
    ϕ::Number,
    ph_term::Bool = false,
    n_max::Integer = -1,
    m_max::Integer = -1
)
    return legendre_schmidt_quasi_normalized!(P, float(ϕ), ph_term, n_max, m_max)
end

function legendre!(
    ::Val{:conv},
    P::AbstractMatrix,
    ϕ::Number,
    ph_term::Bool = false,
    n_max::Integer = -1,
    m_max::Integer = -1
)
    return legendre_conventional!(P, float(ϕ), ph_term, n_max, m_max)
end

"""
    legendre([N,] ϕ::Number, n_max::Integer, m_max::Integer = -1, ph_term::Bool = false)

Compute the associated Legendre function `P_n,m[cos(ϕ)]`. The maximum degree
that will be computed is `n_max` and the maximum order is `m_max`. Notice that
if `m_max` is higher than `n_max` or negative, than it is set to `n_max`.

The optional parameter `N` can be used to select the normalization. The
following values are valid:

- `Val(:full)`: Compute the fully normalized associated Legendre function (see
    [`legendre_fully_normalized`](@ref)).
- `Val(:schmidt)`: Compute the Schmidt quasi-normalized associated Legendre
    function (see [`legendre_schmidt_quasi_normalized`](@ref)).
- `Val(:conv)`: Compute the conventional associated Legendre function (see
    [`legendre_conventional`](@ref)).

If `N` is omitted, then the full normalization will be used (`Val(:full)`).

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.
"""
function legendre(
    ϕ::Number,
    n_max::Integer,
    m_max::Integer = -1,
    ph_term::Bool = false
)
    return legendre_fully_normalized(float(ϕ), n_max, m_max, ph_term)
end

function legendre(
    ::Val{:full},
    ϕ::Number,
    n_max::Integer,
    m_max::Integer = -1,
    ph_term::Bool = false
)
    return legendre_fully_normalized(float(ϕ), n_max, m_max, ph_term)
end

function legendre(
    ::Val{:schmidt},
    ϕ::Number,
    n_max::Integer,
    m_max::Integer = -1,
    ph_term::Bool = false
)
    return legendre_schmidt_quasi_normalized(float(ϕ), n_max, m_max, ph_term)
end

function legendre(
    ::Val{:conv},
    ϕ::Number,
    n_max::Integer,
    m_max::Integer = -1,
    ph_term::Bool = false
)
    return legendre_conventional(float(ϕ), n_max, m_max, ph_term)
end

################################################################################
#                Fully Normalized Associated Legendre Functions
################################################################################

"""
    legendre_fully_normalized!(P::AbstractMatrix, ϕ::Number, ph_term::Bool = false, n_max::Integer = -1, m_max::Integer = -1)

Compute the fully normalized associated Legendre function `P_n,m[cos(ϕ)]`. The
maximum degree and order that will be computed are given by the parameters
`n_max` and `m_max`. If they are negative, then the dimensions of matrix `P`
will be used:

    maximum degree -> number of rows
    maximum order  -> number of columns

The result will be stored at matrix `P`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Remarks

This algorithm was based on **[1]**. Our definition of fully normalized associated
Legendre function can be seen in **[2, p. 546]**. The conversion is obtained by:

                 _                     -
                |  (n-m)! . k . (2n+1)  |      k = 1 if m  = 0
    K_n,m = sqrt| --------------------- |,     k = 2 if m != 0
                |         (n+m)!        |
                 -                     -
    _
    P_n,m = P_n,m * K_n,m,

          _
    where P_n,m is the fully normalized Legendre associated function.

# References

- **[1]** Holmes, S. A. and W. E. Featherstone, 2002. A unified approach to the
    Clenshaw summation and the recursive computation of very high degree and
    order normalised associated Legendre functions. Journal of Geodesy,
    76(5), pp. 279-299. For more info.:
    http://mitgcm.org/~mlosch/geoidcookbook/node11.html

- **[2]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
function legendre_fully_normalized!(
    P::AbstractMatrix,
    ϕ::Number,
    ph_term::Bool = false,
    n_max::Integer = -1,
    m_max::Integer = -1
)
    # Obtain the maximum degree and order that must be computed.
    n_max, m_max = _get_degree_and_order(P, n_max, m_max)

    # Auxiliary variables to improve code performance.
    c = cos(ϕ)
    s = sqrt(1 - c^2)

    s_fact = !ph_term ? +s : -s

    @inbounds for n in 0:n_max
        # Starting values.
        if n == 0
            P[0+1,0+1] = 1
            continue
        elseif n == 1
            P[1+1,0+1] = +sqrt(3) * c

            if m_max > 0
                P[1+1,1+1] = +sqrt(3) * s_fact
            end

            continue
        end

        aux_n = (2n - 1) * (2n + 1)

        for m in 0:n

            if n == m
                P[n+1,n+1] = s_fact * sqrt((2n + 1) / (2n)) * P[n-1+1,n-1+1]
            else
                aux_nm = (n - m) * (n + m)
                a_nm   = sqrt(aux_n / aux_nm)
                b_nm   = sqrt(((2n + 1) * (n + m - 1)*(n - m - 1)) / (aux_nm * (2n - 3)))

                # We assume that the matrix is not initialized. Hence, we must
                # not access elements on the upper triangle.
                if m != n - 1
                    P[n+1,m+1] = a_nm * c * P[n-1+1,m+1] - b_nm * P[n-2+1,m+1]
                else
                    P[n+1,m+1] = a_nm * c * P[n-1+1,m+1]
                end
            end

            # Check if the maximum desired order has been reached.
            m == m_max && break
        end
    end

    return nothing
end

"""
    legendre_fully_normalized(ϕ::T, n_max::Integer, m_max::Integer = -1, ph_term::Bool = false) where T<:AbstractFloat

Compute the fully normalized associated Legendre function `P_n,m[cos(ϕ)]`. The
maximum degree that will be computed is `n_max` and the maximum order is
`m_max`. Notice that if `m_max` is higher than `n_max` or negative, than it is
set to `n_max`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.

# Remarks

This algorithm was based on [1]. Our definition of fully normalized associated
Legendre function can be seen in [2, p. 546]. The conversion is obtained by:

                 _                     -
                |  (n-m)! . k . (2n+1)  |      k = 1 if m  = 0
    K_n,m = sqrt| --------------------- |,     k = 2 if m != 0
                |         (n+m)!        |
                 -                     -
    _
    P_n,m = P_n,m * K_n,m,

          _
    where P_n,m is the fully normalized Legendre associated function.

# References

- **[1]** Holmes, S. A. and W. E. Featherstone, 2002. A unified approach to the
    Clenshaw summation and the recursive computation of very high degree and
    order normalised associated Legendre functions. Journal of Geodesy,
    76(5), pp. 279-299. For more info.:
    http://mitgcm.org/~mlosch/geoidcookbook/node11.html

- **[2]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
function legendre_fully_normalized(
    ϕ::T,
    n_max::Integer,
    m_max::Integer = -1,
    ph_term::Bool = false
) where T<:AbstractFloat
    (n_max < 0) && throw(ArgumentError("n_max must be positive."))

    if ((m_max < 0) || (m_max > n_max))
        m_max = n_max
    end

    P = zeros(T, n_max + 1, m_max + 1)
    legendre_fully_normalized!(P, ϕ, ph_term)

    return P
end

################################################################################
#            Schmidt Quasi-Normalized Associated Legendre Functions
################################################################################

"""
    legendre_schmidt_quasi_normalized!(P::AbstractMatrix, ϕ::Number, ph_term::Bool = false, n_max::Integer = -1, m_max::Integer = -1)

Compute the Schmidt quasi-normalized associated Legendre function
`P_n,m[cos(ϕ)]` [3,4]. The maximum degree and order that will be computed are
given by the parameters `n_max` and `m_max`. If they are negative, then the
dimensions of matrix `P` will be used:

    maximum degree -> number of rows
    maximum order  -> number of columns

The result will be stored at matrix `P`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Remarks

This algorithm was based on **[1, 2]**. The conversion is obtained by:

                 _           -
                |     (n-m)!  |    k = 1 if m  = 0
    K_n,m = sqrt| k. -------- |,   k = 2 if m != 0
                |     (n+m)!  |
                 -           -

    =
    P_n,m = P_n,m * K_n,m,

          =
    where P_n,m is the quasi-normalized normalized Legendre associated function.

# References

- **[1]** Schmidt, A (1917). Erdmagnetismus, Enzykl. Math. Wiss., 6, pp.
    265–396.

- **[2]** Winch, D. E., Ivers, D. J., Turner, J. P. R., Stening R. J (2005).
    Geomagnetism and Schmidt quasi-normalization. Geophysical Journal
    International, 160(2), pp. 487-504.
"""
function legendre_schmidt_quasi_normalized!(
    P::AbstractMatrix,
    ϕ::Number,
    ph_term::Bool = false,
    n_max::Integer = -1,
    m_max::Integer = -1
)
    # Obtain the maximum degree and order that must be computed.
    n_max, m_max = _get_degree_and_order(P, n_max, m_max)

    # Auxiliary variables to improve code performance.
    c = cos(ϕ)
    s = sqrt(1 - c^2)

    s_fact = !ph_term ? +s : -s

    @inbounds for n in 0:n_max
        # Starting values.
        if n == 0
            P[0+1,0+1] = 1
            continue

        elseif n == 1
            P[1+1,0+1] = +c

            if m_max > 0
                P[1+1,1+1] = +s_fact
            end

            continue
        end

        aux_n = 2n - 1 # -> sqrt( (2n-1)*(2n-1) )

        for m in 0:n

            if m == n
                P[n+1,n+1] = s_fact * sqrt(aux_n / (2n))*P[n-1+1,n-1+1]
            else
                aux_nm = sqrt((n - m) * (n + m))
                a_nm   = aux_n / aux_nm
                b_nm   = sqrt((n + m - 1) * (n - m - 1)) / aux_nm

                # We assume that the matrix is not initialized. Hence, we must not
                # access elements on the upper triangle.
                if m != n - 1
                    P[n+1,m+1] = a_nm * c * P[n-1+1,m+1] - b_nm * P[n-2+1,m+1]
                else
                    P[n+1,m+1] = a_nm * c * P[n-1+1,m+1]
                end
            end

            # Check if the maximum desired order has been reached.
            m == m_max && break
        end
    end

    return nothing
end

"""
    legendre_schmidt_quasi_normalized(ϕ::T, n_max::Integer, m_max::Integer = -1, ph_term::Bool = false) where T<:AbstractFloat

Compute the Schmidt quasi-normalized associated Legendre function
`P_n,m[cos(ϕ)]`. The maximum degree that will be computed is `n_max` and the
maximum order is `m_max`. Notice that if `m_max` is higher than `n_max` or
negative, than it is set to `n_max`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.

# Remarks

This algorithm was based on **[1, 2]**. The conversion is obtained by:

                 _           -
                |     (n-m)!  |    k = 1 if m  = 0
    K_n,m = sqrt| k. -------- |,   k = 2 if m != 0
                |     (n+m)!  |
                 -           -

    =
    P_n,m = P_n,m * K_n,m,

          =
    where P_n,m is the quasi-normalized normalized Legendre associated function.

# References

- **[1]** Schmidt, A (1917). Erdmagnetismus, Enzykl. Math. Wiss., 6, pp.
    265–396.

- **[2]** Winch, D. E., Ivers, D. J., Turner, J. P. R., Stening R. J (2005).
    Geomagnetism and Schmidt quasi-normalization. Geophysical Journal
    International, 160(2), pp. 487-504.
"""
function legendre_schmidt_quasi_normalized(
    ϕ::T,
    n_max::Integer,
    m_max::Integer = -1,
    ph_term::Bool = false
) where T<:AbstractFloat
    (n_max < 0) && throw(ArgumentError("n_max must be positive."))

    if ((m_max < 0) || (m_max > n_max))
        m_max = n_max
    end

    P = zeros(T, n_max + 1, m_max + 1)
    legendre_schmidt_quasi_normalized!(P, ϕ, ph_term)

    return P
end

################################################################################
#                  Conventional Associated Legendre Function
################################################################################

"""
    legendre_conventional!(P::AbstractMatrix, ϕ::Number, ph_term::Bool = false, n_max::Integer = -1, m_max::Integer = -1)

Compute the conventional associated Legendre function `P_n,m[cos(ϕ)]`. The
maximum degree and order that will be computed are given by the parameters
`n_max` and `m_max`. If they are negative, then the dimensions of matrix `P`
will be used:

    maximum degree -> number of rows
    maximum order  -> number of columns

The result will be stored at matrix `P`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.
"""
function legendre_conventional!(
    P::AbstractMatrix,
    ϕ::Number,
    ph_term::Bool = false,
    n_max::Integer = -1,
    m_max::Integer = -1
)
    # Obtain the maximum degree and order that must be computed.
    n_max, m_max = _get_degree_and_order(P, n_max, m_max)

    # Auxiliary variables to improve code performance.
    c = cos(ϕ)
    s = sqrt(1 - c^2)

    s_fact = !ph_term ? +s : -s

    @inbounds for n in 0:n_max
        # Starting values.
        if n == 0
            P[0+1,0+1] = 1
            continue
        elseif n == 1
            P[1+1,0+1] = +c

            if m_max > 0
                P[1+1,1+1] = +s_fact
            end

            continue
        end

        aux_n = 2n - 1 # -> sqrt( (2n-1)*(2n-1) )

        for m in 0:n

            if n == m
                P[n+1,n+1] = s_fact * aux_n * P[n-1+1,n-1+1]
            else
                aux_nm = n-m # -> sqrt( (n-m)*(n-m) )
                a_nm   = aux_n / aux_nm
                b_nm   = (n + m - 1) / aux_nm  # -> sqrt( (n+m-1)*(n+m-1) ) / aux_nm

                # We assume that the matrix is not initialized. Hence, we must
                # not access elements on the upper triangle.
                if m != n-1
                    P[n+1,m+1] = a_nm * c * P[n-1+1,m+1] - b_nm * P[n-2+1,m+1]
                else
                    P[n+1,m+1] = a_nm * c *P[n-1+1,m+1]
                end
            end

            # Check if the maximum desired order has been reached.
            m == m_max && break
        end
    end

    return nothing
end

"""
    legendre_conventional(ϕ::T, n_max::Integer, m_max::Integer = -1, ph_term::Bool = false) where T<:AbstractFloat

Compute the conventional associated Legendre function `P_n,m[cos(ϕ)]`. The
maximum degree that will be computed is `n_max` and the maximum order is
`m_max`. Notice that if `m_max` is higher than `n_max` or negative, than it is
set to `n_max`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.
"""
function legendre_conventional(
    ϕ::T,
    n_max::Integer,
    m_max::Integer = -1,
    ph_term::Bool = false
) where T<:AbstractFloat
    (n_max < 0) && throw(ArgumentError("n_max must be positive."))

    if ((m_max < 0) || (m_max > n_max))
        m_max = n_max
    end

    P = zeros(T, n_max + 1, m_max + 1)
    legendre_conventional!(P, ϕ, ph_term)

    return P
end

################################################################################
#                                   Private
################################################################################

"""
    _get_degree_and_order(P, n_max, m_max)

Return the maximum degree and order to compute the Legendre associated functions
given the matrix `P` and the configuration values `n_max` and `m_max`.
"""
@inline function _get_degree_and_order(P, n_max, m_max)
    # Get the size of the matrix.
    rows, cols = size(P)

    # If the order or degree is less than 0, then the user wants to use all the
    # available memory.
    if n_max < 0
        n_max = rows - 1
    end

    if m_max < 0
        m_max = (cols <= rows) ? cols - 1 : n_max
    end

    # Make sure that the degree and order fits the matrix.
    if n_max > rows - 1
        n_max = rows - 1
    end

    if ((m_max > cols - 1) || (m_max > n_max))
        m_max = min(cols - 1, n_max)
    end

    return n_max, m_max
end
