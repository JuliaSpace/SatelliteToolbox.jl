#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions related to the associated Legendre functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export denormalize_legendre, legendre!, legendre

export legendre_fully_normalized!, legendre_fully_normalized
export legendre_schmidt_quasi_normalized!, legendre_schmidt_quasi_normalized
export legendre_conventional!, legendre_conventional

"""
    function legendre!([N,] P::AbstractMatrix, ϕ::Number, ph_term::Bool = false)

Compute the associated Legendre function `P_n,m[cos(ϕ)]`. The maximum degree
and order that will be computed are given by the dimensions of matrix `P`.

The result will be stored at matrix `P`.

The optional parameter `N` can be used to select the normalization. The
following values are valid:

* `Val{:full}`: Compute the fully normalized associated Legendre function (see
  `legendre_fully_normalized!`).
* `Val{:schmidt}`: Compute the Schmidt quasi-normalized associated Legendre
  function (see `legendre_schmidt_quasi_normalized!`).
* `Val{:conv}`: Compute the conventional associated Legendre function (see
  `legendre_conventional!`).

If `N` is omitted, then the full normalization will be used.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

"""
legendre!(P::AbstractMatrix, ϕ::Number, ph_term::Bool = false) =
    legendre_fully_normalized!(P, float(ϕ), ph_term)

legendre!(::Type{Val{:full}}, P::AbstractMatrix, ϕ::Number, ph_term::Bool = false) =
    legendre_fully_normalized!(P, float(ϕ), ph_term)

legendre!(::Type{Val{:schmidt}}, P::AbstractMatrix, ϕ::Number, ph_term::Bool = false) =
    legendre_schmidt_quasi_normalized!(P, float(ϕ), ph_term)

legendre!(::Type{Val{:conv}}, P::AbstractMatrix, ϕ::Number, ph_term::Bool = false) =
    legendre_conventional!(P, float(ϕ), ph_term)

"""
    function legendre([N,] ϕ::Number, n_max::Number, m_max::Number = -1, ph_term::Bool = false)

Compute the associated Legendre function `P_n,m[cos(ϕ)]`. The maximum degree
that will be computed is `n_max` and the maximum order is `m_max`. Notice that
if `m_max` is higher than `n_max` or negative, than it is set to `n_max`.

The optional parameter `N` can be used to select the normalization. The
following values are valid:

* `Val{:full}`: Compute the fully normalized associated Legendre function (see
  `legendre_fully_normalized`).
* `Val{:schmidt}`: Compute the Schmidt quasi-normalized associated Legendre
  function (see `legendre_schmidt_quasi_normalized`).
* `Val{:conv}`: Compute the conventional associated Legendre function (see
  `legendre_conventional`).

If `N` is omitted, then the full normalization will be used (`Val{:full}`).

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.

"""
legendre(ϕ::Number, n_max::Number, m_max::Number = -1, ph_term::Bool = false) =
    legendre_fully_normalized(float(ϕ), n_max, m_max, ph_term)

legendre(::Type{Val{:full}}, ϕ::Number, n_max::Number, m_max::Number = -1,
         ph_term::Bool = false) =
    legendre_fully_normalized(float(ϕ), n_max, m_max, ph_term)

legendre(::Type{Val{:schmidt}}, ϕ::Number, n_max::Number, m_max::Number = -1,
         ph_term::Bool = false) =
    legendre_schmidt_quasi_normalized(float(ϕ), n_max, m_max, ph_term)

legendre(::Type{Val{:conv}}, ϕ::Number, n_max::Number, m_max::Number = -1,
         ph_term::Bool = false) =
    legendre_conventional(float(ϕ), n_max, m_max, ph_term)

################################################################################
#                Fully Normalized Associated Legendre Functions
################################################################################

"""
    function legendre_fully_normalized!(P::AbstractMatrix, ϕ::Number, ph_term::Bool = false)

Compute the fully normalized associated Legendre function `P_n,m[cos(ϕ)]`.
The maximum degree and order that will be computed are given by the dimensions
of matrix `P`:

    maximum degree -> number of rows
    maximum order  -> number of columns

The result will be stored at matrix `P`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

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

"""
function legendre_fully_normalized!(P::AbstractMatrix, ϕ::Number,
                                    ph_term::Bool = false)
    (rows, cols) = size(P)

    # Obtain the maximum degree and order that must be computed.
    n_max = rows - 1
    m_max = cols <= rows ? cols - 1 : n_max

    # Auxiliary variables to improve code performance.
    c = cos(ϕ)
    s = sqrt(1-c^2)

    s_fact = !ph_term ? +s : -s

    @inbounds for n = 0:n_max
        # Starting values.
        if n == 0
            P[0+1,0+1] = 1
            continue
        elseif n == 1
            P[1+1,0+1] = +sqrt(3)*c

            if cols > 1
                P[1+1,1+1] = +sqrt(3)*s_fact
            end

            continue
        end

        aux_n = (2n-1)*(2n+1)

        for m = 0:n

            if n == m
                P[n+1,n+1] = s_fact*sqrt( (2n+1)/(2n) )*P[n-1+1,n-1+1]
                continue
            end

            aux_nm = (n-m)*(n+m)
            a_nm   = sqrt( aux_n / aux_nm )
            b_nm   = sqrt( ( (2n+1)*(n+m-1)*(n-m-1) ) / ( aux_nm*(2n-3) ) )

            # We assume that the matrix is not initialized. Hence, we must not
            # access elements on the upper triangle.
            if m != n-1
                P[n+1,m+1] = a_nm*c*P[n-1+1,m+1] - b_nm*P[n-2+1,m+1]
            else
                P[n+1,m+1] = a_nm*c*P[n-1+1,m+1]
            end

            # Check if the maximum desired order has been reached.
            m == m_max && break
        end
    end

    nothing
end

"""
    function legendre_fully_normalized(ϕ::T, n_max::Number, m_max::Number = -1, ph_term::Bool = false) where T<:AbstractFloat

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

"""
function legendre_fully_normalized(ϕ::T, n_max::Number, m_max::Number = -1,
                                   ph_term::Bool = false) where T<:AbstractFloat

    (n_max < 0) && throw(ArgumentError("n_max must be positive."))

    ( (m_max < 0) || (m_max > n_max) ) && (m_max = n_max)

    P = zeros(T, n_max+1, m_max+1)
    legendre_fully_normalized!(P, ϕ, ph_term)
    return P
end

################################################################################
#            Schmidt Quasi-Normalized Associated Legendre Functions
################################################################################

"""
    function legendre_schmidt_quasi_normalized!(P::AbstractMatrix, ϕ::Number, ph_term::Bool = false)

Compute the Schmidt quasi-normalized associated Legendre function
`P_n,m[cos(ϕ)]` [3,4]. The maximum degree and order that will be computed are
given by the dimensions of matrix `P`:

    maximum degree -> number of rows
    maximum order  -> number of columns

The result will be stored at matrix `P`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Remarks

This algorithm was based on [3,4]. The conversion is obtained by:

                 _           -
                |     (n-m)!  |    k = 1 if m  = 0
    K_n,m = sqrt| k. -------- |,   k = 2 if m != 0
                |     (n+m)!  |
                 -           -

    =
    P_n,m = P_n,m * K_n,m,

          =
    where P_n,m is the quasi-normalized normalized Legendre associated function.

"""
function legendre_schmidt_quasi_normalized!(P::AbstractMatrix, ϕ::Number,
                                            ph_term::Bool = false)
    (rows, cols) = size(P)

    # Obtain the maximum degree and order that must be computed.
    n_max = rows - 1
    m_max = cols <= rows ? cols - 1 : n_max

    # Auxiliary variables to improve code performance.
    c = cos(ϕ)
    s = sqrt(1-c^2)

    s_fact = !ph_term ? +s : -s

    @inbounds for n = 0:n_max
        # Starting values.
        if n == 0
            P[0+1,0+1] = 1
            continue

        elseif n == 1
            P[1+1,0+1] = +c

            if cols > 1
                P[1+1,1+1] = +s_fact
            end

            continue
        end

        aux_n = 2n-1 # -> sqrt( (2n-1)*(2n-1) )

        for m = 0:n

            if m == n
                P[n+1,n+1] = s_fact*sqrt( aux_n/(2n) )*P[n-1+1,n-1+1]
                continue
            end

            aux_nm = sqrt( (n-m)*(n+m) )
            a_nm   = aux_n / aux_nm
            b_nm   = sqrt( (n+m-1)*(n-m-1) ) / aux_nm

            # We assume that the matrix is not initialized. Hence, we must not
            # access elements on the upper triangle.
            if m != n-1
                P[n+1,m+1] = a_nm*c*P[n-1+1,m+1] - b_nm*P[n-2+1,m+1]
            else
                P[n+1,m+1] = a_nm*c*P[n-1+1,m+1]
            end

            # Check if the maximum desired order has been reached.
            m == m_max && break
        end
    end

    nothing
end

"""
    function legendre_schmidt_quasi_normalized(ϕ::T, n_max::Number, m_max::Number = -1, ph_term::Bool = false) where T<:AbstractFloat

Compute the Schmidt quasi-normalized associated Legendre function
`P_n,m[cos(ϕ)]`. The maximum degree that will be computed is `n_max` and the
maximum order is `m_max`. Notice that if `m_max` is higher than `n_max` or
negative, than it is set to `n_max`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.

# Remarks

This algorithm was based on [3,4]. The conversion is obtained by:

                 _           -
                |     (n-m)!  |    k = 1 if m  = 0
    K_n,m = sqrt| k. -------- |,   k = 2 if m != 0
                |     (n+m)!  |
                 -           -

    =
    P_n,m = P_n,m * K_n,m,

          =
    where P_n,m is the quasi-normalized normalized Legendre associated function.

"""
function legendre_schmidt_quasi_normalized(ϕ::T, n_max::Number,
                                           m_max::Number = -1,
                                           ph_term::Bool = false) where T<:AbstractFloat

    (n_max < 0) && throw(ArgumentError("n_max must be positive."))

    ( (m_max < 0) || (m_max > n_max) ) && (m_max = n_max)

    P = zeros(T, n_max+1, m_max+1)
    legendre_schmidt_quasi_normalized!(P, ϕ, ph_term)
    return P
end

################################################################################
#                  Conventional Associated Legendre Function
################################################################################

"""
    function legendre_conventional!(P::AbstractMatrix, ϕ::Number, ph_term::Bool = false)

Compute the conventional associated Legendre function `P_n,m[cos(ϕ)]`. The
maximum degree and order that will be computed are given by the dimensions of
matrix `P`:

    maximum degree -> number of rows
    maximum order  -> number of columns

The result will be stored at matrix `P`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

"""
function legendre_conventional!(P::AbstractMatrix, ϕ::Number, ph_term::Bool = false)
    (rows, cols) = size(P)

    # Obtain the maximum degree and order that must be computed.
    n_max = rows - 1
    m_max = cols <= rows ? cols - 1 : n_max

    # Auxiliary variables to improve code performance.
    c = cos(ϕ)
    s = sqrt(1-c^2)

    s_fact = !ph_term ? +s : -s

    @inbounds for n = 0:n_max
        # Starting values.
        if n == 0
            P[0+1,0+1] = 1
            continue
        elseif n == 1
            P[1+1,0+1] = +c

            if cols > 1
                P[1+1,1+1] = +s_fact
            end

            continue
        end

        aux_n = 2n-1 # -> sqrt( (2n-1)*(2n-1) )

        for m = 0:n

            if n == m
                P[n+1,n+1] = s_fact*aux_n*P[n-1+1,n-1+1]
                continue
            end

            aux_nm = n-m # -> sqrt( (n-m)*(n-m) )
            a_nm   = aux_n / aux_nm
            b_nm   = (n+m-1) / aux_nm # -> sqrt( (n+m-1)*(n+m-1) ) / aux_nm

            # We assume that the matrix is not initialized. Hence, we must not
            # access elements on the upper triangle.
            if m != n-1
                P[n+1,m+1] = a_nm*c*P[n-1+1,m+1] - b_nm*P[n-2+1,m+1]
            else
                P[n+1,m+1] = a_nm*c*P[n-1+1,m+1]
            end

            # Check if the maximum desired order has been reached.
            m == m_max && break
        end
    end

    nothing
end

"""
    function legendre_conventional(ϕ::T, n_max::Number, m_max::Number = -1, ph_term::Bool = false) where T<:AbstractFloat

Compute the conventional associated Legendre function `P_n,m[cos(ϕ)]`. The
maximum degree that will be computed is `n_max` and the maximum order is
`m_max`. Notice that if `m_max` is higher than `n_max` or negative, than it is
set to `n_max`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.

"""
function legendre_conventional(ϕ::T, n_max::Number, m_max::Number = -1,
                               ph_term::Bool = false) where T<:AbstractFloat

    (n_max < 0) && throw(ArgumentError("n_max must be positive."))

    ( (m_max < 0) || (m_max > n_max) ) && (m_max = n_max)

    P = zeros(T, n_max+1, m_max+1)
    legendre_conventional!(P, ϕ, ph_term)
    return P
end
