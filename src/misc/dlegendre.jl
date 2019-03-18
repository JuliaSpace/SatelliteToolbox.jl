#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions related to the first-order derivative of associated Legendre
#   functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Du, J., Chen, C., Lesur, V., and Wang, L (2015). Non-singular spherical
#       harmonic expressions of geomagnetic vector and gradient tensor fields in
#       the local north-oriented reference frame. Geoscientific Model
#       Development, 8, pp. 1979-1990.
#
#   [2] Ilk, K. H.: Ein eitrag zur Dynamik ausgedehnter
#       Körper-Gravitationswechselwirkung, Deutsche Geodätische Kommis- sion.
#       Reihe C, Heft Nr. 288, München, 1983.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export dlegendre!, dlegendre
export dlegendre_fully_normalized!, dlegendre_fully_normalized
export dlegendre_schmidt_quasi_normalized!, dlegendre_schmidt_quasi_normalized
export dlegendre_conventional!, dlegendre_conventional

"""
    function dlegendre!([N,] dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix, ph_term::Bool = false)

Compute the first-order derivative of the associated Legendre function
`P_n,m[x]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

The derivatives will be stored in the matrix `dP`. Hence, the maximum degree and
order that will be computed are given by the dimensions of this matrix.

This algorithm needs the matrix `P` with the associated Legendre function. This
can be computed using the function `legendre`. Notice that this matrix must be
computed using the same normalization (see `N`) as the one selected here.

The optional parameter `N` can be used to select the normalization. The
following values are valid:

* `Val{:full}`: Compute the fully normalized associated Legendre function (see
                `dlegendre_fully_normalized!`).
* `Val{:schmidt}`: Compute the Schmidt quasi-normalized associated Legendre
                   function (see `dlegendre_schmidt_quasi_normalized!`).
* `Val{:conv}`: Compute the conventional associated Legendre function (see
                `dlegendre_conventional!`).

If `N` is omitted, then the full normalization will be used (`Val{:full}`).

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

"""
dlegendre!(dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix, ph_term::Bool = false) =
    dlegendre_fully_normalized!(dP, ϕ, P, ph_term)

dlegendre!(::Type{Val{:full}}, dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix,
           ph_term::Bool = false) =
    dlegendre_fully_normalized!(dP, ϕ, P, ph_term)

dlegendre!(::Type{Val{:schmidt}}, dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix,
           ph_term::Bool = false) =
    dlegendre_schmidt_quasi_normalized!(dP, ϕ, P, ph_term)

dlegendre!(::Type{Val{:conv}}, dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix,
           ph_term::Bool = false) =
    dlegendre_conventional!(dP, ϕ, P, ph_term)

"""
    function dlegendre([N,] ϕ::Number, n_max::Integer, m_max::Integer = -1, ph_term::Bool = false)

Compute the first-order derivative of the associated Legendre function
`P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

The maximum degree that will be computed is `n_max` and the maximum order is
`m_max`. Notice that if `m_max` is higher than `n_max` or negative, than it is
set to `n_max`.

The optional parameter `N` can be used to select the normalization. The
following values are valid:

* `Val{:full}`: Compute the fully normalized associated Legendre function (see
                `legendre_fully_normalized`).
* `Val{:schmidt}`: Compute the Schmidt quasi-normalized associated Legendre
                   function (see `legendre_schmidt_quasi_normalized`).
* `Val{:conv}`: Compute the conventional associated Legendre function (see
                `dlegendre_conventional!`).

If `N` is omitted, then the full normalization will be used (`Val{:full}`).

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the first-order derivative of the Legendre associated functions
`P_n,m[cos(ϕ)]`.

"""
dlegendre(ϕ::Number, n_max::Integer, m_max::Integer = -1, ph_term::Bool = false) =
    dlegendre_fully_normalized(float(ϕ), n_max, m_max, ph_term)

dlegendre(::Type{Val{:full}}, ϕ::Number, n_max::Integer, m_max::Integer = -1,
          ph_term::Bool = false) =
    dlegendre_fully_normalized(float(ϕ), n_max, m_max, ph_term)

dlegendre(::Type{Val{:schmidt}}, ϕ::Number, n_max::Integer, m_max::Integer = -1,
          ph_term::Bool = false) =
    dlegendre_schmidt_quasi_normalized(float(ϕ), n_max, m_max, ph_term)

dlegendre(::Type{Val{:conv}}, ϕ::Number, n_max::Integer, m_max::Integer = -1,
          ph_term::Bool = false) =
    dlegendre_conventional(float(ϕ), n_max, m_max, ph_term)

################################################################################
#                Fully Normalized Associated Legendre Functions
################################################################################

"""
    function dlegendre_fully_normalized!(dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix, ph_term::Bool = false)

Compute the first-order derivative of the fully normalized associated Legendre
function `P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

The derivatives will be stored in the matrix `dP`. Hence, the maximum degree and
order that will be computed are given by the dimensions of this matrix.

This algorithm needs the matrix `P` with the fully normalized associated
Legendre function. This can be computed using the function
`legendre_fully_normalized`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Remarks

The user is responsible to pass a matrix `P` with the correct values. For
example, if `ph_term` is `true`, then `P` must also be computed with `ph_term`
set to `true`.

"""
function dlegendre_fully_normalized!(dP::AbstractMatrix, ϕ::Number,
                                     P::AbstractMatrix, ph_term::Bool = false)

    # Get the maximum degree and order of `P`.
    (rows_P,  cols_P) = size(P)
    n_max_P = rows_P - 1
    m_max_P = cols_P <= rows_P ? cols_P - 1 : n_max_P

    # Compute the maximum degree and order of `dP`.
    (rows, cols) = size(dP)
    n_max = min(rows - 1, n_max_P)
    m_max = min(m_max_P - 1, cols <= rows ? cols - 1 : n_max)

    # The derivative is compute using the following equation [1, p. 1981]:
    #
    #   dP(n,m)
    #   ------- = a_nm.P(n,m-1) + b_nm.P(n,m+1) ,
    #     dϕ
    #
    #   a_nm = +0.5*sqrt(n+m)*sqrt(n-m+1)*sqrt(C_{m}/C_{m-1})
    #   b_nm = -0.5*sqrt(n+m+1)*sqrt(n-m)*sqrt(C_{m}/C_{m+1})
    #
    #           | 1, m  = 0
    #   C_{m} = |
    #           | 2, m != 0
    #
    # NOTE: The conversion of the coefficients between the full normalization
    # and the Schmidt normalization is performed using:
    #
    #   sqrt(2n+1) ,
    #
    # which depends only on `n`. Since the derivative equation only has terms
    # related to the order `n`, then the same algorithm will work for both full
    # normalization and Schmidt normalization.
    #
    # TODO: This algorithm is based on eq. Z.1.44 of [2]. However, it was
    # verified that it does not provide the correct sign when ϕ ∈ [π, 2π]. This
    # makes sense because the algorithm uses only the values of the
    # coefficients, which are equal for ϕ and -ϕ. However, the derivative w.r.t.
    # ϕ does change. This hack was used so that the values are correct, but
    # further verification is needed.
    #
    # In fact, in [2, p. 119], it is mentioned that `0 <= ϕ <= π`.  However,
    # further clarification is required.

    ϕ    = mod(ϕ,2*π)
    fact = (ϕ > π) ? -1 : 1

    ph_term && (fact *= -1)

    dP[0+1,0+1] = 0

    m_max < 0 && return nothing

    @inbounds for n = 1:n_max
        for m = 0:n
            if m == 0
                aux  = sqrt(n*(n+1)/2)
                a_nm = aux/2
                b_nm = -a_nm

                # Notice that [1, p. 1985]:
                #
                #                  m
                #   P_(n,-m) = (-1).P_(n,m) ,
                #
                dP[n+1,0+1] = -a_nm*P[n+1,1+1] + b_nm*P[n+1,0+1+1]

            # We should consider the case `m == 1` separately from `n == m`
            # because of the coefficient `C_{m}`.
            elseif m == 1
                a_nm = sqrt(2n*(n+1))/2
                dP[n+1,1+1] = a_nm*P[n+1,1-1+1]

                # Only compute `b_nm` if `n > 1`. Otherwise, we could access an
                # invalid memory region if `P` is 2x2.
                if n > 1
                    b_nm = -sqrt((n+2)*(n-1))/2
                    dP[n+1,1+1] += b_nm*P[n+1,1+1+1]
                end

            elseif n != m
                a_nm = +sqrt((n+m)*(n-m+1))/2
                b_nm = -sqrt((n+m+1)*(n-m))/2

                dP[n+1,m+1] = a_nm*P[n+1,m-1+1] + b_nm*P[n+1,m+1+1]
            else
                a_nm = +sqrt((n+m)*(n-m+1))/2

                dP[n+1,m+1] = a_nm*P[n+1,m-1+1]
            end

            dP[n+1,m+1] *= fact

            # Check if the maximum desired order has been reached.
            #
            # Notice that we can compute the term `dP[m_max + 1, m_max + 1]` if
            # there is space in `dP`.
            if m >= m_max
                if (n != m_max + 1) || (cols - 1 < n)
                    break
                end
            end
        end
    end

    return nothing
end

"""
    function dlegendre_fully_normalized(ϕ::T, n_max::Integer, m_max::Integer = -1, ph_term::Bool = false) where T<:AbstractFloat

Compute the first-order derivative of the Schmidt fully normalized associated
Legendre function `P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

The maximum degree that will be computed is `n_max` and the maximum order is
`m_max`. Notice that if `m_max` is higher than `n_max` or negative, than it is
set to `n_max`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the first-order derivative of the Legendre associated functions
`P_n,m[cos(ϕ)]`.

"""
function dlegendre_fully_normalized(ϕ::T, n_max::Integer, m_max::Integer = -1,
                                    ph_term::Bool = false) where T<:AbstractFloat

    ( (m_max < 0) || (m_max > n_max) ) && (m_max = n_max)

    # Check if we need to compute and additional degree in `P` to provide the
    # desire order in `dP`.
    if n_max == m_max
        n_max_P = m_max_P = n_max
    else
        n_max_P = n_max
        m_max_P = m_max + 1
    end

    # First, compute the matrix with the associated Legendre functions.
    P = legendre_fully_normalized(ϕ, n_max_P, m_max_P, ph_term)

    # Now, compute and return the time-derivative of the associated Legendre
    # functions.
    dP = zeros(T, n_max + 1, m_max + 1)
    dlegendre_fully_normalized!(dP, ϕ, P, ph_term)
    return dP, P
end

################################################################################
#            Schmidt Quasi-Normalized Associated Legendre Functions
################################################################################

"""
    function dlegendre_schmidt_quasi_normalized!(dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix, ph_term::Bool = false)

Compute the first-order derivative of the Schmidt quasi-normalized associated
Legendre function `P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

The derivatives will be stored in the matrix `dP`. Hence, the maximum degree and
order that will be computed are given by the dimensions of this matrix.

This algorithm needs the matrix `P` with the Schmidt quasi-normalized associated
Legendre function. This can be computed using the function
`legendre_schmidt_quasi_normalized`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Remarks

The user is responsible to pass a matrix `P` with the correct values. For
example, if `ph_term` is `true`, then `P` must also be computed with `ph_term`
set to `true`.

"""
dlegendre_schmidt_quasi_normalized!(dP::AbstractMatrix, ϕ::Number,
                                    P::AbstractMatrix, ph_term::Bool = false) =

    # The algorithm to compute the first-order derivative using Schmidt
    # normalizations is precisely the same as the one that computes using full
    # normalization.
    dlegendre_fully_normalized!(dP, ϕ, P, ph_term)

"""
    function dlegendre_schmidt_quasi_normalized(ϕ::T, n_max::Integer, m_max::Integer = -1, ph_term::Bool = false) where T<:AbstractFloat

Compute the first-order derivative of the Schmidt quasi-normalized associated
Legendre function `P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

The maximum degree that will be computed is `n_max` and the maximum order is
`m_max`. Notice that if `m_max` is higher than `n_max` or negative, than it is
set to `n_max`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the first-order derivative of the Legendre associated functions
`P_n,m[cos(ϕ)]`.

"""
function dlegendre_schmidt_quasi_normalized(ϕ::T, n_max::Integer,
                                            m_max::Integer = -1,
                                            ph_term::Bool = false) where T<:AbstractFloat

    ( (m_max < 0) || (m_max > n_max) ) && (m_max = n_max)

    # Check if we need to compute and additional degree in `P` to provide the
    # desire order in `dP`.
    if n_max == m_max
        n_max_P = m_max_P = n_max
    else
        n_max_P = n_max
        m_max_P = m_max + 1
    end

    # First, compute the matrix with the associated Legendre functions.
    P = legendre_schmidt_quasi_normalized(ϕ, n_max_P, m_max_P, ph_term)

    # Now, compute and return the time-derivative of the associated Legendre
    # functions.
    dP = zeros(T, n_max + 1, m_max + 1)
    dlegendre_schmidt_quasi_normalized!(dP, ϕ, P, ph_term)
    return dP, P
end

################################################################################
#                  Conventional Associated Legendre Function
################################################################################

"""
    function dlegendre_conventional!(dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix, ph_term::Bool = false)

Compute the first-order derivative of the conventional associated Legendre
function `P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

The derivatives will be stored in the matrix `dP`. Hence, the maximum degree and
order that will be computed are given by the dimensions of this matrix.

This algorithm needs the matrix `P` with the conventional associated Legendre
function. This can be computed using the function `legendre_conventional`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Remarks

The user is responsible to pass a matrix `P` with the correct values. For
example, if `ph_term` is `true`, then `P` must also be computed with `ph_term`
set to `true`.

"""
function dlegendre_conventional!(dP::AbstractMatrix, ϕ::Number, P::AbstractMatrix,
                                 ph_term::Bool = false)

    # Get the maximum degree and order of `P`.
    (rows_P,  cols_P) = size(P)
    n_max_P = rows_P - 1
    m_max_P = cols_P <= rows_P ? cols_P - 1 : n_max_P

    # Compute the maximum degree and order of `dP`.
    (rows, cols) = size(dP)
    n_max = min(rows - 1, n_max_P)
    m_max = min(m_max_P - 1, cols <= rows ? cols - 1 : n_max)

    # The derivative is computed using the following equation [1, p. 1981]:
    #
    #   dP(n,m)
    #   ------- = 0.5.( (n+m)(n-m+1)P(n,m-1) - P(n,m+1) ) ,
    #     dθ
    #
    # TODO: This algorithm is based on eq. Z.1.44 of [2]. However, it was
    # verified that it does not provide the correct sign when ϕ ∈ [π, 2π]. This
    # makes sense because the algorithm uses only the values of the
    # coefficients, which are equal for ϕ and -ϕ. However, the derivative w.r.t.
    # ϕ does change. This hack was used so that the values are correct, but
    # further verification is needed.
    #
    # In fact, in [2, p. 119], it is mentioned that `0 <= ϕ <= π`.  However,
    # further clarification is required.

    ϕ    = mod(ϕ,2*π)
    fact = (ϕ > π) ? -1 : 1

    ph_term && (fact *= -1)

    dP[0+1,0+1] = 0

    m_max < 0 && return nothing

    @inbounds for n = 1:n_max
        for m = 0:n
            if m == 0
                dP[n+1,m+1] = -P[n+1,1+1]
            elseif n != m
                dP[n+1,m+1] = ( (n+m)*(n-m+1)*P[n+1,m-1+1] - P[n+1,m+1+1] )/2
            else
                dP[n+1,m+1] = ( (n+m)*(n-m+1)*P[n+1,m-1+1] )/2
            end

            dP[n+1,m+1] *= fact

            # Check if the maximum desired order has been reached.
            #
            # Notice that we can compute the term `dP[m_max + 1, m_max + 1]` if
            # there is space in `dP`.
            if m >= m_max
                if (n != m_max + 1) || (cols - 1 < n)
                    break
                end
            end
        end
    end

    nothing
end

"""
    function dlegendre_conventional(ϕ::Number, n_max::Integer, m_max::Integer = -1, ph_term::Bool = false)

Compute the first-order derivative of the conventional associated Legendre
function `P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

The maximum degree that will be computed is `n_max` and the maximum order is
`m_max`. Notice that if `m_max` is higher than `n_max` or negative, than it is
set to `n_max`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the first-order derivative of the Legendre associated functions
`P_n,m[cos(ϕ)]`.

"""
function dlegendre_conventional(ϕ::T, n_max::Integer, m_max::Integer = -1,
                                ph_term::Bool = false) where T<:AbstractFloat

    ( (m_max < 0) || (m_max > n_max) ) && (m_max = n_max)

    # Check if we need to compute and additional degree in `P` to provide the
    # desire order in `dP`.
    if n_max == m_max
        n_max_P = m_max_P = n_max
    else
        n_max_P = n_max
        m_max_P = m_max + 1
    end

    # First, compute the matrix with the associated Legendre functions.
    P = legendre_conventional(ϕ, n_max_P, m_max_P, ph_term)

    # Now, compute and return the time-derivative of the associated Legendre
    # functions.
    dP = zeros(T, n_max + 1, m_max + 1)
    dlegendre_conventional!(dP, ϕ, P, ph_term)
    return dP, P
end
