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
order that will be computed are given by the dimensions of this matrix. Notice,
however, that `dP` must be a square matrix.

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
dlegendre!(dP::AbstractMatrix,
           ϕ::Number,
           P::AbstractMatrix,
           ph_term::Bool = false) =
    dlegendre_fully_normalized!(dP, ϕ, P, ph_term)

dlegendre!(::Type{Val{:full}},
           dP::AbstractMatrix,
           ϕ::Number,
           P::AbstractMatrix,
           ph_term::Bool = false) =
    dlegendre_fully_normalized!(dP, ϕ, P, ph_term)

dlegendre!(::Type{Val{:schmidt}},
           dP::AbstractMatrix,
           ϕ::Number,
           P::AbstractMatrix,
           ph_term::Bool = false) =
    dlegendre_schmidt_quasi_normalized!(dP, ϕ, P, ph_term)

dlegendre!(::Type{Val{:conv}},
           dP::AbstractMatrix,
           ϕ::Number,
           P::AbstractMatrix,
           ph_term::Bool = false) =
    dlegendre_conventional!(dP, ϕ, P, ph_term)

"""
    function dlegendre([N,] ϕ::Number, P::Matrix, ph_term::Bool)

Compute the first-order derivative of the associated Legendre function
`P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

This algorithm needs the matrix `P` with the associated Legendre function. This
can be computed using the function `legendre`. The maximum order of the computed
derivatives will be selected according to the dimensions of the matrix `P`.
Notice that this matrix must be computed using the same normalization (see `T`)
as the one selected here.

The optional parameter `N` can be used to select the normalization. The
following values are valid:

* `Val{:full}`: Compute the fully normalized associated Legendre function (see
                `dlegendre_fully_normalized`).
* `Val{:schmidt}`: Compute the Schmidt quasi-normalized associated Legendre
                   function (see `dlegendre_schmidt_quasi_normalized`).
* `Val{:conv}`: Compute the conventional associated Legendre function (see
                `dlegendre_conventional!`).

If `N` is omitted, then the full normalization will be used (`Val{:full}`).

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the first-order derivative of the Legendre associated functions
`P_n,m[cos(ϕ)]`.

"""
dlegendre(ϕ::Number, P::Matrix, ph_term::Bool) =
    dlegendre_fully_normalized(ϕ, P, ph_term)

dlegendre(::Type{Val{:full}}, ϕ::Number, P::Matrix, ph_term::Bool) =
    dlegendre_fully_normalized(ϕ, P, ph_term)

dlegendre(::Type{Val{:schmidt}}, ϕ::Number, P::Matrix, ph_term::Bool) =
    dlegendre_schmidt_quasi_normalized(ϕ, P, ph_term)

dlegendre(::Type{Val{:conv}}, ϕ::Number, P::Matrix, ph_term::Bool) =
    dlegendre_conventional(ϕ, P, ph_term)

"""
    function dlegendre([N,] ϕ::Number, n_max::Number, ph_term::Bool = false)

Compute the first-order derivative of the associated Legendre function
`P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

The maximum degree that will be computed is `n_max`.

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
dlegendre(ϕ::Number, n_max::Number, ph_term::Bool = false) =
    dlegendre_fully_normalized(ϕ, n_max, ph_term)

dlegendre(::Type{Val{:full}},
          ϕ::Number,
          n_max::Number,
          ph_term::Bool = false) =
    dlegendre_fully_normalized(ϕ, n_max, ph_term)

dlegendre(::Type{Val{:schmidt}},
          ϕ::Number,
          n_max::Number,
          ph_term::Bool = false) =
    dlegendre_schmidt_quasi_normalized(ϕ, n_max, ph_term)

dlegendre(::Type{Val{:conv}},
          ϕ::Number,
          n_max::Number,
          ph_term::Bool = false) =
    dlegendre_conventional(ϕ, n_max, ph_term)

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
order that will be computed are given by the dimensions of this matrix. Notice,
however, that `dP` must be a square matrix.

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
function dlegendre_fully_normalized!(dP::AbstractMatrix,
                                     ϕ::Number,
                                     P::AbstractMatrix,
                                     ph_term::Bool = false)
    (rows, cols) = size(dP)

    # Check if the matrix P is a square matrix.
    (rows != cols) && throw(ArgumentError("dP must be a square matrix."))

    # The matrix must be, at least, 2 rows and 2 columns.
    (rows < 2) && throw(ArgumentError("dP must have at least 2 rows."))

    # The matrix dP must not have more rows or columns than P.
    (rows > size(P,1)) && throw(ArgumentError("dP must not have more rows than P."))
    (cols > size(P,2)) && throw(ArgumentError("dP must not have more columns than P."))

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

    @inbounds for n = 1:rows-1
        for m = 0:n
            if m == 0
                aux  = sqrt(n*(n+1)/2)
                a_nm = +0.5*aux
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
                a_nm = +0.5*sqrt(2n*(n+1))
                dP[n+1,1+1] = a_nm*P[n+1,1-1+1]

                # Only compute `b_nm` if `n > 1`. Otherwise, we could access an
                # invalid memory region if `P` is 2x2.
                if n > 1
                    b_nm = -0.5*sqrt((n+2)*(n-1))
                    dP[n+1,1+1] += b_nm*P[n+1,1+1+1]
                end

            elseif n != m
                a_nm = +0.5*sqrt((n+m)*(n-m+1))
                b_nm = -0.5*sqrt((n+m+1)*(n-m))

                dP[n+1,m+1] = a_nm*P[n+1,m-1+1] + b_nm*P[n+1,m+1+1]
            else
                a_nm = +0.5*sqrt((n+m)*(n-m+1))

                dP[n+1,m+1] = a_nm*P[n+1,m-1+1]
            end

            dP[n+1,m+1] *= fact
        end
    end

    nothing
end

"""
    function dlegendre_fully_normalized(ϕ::Number, P::AbstractMatrix, ph_term = false)

Compute the first-order derivative of the Schmidt fully normalized associated
Legendre function `P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

This algorithm needs the matrix `P` with the Schmidt fully normalized associated
Legendre function. This can be computed using the function
`legendre_fully_normalized`. The maximum order of the computed
derivatives will be selected according to the dimensions of the matrix `P`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the first-order derivative of the Legendre associated functions
`P_n,m[cos(ϕ)]`.

# Remarks

The user is responsible to pass a matrix `P` with the correct values. For
example, if `ph_term` is `true`, then `P` must also be computed with `ph_term`
set to `true`.

"""
function dlegendre_fully_normalized(ϕ::Number, P::AbstractMatrix, ph_term = false)
    dP = zero(P)
    dlegendre_fully_normalized!(dP,ϕ,P,ph_term)
    dP
end

"""
    function dlegendre_fully_normalized(ϕ::Number, n_max::Number, ph_term::Bool = false)

Compute the first-order derivative of the Schmidt fully normalized associated
Legendre function `P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

The maximum degree that will be computed is `n_max`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the first-order derivative of the Legendre associated functions
`P_n,m[cos(ϕ)]`.

"""
function dlegendre_fully_normalized(ϕ::Number,
                                    n_max::Number,
                                    ph_term::Bool = false)
    # First, compute the matrix with the associated Legendre functions.
    P = legendre_fully_normalized(ϕ, n_max, ph_term)

    # Now, compute and return the time-derivative of the associated Legendre
    # functions.
    dlegendre_fully_normalized(ϕ, P, ph_term)
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
order that will be computed are given by the dimensions of this matrix. Notice,
however, that `dP` must be a square matrix.

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
dlegendre_schmidt_quasi_normalized!(dP::AbstractMatrix,
                                    ϕ::Number,
                                    P::AbstractMatrix,
                                    ph_term::Bool = false) =

    # The algorithm to compute the first-order derivative using Schmidt
    # normalizations is precisely the same as the one that computes using full
    # normalization.
    dlegendre_fully_normalized!(dP, ϕ, P, ph_term)

"""
    function dlegendre_schmidt_quasi_normalized(ϕ::Number, P::AbstractMatrix, ph_term = false)

Compute the first-order derivative of the Schmidt quasi-normalized associated
Legendre function `P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

This algorithm needs the matrix `P` with the Schmidt quasi-normalized associated
Legendre function. This can be computed using the function
`legendre_schmidt_quasi_normalized`. The maximum order of the computed
derivatives will be selected according to the dimensions of the matrix `P`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the first-order derivative of the Legendre associated functions
`P_n,m[cos(ϕ)]`.

# Remarks

The user is responsible to pass a matrix `P` with the correct values. For
example, if `ph_term` is `true`, then `P` must also be computed with `ph_term`
set to `true`.

"""
function dlegendre_schmidt_quasi_normalized(ϕ::Number,
                                            P::AbstractMatrix,
                                            ph_term = false)
    dP = zero(P)
    dlegendre_schmidt_quasi_normalized!(dP, ϕ, P, ph_term)
    dP
end

"""
    function dlegendre_schmidt_quasi_normalized(ϕ::Number, n_max::Number, ph_term::Bool = false)

Compute the first-order derivative of the Schmidt quasi-normalized associated
Legendre function `P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

The maximum degree that will be computed is `n_max`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the first-order derivative of the Legendre associated functions
`P_n,m[cos(ϕ)]`.

"""
function dlegendre_schmidt_quasi_normalized(ϕ::Number,
                                            n_max::Number,
                                            ph_term::Bool = false)
    # First, compute the matrix with the associated Legendre functions.
    P = legendre_schmidt_quasi_normalized(ϕ, n_max, ph_term)

    # Now, compute and return the time-derivative of the associated Legendre
    # functions.
    dlegendre_schmidt_quasi_normalized(ϕ, P, ph_term)
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
order that will be computed are given by the dimensions of this matrix. Notice,
however, that `dP` must be a square matrix.

This algorithm needs the matrix `P` with the conventional associated Legendre
function. This can be computed using the function `legendre_conventional`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Remarks

The user is responsible to pass a matrix `P` with the correct values. For
example, if `ph_term` is `true`, then `P` must also be computed with `ph_term`
set to `true`.

"""
function dlegendre_conventional!(dP::AbstractMatrix,
                                 ϕ::Number,
                                 P::AbstractMatrix,
                                 ph_term::Bool = false)
    (rows, cols) = size(dP)

    # Check if the matrix P is a square matrix.
    (rows != cols) && throw(ArgumentError("dP must be a square matrix."))

    # The matrix must be, at least, 2 rows and 2 columns.
    (rows < 2) && throw(ArgumentError("dP must have at least 2 rows."))

    # The matrix dP must not have more rows or columns than P.
    (rows > size(P,1)) && throw(ArgumentError("dP must not have more rows than P."))
    (cols > size(P,2)) && throw(ArgumentError("dP must not have more columns than P."))

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

    @inbounds for n = 1:rows-1
        for m = 0:n
            if m == 0
                dP[n+1,m+1] = -P[n+1,1+1]
            elseif n != m
                dP[n+1,m+1] = 0.5*( (n+m)*(n-m+1)*P[n+1,m-1+1] - P[n+1,m+1+1] )
            else
                dP[n+1,m+1] = 0.5*( (n+m)*(n-m+1)*P[n+1,m-1+1] )
            end

            dP[n+1,m+1] *= fact
        end
    end

    nothing
end

"""
    function dlegendre_conventional(ϕ::Number, P::AbstractMatrix, ph_term = false)

Compute the first-order derivative of the conventional associated Legendre
function `P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

This algorithm needs the matrix `P` with the conventional associated Legendre
function. This can be computed using the function
`legendre_schmidt_quasi_normalized`. The maximum order of the computed
derivatives will be selected according to the dimensions of the matrix `P`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the first-order derivative of the Legendre associated functions
`P_n,m[cos(ϕ)]`.

# Remarks

The user is responsible to pass a matrix `P` with the correct values. For
example, if `ph_term` is `true`, then `P` must also be computed with `ph_term`
set to `true`.

"""
function dlegendre_conventional(ϕ::Number, P::AbstractMatrix, ph_term = false)
    dP = zero(P)
    dlegendre_conventional!(dP, ϕ, P, ph_term)
    dP
end

"""
    function dlegendre_conventional(ϕ::Number, n_max::Number, ph_term::Bool = false)

Compute the first-order derivative of the conventional associated Legendre
function `P_n,m[cos(ϕ)]` w.r.t. `ϕ` [rad]:

    dP_n,m[cos(ϕ)]
    --------------
          dϕ

The maximum degree that will be computed is `n_max`.

If `ph_term` is set to `true`, then the Condon-Shortley phase term `(-1)ᵐ` will
be included. If `ph_term` is not present, then it defaults to `false`.

# Returns

A matrix with the first-order derivative of the Legendre associated functions
`P_n,m[cos(ϕ)]`.

"""
function dlegendre_conventional(ϕ::Number,
                                n_max::Number,
                                ph_term::Bool = false)
    # First, compute the matrix with the associated Legendre functions.
    P = legendre_conventional(ϕ, n_max, ph_term)

    # Now, compute and return the time-derivative of the associated Legendre
    # functions.
    dlegendre_conventional(ϕ, P, ph_term)
end
