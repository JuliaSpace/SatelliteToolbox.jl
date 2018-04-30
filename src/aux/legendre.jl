#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
#       order normalised associated Legendre functions Journal of Geodesy,
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
#
# Changelog
#
# 2018-04-30: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   - Improve the structure of the source code.
#   - Add function to compute the Schmidt quasi-normalized associated Legendre
#     functions.
#   - **BREAKING**: The functions now compute the associated Legendre function
#     `P_n,m[cos(ϕ)]` instead of `P_n,m[sin(ϕ)]`.
#
# 2018-04-06: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export denormalize_legendre, legendre!, legendre

export legendre_fully_normalized!, legendre_fully_normalized
export legendre_schmidt_quasi_normalized!, legendre_schmidt_quasi_normalized

"""
### function denormalize_legendre(P_lm::Number, n::Number, m::Number)

Compute the conventional associated Legendre function from a fully normalized
value `P_nm` that has order `n` and degree `m`.

##### Args

* P_nm: Fully normalized Legendre associated function with degree `n` and order
        `m`.
* n: Degree of `P_nm`.
* m: Order of `P_nm`.

##### Returns

The conventional associated Legendre function.

"""

function denormalize_legendre(P_nm::Number, n::Number, m::Number)
    (n-m < 0) && throw(ArgumentError("n must be equal or bigger than m."))

    k = (m == 0) ? 1 : 2

    P_nm/sqrt(k*(2n+1)*factorial(n-m)/factorial(n+m))
end

"""
### function legendre!([N,] P::Matrix, ϕ::Number, ph_term::Bool = false)

Compute the associated Legendre function `P_n,m[cos(ϕ)]`. The maximum degree
and order that will be computed are given by the dimensions of matrix `P`.
Notice, however, that `P` must be a square matrix.

The result will be stored at matrix `P`.

The optional parameter `N` can be used to select the normalization. The
following values are valid:

* `Val{:full}`: Compute the fully normalized associated Legendre function (see
  `legendre_fully_normalized!`).
* `Val{:schmidt}`: Compute the Schmidt quasi-normalized associated Legendre
  function (see `legendre_schmidt_quasi_normalized!`).

If `N` is omitted, then the full normalization will be used.

##### Args

* N: (OPTIONAL) Choose the normalization (**DEFAULT**: none).
* P: Matrix that will store the computed associated Legendre function.
* ϕ: Angle [rad].
* ph_term: (OPTIONAL) Include the Condon-Shortley phase term `(-1)^m`
           (**DEFAULT** = `false`).

"""
legendre!(P::Matrix, ϕ::Number, ph_term::Bool = false) =
    legendre_fully_normalized!(P, ϕ, ph_term)

function legendre!(::Type{Val{:full}},
                   P::Matrix,
                   ϕ::Number,
                   ph_term::Bool = false)
    legendre_fully_normalized!(P, ϕ, ph_term)
end

function legendre!(::Type{Val{:schmidt}},
                   P::Matrix,
                   ϕ::Number,
                   ph_term::Bool = false)
    legendre_schmidt_quasi_normalized!(P, ϕ, ph_term)
end

"""
### function legendre([N,] ϕ::Number, n_max::Number, ph_term::Bool = false)

Compute the associated Legendre function `P_n,m[cos(ϕ)]`. The maximum degree
that will be computed is `n_max`.

The optional parameter `N` can be used to select the normalization. The
following values are valid:

* `Val{:full}`: Compute the fully normalized associated Legendre function (see
  `legendre_fully_normalized!`).
* `Val{:schmidt}`: Compute the Schmidt quasi-normalized associated Legendre
  function (see `legendre_schmidt_quasi_normalized!`).

If `N` is omitted, then the full normalization will be used.

##### Args

* N: (OPTIONAL) Choose the normalization (**DEFAULT**: none).
* ϕ: Angle [rad].
* n_max: The maximum degree that will be computed.
* ph_term: (OPTIONAL) Include the Condon-Shortley phase term `(-1)^m`
           (**DEFAULT** = `false`).

##### Returns

A square matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.

"""

legendre(ϕ::Number, n_max::Number, ph_term::Bool = false) =
    legendre_fully_normalized(ϕ, n_max, ph_term)

function legendre(::Type{Val{:full}},
                  ϕ::Number,
                  n_max::Number,
                  ph_term::Bool = false)
    legendre_fully_normalized(ϕ, n_max, ph_term)
end

function legendre(::Type{Val{:schmidt}},
                  ϕ::Number,
                  n_max::Number,
                  ph_term::Bool = false)
    legendre_schmidt_quasi_normalized(ϕ, n_max, ph_term)
end

"""
### function legendre_fully_normalized!(P::Matrix, ϕ::Number, ph_term::Bool = false)

Compute the fully normalized associated Legendre function `P_n,m[cos(ϕ)]`.
The maximum degree and order that will be computed are given by the dimensions
of matrix `P`. Notice, however, that `P` must be a square matrix.

The result will be stored at matrix `P`.

##### Args

* P: Matrix that will store the computed associated Legendre function.
* ϕ: Angle [rad].
* ph_term: (OPTIONAL) Include the Condon-Shortley phase term `(-1)^m`
           (**DEFAULT** = `false`).

##### Remarks

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

function legendre_fully_normalized!(P::Matrix, ϕ::Number, ph_term::Bool = false)
    (rows, cols) = size(P)

    # Check if the matrix P is a square matrix.
    (rows != cols) && throw(ArgumentError("P must be a square matrix."))

    # The matrix must be, at least, 2 rows and 2 columns.
    (rows < 2) && throw(ArgumentError("P must have at least 2 rows."))

    # Auxiliary variables to improve code performance.
    s = sin(ϕ)
    c = cos(ϕ)

    # Starting values.
    P[0+1,0+1] = 1
    P[1+1,0+1] = +sqrt(3)*c
    P[1+1,1+1] = -sqrt(3)*s

    (ph_term) && (P[1+1,1+1] *= -1)

    for n = 2:rows-1
        for m = 0:n-1
            a_nm = sqrt( ( (2n-1)*(2n+1) ) / ( (n-m)*(n+m) ) )
            b_nm = sqrt( ( (2n+1)*(n+m-1)*(n-m-1) ) / ( (n-m)*(n+m)*(2n-3) ) )
            P[n+1,m+1] = a_nm*c*P[n-1+1,m+1] - b_nm*P[n-2+1,m+1]
        end

        if ph_term
            P[n+1,n+1] = +s*sqrt( (2n+1)/(2n) )*P[n-1+1,n-1+1]
        else
            P[n+1,n+1] = -s*sqrt( (2n+1)/(2n) )*P[n-1+1,n-1+1]
        end
    end

    nothing
end

"""
### function legendre_fully_normalized(ϕ::Number, n_max::Number, ph_term::Bool = false)

Compute the fully normalized associated Legendre function `P_n,m[cos(ϕ)]`. The
maximum degree that will be computed is `n_max`.

##### Args

* ϕ: Angle [rad].
* n_max: The maximum degree that will be computed.
* ph_term: (OPTIONAL) Include the Condon-Shortley phase term `(-1)^m`
           (**DEFAULT** = `false`).

##### Returns

A square matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.

##### Remarks

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

function legendre_fully_normalized(ϕ::Number,
                                   n_max::Number,
                                   ph_term::Bool = false)
    (n_max < 1) && throw(ArgumentError("n_max must be at least 1."))

    P = zeros(n_max+1, n_max+1)
    legendre_fully_normalized!(P, ϕ, ph_term)
    P
end

"""
### function legendre_schmidt_quasi_normalized!(P::Matrix, ϕ::Number, ph_term::Bool = false)

Compute the Schmidt quasi-normalized associated Legendre function
`P_n,m[cos(ϕ)]` [3,4].  The maximum degree and order that will be computed are
given by the dimensions of matrix `P`. Notice, however, that `P` must be a
square matrix.

The result will be stored at matrix `P`.

##### Args

* P: Matrix that will store the computed associated Legendre function.
* ϕ: Angle [rad].
* ph_term: (OPTIONAL) Include the Condon-Shortley phase term `(-1)^m`
           (**DEFAULT** = `false`).

##### Remarks

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

function legendre_schmidt_quasi_normalized!(P::Matrix,
                                            ϕ::Number,
                                            ph_term::Bool = false)
    (rows, cols) = size(P)

    # Check if the matrix P is a square matrix.
    (rows != cols) && throw(ArgumentError("P must be a square matrix."))

    # The matrix must be, at least, 2 rows and 2 columns.
    (rows < 2) && throw(ArgumentError("P must have at least 2 rows."))

    # Auxiliary variables to improve code performance.
    s = sin(ϕ)
    c = cos(ϕ)

    # Starting values.
    P[0+1,0+1] = 1
    P[1+1,0+1] = +c
    P[1+1,1+1] = -s

    (ph_term) && (P[1+1,1+1] *= -1)

    for n = 2:rows-1
        for m = 0:n-1
            aux = (n-m)*(n+m)
            a_nm = sqrt( ( (2n-1)*(2n-1) ) / aux )
            b_nm = sqrt( ( (n+m-1)*(n-m-1) ) / aux )
            P[n+1,m+1] = a_nm*c*P[n-1+1,m+1] - b_nm*P[n-2+1,m+1]
        end

        if ph_term
            P[n+1,n+1] = +s*sqrt( (2n-1)/(2n) )*P[n-1+1,n-1+1]
        else
            P[n+1,n+1] = -s*sqrt( (2n-1)/(2n) )*P[n-1+1,n-1+1]
        end
    end

    nothing
end

"""
### function legendre_schmidt_quasi_normalized(ϕ::Number, n_max::Number, ph_term::Bool = false)

Compute the fully normalized associated Legendre function `P_n,m[cos(ϕ)]`.
The maximum degree that will be computed is `n_max`.

##### Args

* ϕ: Angle [rad].
* n_max: The maximum degree that will be computed.
* ph_term: (OPTIONAL) Include the Condon-Shortley phase term `(-1)^m`
           (**DEFAULT** = `false`).

##### Returns

A square matrix with the Legendre associated functions `P_n,m[cos(ϕ)]`.

##### Remarks

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

function legendre_schmidt_quasi_normalized(ϕ::Number,
                                           n_max::Number,
                                           ph_term::Bool = false)
    (n_max < 1) && throw(ArgumentError("n_max must be at least 1."))

    P = zeros(n_max+1, n_max+1)
    legendre_schmidt_quasi_normalized!(P, ϕ, ph_term)
    P
end
