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
#   Earth Gravitational Model (EGM).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/egm96.html
#   [1] http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/first_release.html
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-04-06: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export read_egm_coefs, EGM_Coefs, compute_g, compute_U

################################################################################
#                                    Types
################################################################################

"""
Structure to store the EGM coefficients.

"""

struct EGM_Coefs
    C::SparseMatrixCSC
    S::SparseMatrixCSC
    μ::Number
    R0::Number
end

# Serialization of the arguments in EGM_Coefs.
function getindex(egm_coefs::EGM_Coefs, ::Colon)
    egm_coefs.C, egm_coefs.S, egm_coefs.μ, egm_coefs.R0
end

################################################################################
#                                  Functions
################################################################################

"""
### function compute_g(egm_coefs::EGM_Coefs, r::Vector, n_max::Number)

Compute the gravity acceleration at position `r` using the EGM coefficients
`egm_coefs`. The maximum degree that will be used while computing the spherical
harmonics will be `n_max`.

##### Args

* egm_coefs: EGM coefficients.
* r: Position in ITRF in which the gravity will be computed [m].
* n_max: Maximum degree when computing the spherical harmonics.

##### Returns

A vector with the gravity acceleration represented in ITRF (Earth body-fixed
frame).

"""

function compute_g(egm_coefs::EGM_Coefs, r::Vector, n_max::Number)
    # Unpack EGM coefficients.
    C, S, μ, R0 = egm_coefs[:]

    # Check if n_max is bigger than the maximum degree available.
    if size(C,1)-1 < n_max
        n_max = size(C,1)-1
    end

    # Get the geocentric latitude and longitude.
    # ==========================================
    r_gc = norm(r)
    ρ_gc = sqrt(r[1]^2 + r[2]^2)
    ϕ_gc = atan2(r[3], ρ_gc)
    λ_gc = atan2(r[2], r[1])

    # Auxiliary variables
    # ===================

    # Check if we are at the poles.
    in_pole = ( abs(ρ_gc) < 1e-12 )

    # Tangent of the geocentric latitude.
    tan_ϕ = tan(ϕ_gc)

    # First derivative of the non-spherical portion of the grav. potential.
    # ====================================================================
    dUr = 0.0  # Derivative w.r.t. the radius.
    dUϕ = 0.0  # Derivative w.r.t. the geocentric latitude.
    dUλ = 0.0  # Derivative w.r.t. the geocentric longitude.

    # Compute the associated Legendre functions `P_n,m[sin(ϕ_gc)]`.
    P = legendre(ϕ_gc, n_max)

    # Compute the derivatives.
    for n = 2:n_max
        rn = (R0/r_gc)^n

        aux_dUr = 0.0
        aux_dUϕ = 0.0
        aux_dUλ = 0.0

        for m = 0:n
            sin_mλ = sin(m*λ_gc)
            cos_mλ = cos(m*λ_gc)

            CcSs_nm = C[n+1,m+1]*cos_mλ + S[n+1,m+1]*sin_mλ

            aux_dUr +=  P[n+1,m+1]*CcSs_nm
            aux_dUλ += m*P[n+1,m+1]*(S[n+1,m+1]*cos_mλ - C[n+1,m+1]*sin_mλ)

            if n == m
                # If n == m, then P_n,m+1 is 0.
                aux_dUϕ += -m*tan_ϕ*P[n+1,m+1]*CcSs_nm
            else
                # The factor `norm_coef` is not present in [1, p.  550].
                # However, this is necessary because the spherical harmonics of
                # EGM considers the fully normalized associated Legendre
                # functions, whereas [1] considers the conventional ones.
                norm_coef =
                    (m == 0) ? sqrt((n+m+1)*(n-m)/2) : sqrt((n+m+1)*(n-m))

                aux_dUϕ += (norm_coef*P[n+1,m+1+1]-m*tan_ϕ*P[n+1,m+1])*CcSs_nm
            end
        end

        dUr += rn*(n+1)*aux_dUr
        dUϕ += rn*aux_dUϕ
        dUλ += rn*aux_dUλ
    end

    dUr *= -μ/r_gc^2
    dUϕ *= +μ/r_gc
    dUλ *= +μ/r_gc

    # Compute the acceleration represented in the ITRF.
    #
    # If we are at the poles, then we must avoid divisions by 0.
    a_ITRF = !(in_pole) ?
        [ (dUr/r_gc - r[3]/(r_gc^2*ρ_gc)*dUϕ)*r[1] - (dUλ/ρ_gc^2)*r[2] - μ*r[1]/r_gc^3;
          (dUr/r_gc - r[3]/(r_gc^2*ρ_gc)*dUϕ)*r[2] - (dUλ/ρ_gc^2)*r[1] - μ*r[2]/r_gc^3;
          (dUr/r_gc - μ/r_gc^3)*r[3] + ρ_gc/r_gc^2*dUϕ ] :
        [ 0.0;
          0.0;
          (dUr/r_gc - μ/r_gc^3)*r[3] ]
end

"""
### function compute_U(egm_coefs::EGM_Coefs, ϕ_gc::Number, λ_gc::Number, r_gc::Number, n_max::Number)

Compute the gravity potential using the EGM coefficients `egm_coefs` at the
geocentric latitude `ϕ_gc`, geocentric longitude `λ_gc`, and geocentric radius
`r_gc`. The maximum degree that will be used while computing the spherical
harmonics will be `n_max`.

##### Args

* egm_coefs: EGM coefficients.
* ϕ_gc: Geocentric latitude [rad].
* λ_gc: Geocentric longitude [rad].
* r_gc: Geocentric radius [m].
* n_max: Maximum degree when computing the spherical harmonics.

##### Returns

The gravitational potential at specified location [J/kg].

"""

function compute_U(egm_coefs::EGM_Coefs,
                   ϕ_gc::Number,
                   λ_gc::Number,
                   r_gc::Number,
                   n_max::Number)
    # Unpack EGM coefficients.
    C, S, μ, R0 = egm_coefs[:]

    # Check if n_max is bigger than the maximum degree available.
    if size(C,1)-1 < n_max
        n_max = size(C,1)-1
    end

    U = 1.0

    # Compute the associated Legendre functions.
    P = legendre(ϕ_gc, 360)

    for n = 2:n_max
        aux_U = 0.0

        for m = 0:n
            sin_mλ = sin(m*λ_gc)
            cos_mλ = cos(m*λ_gc)

            aux_U += P[n+1,m+1]*(C[n+1,m+1]*cos_mλ + S[n+1,m+1]*sin_mλ)
        end

        U += (R0/r)^n*aux_U
    end

    U *= μ/r
end

"""
### function read_egm_coefs(filename::String)

Read the file `filename` with the EGM coefficients and create the structure
`EGM_Coefs` with them.

##### Args

* filename: The path to the file with the coefficients.
* μ: Earth gravitational constant with atmosphere [m³/s²].
* R0: Semi-major axis of the reference elipsoid [m].

##### Returns

A structure `EGM_Coefs` with the coefficients.

##### Remarks

The original file with EGM2008 coefficients use a Fortan formated float:

    -0.484169317366974D-03

This must be changed to a format understandable by Julia:

    -0.484169317366974e-03

"""

function read_egm_coefs(filename::String, μ::Number, R0::Number)
    raw = readdlm(filename)

    I = map(x->Int(x)+1, raw[:,1])
    J = map(x->Int(x)+1, raw[:,2])

    C = sparse(I, J, raw[:,3])
    S = sparse(I, J, raw[:,4])

    EGM_Coefs(C,S,μ,R0)
end

"""
### function read_egm_coefs(egm_version::Symbol)

Read the EGM coefficients from the bundled files. The variable `egm_version`
specify which EGM version must be used. The possible values are:

* `:EGM96`: Use EGM96 coefficients.
* `:EGM2008`: Use EGM2008 coefficients.

##### Args

* egm_version: (OPTIONAL) Select EGM version (**DEFAULT**: `:EGM2008`).

##### Returns

A structure `EGM_Coefs` with the coefficients.

"""

function read_egm_coefs(egm_version::Symbol = :EGM2008)
    dir = @__DIR__

    if egm_version == :EGM96
        filename = dir * "/coefficients/egm96"
        return read_egm_coefs(filename, 3986004.418e8, 6378137.0)
    elseif egm_version == :EGM2008
        filename = dir * "/coefficients/EGM2008_to360_ZeroTide"
        return read_egm_coefs(filename, 3986004.415e8, 6378136.3)
    else
        throw(ArgumentError("Unknown EGM version. Possible options are :EGM96 or :EGM2008."))
    end
end
