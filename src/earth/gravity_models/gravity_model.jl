#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   General functions to compute the gravity field potential and gravity
#   acceleration.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#   [2] http://icgem.gfz-potsdam.de/home
#   [3] Delgado, M. R (2008). Aspherical Gravitational Field I: Modeling the
#       Space Environment. Universidad Poletécnica de Madrid, Madrid, Spain.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export create_gravity_model_coefs, compute_dU, compute_U, compute_g

################################################################################
#                               Public Functions
################################################################################

"""
    function create_gravity_model_coefs(icgem::ICGEM)

Return an instance of the structure `GravityModel_Coefs` based on the
information obtained from an ICGEM file in `icgem` (see `parse_icgem`).

"""
function create_gravity_model_coefs(icgem::ICGEM)

    legendre_norm = (icgem.norm == :fully_normalized) ? :full : :conv

    # Create and return the gravity model coefficients.
    gm_coefs::GravityModel_Coefs{Float64} =
        GravityModel_Coefs(name          = icgem.modelname,
                           μ             = icgem.gravity_constant,
                           R0            = icgem.radius,
                           n_max         = icgem.max_degree,
                           C             = copy(icgem.Clm),
                           S             = copy(icgem.Slm),
                           legendre_norm = legendre_norm)
end

"""
    function compute_dU(gm_coefs::GravityModel_Coefs{T}, r::AbstractVector, n_max::Number = -1, m_max::Number = -1) where T<:Number

Compute the derivatives w.r.t. the spherical coordinates of the gravitational
field (`∂U/∂r`, `∂U/∂ϕ`, `∂U/∂λ`) defined by the coefficients `gm_coefs` (see
`GravityModel_Coefs`) at the position `r` [m] in ITRF. The maximum degree that
will be used while computing the spherical harmonics will be `n_max` and the
maximum order is `m_max`.

If `n_max` is negative, then the maximum available degree will be used. If
`n_max` is omitted, then it defaults to 0.

If `m_max` is negative or if it is greater than `n_max`, then it will be set to
`n_max`. If `m_max` is omitted, then it defaults to 0.

!!! info

    By convention, the result with `n_max` 0 and 1 will be the same.

# Returns

* The derivative of the gravitational field w.r.t. the radius (`∂U/∂r`).
* The derivative of the gravitational field w.r.t. the latitude (`∂U/∂ϕ`).
* The derivative of the gravitational field w.r.t. the longitude (`∂U/∂λ`).

# Remarks

In this case, `ϕ` is the geocentric latitude and `λ` is the geocentric
longitude.

"""
function compute_dU(gm_coefs::GravityModel_Coefs{T}, r::AbstractVector,
                    n_max::Number = -1, m_max::Number = -1) where T<:Number

    # Unpack gravity model coefficients.
    #
    # We must not unpack all the variables because `n_max` will be overwritten.
    @unpack μ, R0, C, S, legendre_norm = gm_coefs

    # Check n_max value.
    if (n_max < 0) || (gm_coefs.n_max < n_max)
        n_max = gm_coefs.n_max
    end

    # Check m_max value.
    if (m_max < 0) || (m_max > n_max)
        m_max = n_max
    end

    # Get the geocentric latitude and longitude
    # =========================================
    r_gc = norm(r)
    ρ_gc = sqrt(r[1]^2 + r[2]^2)
    ϕ_gc = atan(r[3], ρ_gc)
    λ_gc = atan(r[2], r[1])

    # Auxiliary variables
    # ===================

    # Sine and cosine of the geocentric longitude.
    #
    # This values were be used in the algorithm to decrease the computational
    # burden.
    sin_λ,  cos_λ  = sincos(1λ_gc)
    sin_2λ, cos_2λ = sincos(2λ_gc)

    # First derivative of the non-spherical portion of the grav. field
    # ================================================================
    ∂Ur = T(1)  # Derivative w.r.t. the radius.
    ∂Uϕ = T(0)  # Derivative w.r.t. the geocentric latitude.
    ∂Uλ = T(0)  # Derivative w.r.t. the geocentric longitude.

    # Compute the associated Legendre functions `P_n,m[sin(ϕ_gc)]` with the
    # required normalization and its time derivative.
    #
    # Notice that cos(ϕ_gc-pi/2) = sin(ϕ_gc).
    dP, P = dlegendre(Val{legendre_norm}, ϕ_gc-pi/2, n_max, n_max, false)

    # Auxiliary variables.
    rn_fact = R0/r_gc
    rn      = rn_fact

    # Compute the derivatives.
    @inbounds for n = 2:n_max
        rn *= rn_fact # -> rn = (R0/r_gc)^n

        aux_∂Ur = T(0)
        aux_∂Uϕ = T(0)
        aux_∂Uλ = T(0)

        # Sine and cosine with m = 0.
        #
        # This values will be used to update recursively `sin(m*λ_gc)` and
        # `cos(m*λ_gc`), reducing the computational burden.
        sin_mλ   = T(0)      # sin( 0*λ_gc)
        sin_m_1λ = -sin_λ    # sin(-1*λ_gc)
        sin_m_2λ = -sin_2λ   # sin(-2*λ_gc)
        cos_mλ   = T(1)      # cos( 0*λ_gc)
        cos_m_1λ = cos_λ     # cos(-1*λ_gc)
        cos_m_2λ = cos_2λ    # cos(-2*λ_gc)

        for m = 0:n
            # Compute recursively `sin(m*λ_gc)` and `cos(m*λ_gc)`.
            sin_mλ = 2cos_λ*sin_m_1λ-sin_m_2λ
            cos_mλ = 2cos_λ*cos_m_1λ-cos_m_2λ

            CcSs_nm = C[n+1,m+1]*cos_mλ + S[n+1,m+1]*sin_mλ

            aux_∂Ur +=   P[n+1,m+1]*CcSs_nm
            aux_∂Uϕ +=  dP[n+1,m+1]*CcSs_nm
            aux_∂Uλ += m*P[n+1,m+1]*(S[n+1,m+1]*cos_mλ - C[n+1,m+1]*sin_mλ)

            # Check if we reached the maximum desired order.
            m >= m_max && break

            # Update the values for the next step.
            sin_m_2λ = sin_m_1λ
            sin_m_1λ = sin_mλ
            cos_m_2λ = cos_m_1λ
            cos_m_1λ = cos_mλ
        end

        ∂Ur += rn*(n+1)*aux_∂Ur
        ∂Uϕ += rn*aux_∂Uϕ
        ∂Uλ += rn*aux_∂Uλ
    end

    ∂Ur *= -μ/r_gc^2
    ∂Uϕ *= +μ/r_gc
    ∂Uλ *= +μ/r_gc

    ∂Ur, ∂Uϕ, ∂Uλ
end

"""
    function compute_g(gm_coefs::GravityModel_Coefs{T}, r::AbstractVector, n_max::Number = -1, m_max::Number = -1) where T<:Number

Compute the gravitational acceleration (ITRF) \\[m/s²] at position `r` \\[m]
(ITRF) using the coefficients `gm_coefs` (see `GravityModel_Coefs`). The maximum
degree that will be used while computing the spherical harmonics will be
`n_max` and the maximum order it `m_max`.

If `n_max` is negative, then the maximum available degree will be used. If
`n_max` is omitted, then it defaults to 0.

If `m_max` is negative or if it is greater than `n_max`, then it will be set to
`n_max`. If `m_max` is omitted, then it defaults to 0.

!!! info

    By convention, the result with `n_max` 0 and 1 will be the same.

# Remarks

Notice that this function computes the **gravitational acceleration**. Hence,
the acceleration due to Earth rotation rate **is not** included.

"""
function compute_g(gm_coefs::GravityModel_Coefs{T}, r::AbstractVector,
                   n_max::Number = -1, m_max::Number = -1) where T<:Number

    # Compute the partial derivatives of the gravitational field w.r.t. the
    # spherical coordinates.
    ∂Ur, ∂Uϕ, ∂Uλ = compute_dU(gm_coefs, r, n_max, m_max)

    # Auxiliary variables.
    r_gc = norm(r)
    ρ_gc = sqrt(r[1]^2 + r[2]^2)
    ϕ_gc = atan(r[3], ρ_gc)
    λ_gc = atan(r[2], r[1])

    # Compute the acceleration represented in the ITRF
    # ================================================

    # Compute the partial derivatives in spherical coordinate systems:
    #
    #     ∂U        1      ∂U     1   ∂U
    #    ---- , ---------.---- , ---.----
    #     ∂r     r.cos ϕ   ∂λ     r   ∂ϕ
    #
    # Notice that the singularity is not a problem here. When computing
    # `cos(pi/2)` a very small number will be returned and `∂U/∂λ` is 0. Hence,
    # the 2nd component will be 0.

    a_l = SVector{3,T}(∂Ur,
                       ∂Uλ/(r_gc*cos(ϕ_gc)),
                       ∂Uϕ/r_gc)

    # The vector `a_l` is represented in the local UEN (Up-Earth-North)
    # reference frame. Hence, we need to describe the unitary vectors of this
    # frame in the ECEF reference frame. This can be accomplished by the
    # following rotations matrix.

    Del = angle_to_dcm(0, ϕ_gc, -λ_gc, :XYZ)
    a_ITRF = Del*a_l

    a_ITRF
end

"""
    function compute_U(gm_coefs::GravityModel_Coefs{T}, r::AbstractVector, n_max::Number = -1, m_max::Number = -1) where T<:Number

Compute the gravitational potential [J/kg] at `r` (ITRF) \\[m] using the
coefficients `gm_coefs` (see `GravityModel_Coefs`). The maximum degree that will
be used while computing the spherical harmonics will be `n_max` and the maximum
order is `m_max`.

If `n_max` is negative, then the maximum available degree will be used. If
`n_max` is omitted, then it defaults to 0.

If `m_max` is negative or if it is greater than `n_max`, then it will be set to
`n_max`. If `m_max` is omitted, then it defaults to 0.

!!! info

    By convention, the result with `n_max` 0 and 1 will be the same.

"""
function compute_U(gm_coefs::GravityModel_Coefs{T}, r::AbstractVector,
                   n_max::Number = -1, m_max::Number = -1) where T<:Number

    # Unpack gravity model coefficients.
    #
    # We must not unpack all the variables because `n_max` will be overwritten.
    @unpack μ, R0, C, S, legendre_norm = gm_coefs

    # Check n_max value.
    if (n_max < 0) || (gm_coefs.n_max < n_max)
        n_max = gm_coefs.n_max
    end

    # Check m_max value.
    if (m_max < 0) || (m_max > n_max)
        m_max = n_max
    end

    # Auxiliary variables.
    r_gc = norm(r)
    ρ_gc = sqrt(r[1]^2 + r[2]^2)
    ϕ_gc = atan(r[3], ρ_gc)
    λ_gc = atan(r[2], r[1])

    # Sine and cosine of the geocentric longitude.
    #
    # This values were be used in the algorithm to decrease the computational
    # burden.
    sin_λ,  cos_λ  = sincos(1λ_gc)
    sin_2λ, cos_2λ = sincos(2λ_gc)

    # Consider the zero degree term.
    U = T(1)

    # Compute the associated Legendre functions `P_n,m[sin(ϕ_gc)]` with the
    # required normalization.
    #
    # Notice that cos(ϕ_gc-pi/2) = sin(ϕ_gc).
    P = legendre(Val{legendre_norm}, ϕ_gc-pi/2, n_max, n_max, false)

    # Auxiliary variables.
    rn_fact = R0/r_gc
    rn      = rn_fact

    @inbounds for n = 2:n_max
        aux_U = T(0)

        # Sine and cosine with m = 0.
        #
        # This values will be used to update recursively `sin(m*λ_gc)` and
        # `cos(m*λ_gc`), reducing the computational burden.
        sin_mλ   = T(0)      # sin( 0*λ_gc)
        sin_m_1λ = -sin_λ    # sin(-1*λ_gc)
        sin_m_2λ = -sin_2λ   # sin(-2*λ_gc)
        cos_mλ   = T(1)      # cos( 0*λ_gc)
        cos_m_1λ = cos_λ     # cos(-1*λ_gc)
        cos_m_2λ = cos_2λ    # cos(-2*λ_gc)

        for m = 0:n
            # Compute recursively `sin(m*λ_gc)` and `cos(m*λ_gc)`.
            sin_mλ = 2cos_λ*sin_m_1λ-sin_m_2λ
            cos_mλ = 2cos_λ*cos_m_1λ-cos_m_2λ

            aux_U += P[n+1,m+1]*(C[n+1,m+1]*cos_mλ + S[n+1,m+1]*sin_mλ)

            # Check if we reached the maximum desired order.
            m >= m_max && break

            # Update the values for the next step.
            sin_m_2λ = sin_m_1λ
            sin_m_1λ = sin_mλ
            cos_m_2λ = cos_m_1λ
            cos_m_1λ = cos_mλ
        end

        rn *= rn_fact # -> rn = (R0/r_gc)^n
        U  += rn*aux_U
    end

    U *= μ/r_gc
end

