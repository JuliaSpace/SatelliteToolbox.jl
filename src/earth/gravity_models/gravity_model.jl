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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-06-15: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version. Adapted from the old `egm.jl` file.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export compute_dU, compute_U, compute_g, parse_gfc

################################################################################
#                               Public Functions
################################################################################

"""
    function compute_dU(gm_coefs::GravityModel_Coefs{T}, r::AbstractVector, n_max::Number = 0) where T<:Number

Compute the derivatives w.r.t. the spherical coordinates of the gravitational
field (`∂U/∂r`, `∂U/∂ϕ`, `∂U/∂λ`) defined by the coefficients `gm_coefs` at the
position `r` in ITRF.

Notice that the zero degree term (`-μ/r²`) is considered in `∂U/∂r`.

# Args

* `gm_coefs`: Gravity model coefficients (see `GravityModel_Coefs`).
* `r`: Position in ITRF in which the gravity will be computed [m].
* `n_max`: (OPTIONAL) Maximum degree when computing the spherical harmonics
           (**Default**: 0, which uses the maximum degree).

# Returns

* The derivative of the gravitational field w.r.t. the radius (`∂U/∂r`).
* The derivative of the gravitational field w.r.t. the latitude (`∂U/∂ϕ`).
* The derivative of the gravitational field w.r.t. the longitude (`∂U/∂λ`).

# Remarks

In this case, `ϕ` is the geocentric latitude and `λ` is the geocentric
longitude.

"""
function compute_dU(gm_coefs::GravityModel_Coefs{T},
                    r::AbstractVector,
                    n_max::Number = 0) where T<:Number
    # Unpack gravity model coefficients.
    @unpack_GravityModel_Coefs gm_coefs

    # Check if n_max is bigger than the maximum degree available.
    if (n_max == 0) || (size(C,1)-1 < n_max)
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

    # Sine and cosine of the geocentric longitude.
    #
    # This values were be used in the algorithm to decrease the computational
    # burden.
    sin_λ   = sin(1λ_gc)
    cos_λ   = cos(1λ_gc)
    sin_2λ  = sin(2λ_gc)
    cos_2λ  = cos(2λ_gc)

    # First derivative of the non-spherical portion of the grav. field.
    # =================================================================
    ∂Ur = T(1)  # Derivative w.r.t. the radius.
    ∂Uϕ = T(0)  # Derivative w.r.t. the geocentric latitude.
    ∂Uλ = T(0)  # Derivative w.r.t. the geocentric longitude.

    # Compute the associated Legendre functions `P_n,m[sin(ϕ_gc)]` with the
    # required normalization and its time derivative.
    #
    # Notice that cos(ϕ_gc-pi/2) = sin(ϕ_gc).
    P  =  legendre(Val{legendre_norm}, ϕ_gc-pi/2, n_max, false)
    dP = dlegendre(Val{legendre_norm}, ϕ_gc-pi/2, P    , false)

    # Compute the derivatives.
    for n = 2:n_max
        rn = (R0/r_gc)^n

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
            sin_mλ = 2*cos_λ*sin_m_1λ-sin_m_2λ
            cos_mλ = 2*cos_λ*cos_m_1λ-cos_m_2λ

            CcSs_nm = C[n+1,m+1]*cos_mλ + S[n+1,m+1]*sin_mλ

            aux_∂Ur +=   P[n+1,m+1]*CcSs_nm
            aux_∂Uϕ +=  dP[n+1,m+1]*CcSs_nm
            aux_∂Uλ += m*P[n+1,m+1]*(S[n+1,m+1]*cos_mλ - C[n+1,m+1]*sin_mλ)

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
    function compute_g(gm_coefs::GravityModel_Coefs{T}, r::AbstractVector [, n_max::Number]) where T<:Number

Compute the gravity acceleration at position `r` (ITRF) using the coefficients
`gm_coefs`. The maximum degree that will be used while computing the spherical
harmonics will be `n_max`.

# Args

* `gm_coefs`: Gravity model coefficients (see `GravityModel_Coefs`).
* `r`: Position in ITRF in which the gravity will be computed [m].
* `n_max`: (OPTIONAL) Maximum degree when computing the spherical harmonics
           (**Default**: 0, which uses the maximum degree).

# Returns

A vector with the gravity acceleration represented in ITRF (Earth Centered,
Earth fixed frame).

"""
function compute_g(gm_coefs::GravityModel_Coefs{T},
                   r::AbstractVector,
                   n_max::Number = 0) where T<:Number

    # Compute the partial derivatives of the gravitational field w.r.t. the
    # spherical coordinates.
    ∂Ur, ∂Uϕ, ∂Uλ = compute_dU(gm_coefs, r, n_max)

    # Auxiliary variables.
    r_gc     = norm(r)
    ρ_gc     = sqrt(r[1]^2 + r[2]^2)
    ϕ_gc     = atan2(r[3], ρ_gc)
    λ_gc     = atan2(r[2], r[1])

    # Compute the acceleration represented in the ITRF.
    # =================================================

    # Compute the partial derivatives in spherical coordinate systems:
    #
    #     ∂U        1      ∂U     1   ∂U
    #    ---- , ---------.---- , ---.----
    #     ∂r     r.cos ϕ   ∂λ     r   ∂ϕ
    #
    # Notice that the singularity is not a problem here. When computing
    # `cos(pi/2)` a very small number will be returned and `∂U/∂λ` is 0. Hence,
    # the 2nd component will be 0.

    a_l = [∂Ur;
           ∂Uλ/(r_gc*cos(ϕ_gc));
           ∂Uϕ/r_gc]

    # The vector `a_l` is represented in the local UEN (Up-Earth-North)
    # reference frame. Hence, we need to describe the unitary vectors of this
    # frame in the ECEF reference frame. This can be accomplished by the
    # following rotations matrix.

    Del = angle2dcm(0, ϕ_gc, -λ_gc, :XYZ)
    a_ITRF = Del*a_l

    a_ITRF
end

"""
    function compute_U(gm_coefs::GravityModel_Coefs{T}, ϕ_gc::Number, λ_gc::Number, r_gc::Number [, n_max::Number]) where T<:Number

Compute the gravity potential using the coefficients `gm_coefs` at the
geocentric latitude `ϕ_gc`, geocentric longitude `λ_gc`, and geocentric radius
`r_gc`. The maximum degree that will be used while computing the spherical
harmonics will be `n_max`. If `n_max` is less or equal 0, then the maximum
available degree will be used.

# Args

* `gm_coefs`: Gravity model coefficients (see `GravityModel_Coefs`).
* `ϕ_gc`: Geocentric latitude [rad].
* `λ_gc`: Geocentric longitude [rad].
* `r_gc`: Geocentric radius [m].
* `n_max`: (OPTIONAL) Maximum degree when computing the spherical harmonics
           (**Default**: 0, which uses the maximum degree).

# Returns

The gravitational potential at specified location [J/kg].

"""
function compute_U(gm_coefs::GravityModel_Coefs{T},
                   ϕ_gc::Number,
                   λ_gc::Number,
                   r_gc::Number,
                   n_max::Number = 0) where T<:Number
    # Unpack gravity model coefficients.
    @unpack_GravityModel_Coefs gm_coefs

    # Check if n_max is bigger than the maximum degree available.
    if (n_max <= 0) || (size(C,1)-1 < n_max)
        n_max = size(C,1)-1
    end

    # Consider the zero degree term.
    U = T(1)

    # Compute the associated Legendre functions `P_n,m[sin(ϕ_gc)]` with the
    # required normalization.
    #
    # Notice that cos(ϕ_gc-pi/2) = sin(ϕ_gc).
    P = legendre(Val{legendre_norm}, ϕ_gc-pi/2, n_max, false)

    for n = 2:n_max
        aux_U = T(0)

        for m = 0:n
            sin_mλ = sin(m*λ_gc)
            cos_mλ = cos(m*λ_gc)

            aux_U += P[n+1,m+1]*(C[n+1,m+1]*cos_mλ + S[n+1,m+1]*sin_mλ)
        end

        U += (R0/r_gc)^n*aux_U
    end

    U *= μ/r_gc
end

"""
    function parse_gfc(filename::String)

This function parse a `gfc` (*Gravity Field Coefficient*) file `filename`. More
information about this format can be seen in [2].

# Args

* `filename`: File name and path of the `gfc` file.

# Returns

An instance of the structure `GravityModel_Coefs` with the parsed values.

"""
function parse_gfc(filename::String)
    # Open the file and count how many lines the header has.
    file = open(filename, "r")

    header_lines = 0

    # Information to be acquired at the header.
    name  = ""
    μ     = 0.0
    R₀    = 0.0
    n_max = 0

    while !eof(file)
        str = readline(file)
        header_lines += 1

        # Check if we are at the end of the header.
        contains(str, "end_of_head") && break

        # Check if the line contains the required information.
        if contains(str, "product_type")
            aux = split(str)

            (aux[2] != "gravity_field") && error("The gfc file $filename has a wrong format.")

        elseif contains(str, "modelname")
            aux  = split(str)
            name = aux[2]

        elseif contains(str, "earth_gravity_constant")
            aux = split(str)
            μ   = parse(Float64, aux[2])

        elseif contains(str, "radius")
            aux = split(str)
            R₀  = parse(Float64, aux[2])

        elseif contains(str, "max_degree")
            aux   = split(str)
            n_max = parse(Int64, aux[2])
        end
    end

    # Check if all information was obtained from the header.
    if (name == "") || (μ == 0) || (R₀ == 0) || (n_max == 0)
        error("The gfc file $filename has a wrong header.")
    end

    # Close the file.
    close(file)

    # Read and parse the coefficients.
    raw = readdlm(filename; skipstart = header_lines)

    # The file will be parsed using Sparse matrices because it is much easier.
    # Perhaps, there is a method with better performance, but this function is
    # called only once per execution.
    I = map(x->Int(x)+1, raw[:,2])
    J = map(x->Int(x)+1, raw[:,3])

    C = sparse(I, J, convert(Vector{Float64}, raw[:,4]))
    S = sparse(I, J, convert(Vector{Float64}, raw[:,5]))

    # Create and return the gravity model coefficients.
    gm_coefs::GravityModel_Coefs{Float64} =
        GravityModel_Coefs(name          = name,
                           μ             = μ,
                           R0            = R₀,
                           n_max         = n_max,
                           C             = full(C),
                           S             = full(S),
                           legendre_norm = :full)

    gm_coefs
end
