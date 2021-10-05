# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   International Geomagnetic Field Model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
#   [2] https://www.ngdc.noaa.gov/IAGA/vmod/igrf12.f
#   [3] https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field--igrf--model
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export igrf, igrfd

################################################################################
#                                  Functions
################################################################################

"""
    igrfd(date::Number, [r,h]::Number, λ::Number, Ω::Number, T[, P, dP]; max_degree = 13, show_warns = true)

**IGRF Model**

*Current version: v13*

Compute the geomagnetic field vector [nT] at the date `date` [Year A.D.] and
position (`r`, `λ`, `Ω`).

The position representation is defined by `T`. If `T` is `Val(:geocentric)`,
then the input must be **geocentric** coordinates:

1. Distance from the Earth center `r` [m];
2. Geocentric latitude `λ` (-90°, +90°); and
3. Geocentric longitude `Ω` (-180°, +180°).

If `T` is `Val(:geodetic)`, then the input must be **geodetic** coordinates:

1 Altitude above the reference ellipsoid `h` (WGS-84) \\[m];
2. Geodetic latitude `λ` (-90°, +90°); and
3. Geodetic longitude `Ω` (-180°, +180°).

If `T` is omitted, then it defaults to `Val(:geocentric)`.

Notice that the output vector will be represented in the same reference system
selected by the parameter `T` (geocentric or geodetic). The Y-axis of the output
reference system always points East. In case of **geocentric coordinates**, the
Z-axis points toward the center of Earth and the X-axis completes a right-handed
coordinate system. In case of **geodetic coordinates**, the X-axis is tangent to
the ellipsoid at the selected location and points toward North, whereas the
Z-axis completes a right-hand coordinate system.

The optional arguments `P` and `dP` must be two matrices with at least
`max_degree + 1 × max_degree + 1` real numbers. If they are present, then they
will be used to store the Legendre coefficients and their derivatives. In this
case, no allocation will be performed when computing the magnetic field. If they
are not present, then 2 allocations will happen to create them.

# Keywords

- `max_degree::Int`: Maximum degree used in the spherical harmonics when
    computing the geomagnetic field. If it is higher than the available number
    of coefficients in the IGRF matrices, then it will be clamped. If it is
    equal of lower than 0, then it will be set to 1. (**Default** = 13)
- `show_warns::Bool`: Show warnings about the data (**Default** = `true`).

# Remarks

The `date` must be greater or equal to 1900 and less than or equal 2030. Notice
that a warning message is printed for dates greater than 2025.

# Disclaimer

This function is an independent implementation of the IGRF model. It contains a
more readable code than the original one in FORTRAN, because it uses features
available in Julia language.

"""
@inline function igrfd(
    date::Number,
    r::Number,
    λ::Number,
    Ω::Number;
    max_degree::Int = 13,
    show_warns::Bool = true
)
    return igrfd(
        date,
        r,
        λ,
        Ω,
        Val(:geocentric);
        max_degree = max_degree,
        show_warns = show_warns
    )
end

@inline function igrfd(
    date::Number,
    r::Number,
    λ::Number,
    Ω::Number,
    P::AbstractMatrix{T},
    dP::AbstractMatrix{T};
    max_degree::Int = 13,
    show_warns::Bool = true
) where T<:AbstractFloat
    return igrfd(
        date,
        r,
        λ,
        Ω,
        Val(:geocentric),
        P,
        dP;
        max_degree = max_degree,
        show_warns = show_warns
    )
end

@inline function igrfd(
    date::Number,
    rh::S1,
    λ::S2,
    Ω::S3,
    R::T;
    max_degree::Int = 13,
    show_warns::Bool = true
) where {
    S1<:Number,
    S2<:Number,
    S3<:Number,
    T<:Union{Val{:geocentric}, Val{:geodetic}}
}
    S = promote_type(float(S1), float(S2), float(S3))
    # Currently, the maximum degree allowed for IGRF is 13.
    max_degree = clamp(max_degree, 1, 13)

    P  = Matrix{S}(undef, max_degree + 1, max_degree + 1)
    dP = Matrix{S}(undef, max_degree + 1, max_degree + 1)

    return igrfd(
        date,
        rh,
        λ,
        Ω,
        R,
        P,
        dP;
        max_degree = max_degree,
        show_warns = show_warns
    )
end

@inline function igrfd(
    date::Number,
    rh::Number,
    λ::Number,
    Ω::Number,
    R::T,
    P::AbstractMatrix{S},
    dP::AbstractMatrix{S};
    max_degree::Int = 13,
    show_warns::Bool = true
) where {T<:Union{Val{:geocentric}, Val{:geodetic}}, S<:AbstractFloat}
    # Check if the latitude and longitude are valid.
    if ((λ < -90) || (λ > 90))
        error("The latitude must be between -90° and +90° rad.")
    end

    if ((Ω < -180) || (Ω > 180))
        error("The longitude must be between -180° and +180° rad.")
    end

    return igrf(
        date,
        rh,
        deg2rad(λ),
        deg2rad(Ω),
        R,
        P,
        dP;
        max_degree = max_degree,
        show_warns = show_warns
    )
end

"""
    igrf(date::Number, [r,h]::Number, λ::Number, Ω::Number, T[, P, dP]; max_degree = 13, show_warns = true)

**IGRF Model**

*Current version: v13*

Compute the geomagnetic field vector [nT] at the date `date` [Year A.D.] and
position (`r`, `λ`, `Ω`).

The position representation is defined by `T`. If `T` is `Val(:geocentric)`,
then the input must be **geocentric** coordinates:

1. Distance from the Earth center `r` [m];
2. Geocentric latitude `λ` (-π/2, +π/2) \\[rad]; and
3. Geocentric longitude `Ω` (-π, +π) \\[rad].

If `T` is `Val(:geodetic)`, then the input must be **geodetic** coordinates:

1. Altitude above the reference ellipsoid `h` (WGS-84) \\[m];
2. Geodetic latitude `λ` (-π/2, +π/2) \\[rad]; and
3. Geodetic longitude `Ω` (-π, +π) \\[rad].

If `T` is omitted, then it defaults to `Val(:geocentric)`.

Notice that the output vector will be represented in the same reference system
selected by the parameter `T` (geocentric or geodetic). The Y-axis of the output
reference system always points East. In case of **geocentric coordinates**, the
Z-axis points toward the center of Earth and the X-axis completes a right-handed
coordinate system. In case of **geodetic coordinates**, the X-axis is tangent to
the ellipsoid at the selected location and points toward North, whereas the
Z-axis completes a right-hand coordinate system.

The optional arguments `P` and `dP` must be two matrices with at least
`max_degree + 1 × max_degree + 1` real numbers. If they are present, then they
will be used to store the Legendre coefficients and their derivatives. In this
case, no allocation will be performed when computing the magnetic field. If they
are not present, then 2 allocations will happen to create them.

# Keywords

- `max_degree::Int`: Maximum degree used in the spherical harmonics when
    computing the geomagnetic field. If it is higher than the available number
    of coefficients in the IGRF matrices, then it will be clamped. If it is
    equal of lower than 0, then it will be set to 1. (**Default** = 13)
- `show_warns::Bool`: Show warnings about the data (**Default** = `true`).

# Remarks

The `date` must be greater or equal to 1900 and less than or equal 2030. Notice
that a warning message is printed for dates greater than 2025.

# Disclaimer

This function is an independent implementation of the IGRF model. It contains a
more readable code than the original one in FORTRAN, because it uses features
available in Julia language.

"""
function igrf(
    date::Number,
    r::Number,
    λ::Number,
    Ω::Number;
    max_degree::Int = 13,
    show_warns::Bool = true
)
    return igrf(date,
        r,
        λ,
        Ω,
        Val(:geocentric);
        max_degree = max_degree,
        show_warns = show_warns
    )
end

function igrf(
    date::Number,
    r::Number,
    λ::Number,
    Ω::Number,
    P::AbstractMatrix{T},
    dP::AbstractMatrix{T};
    max_degree::Int = 13,
    show_warns::Bool = true
) where T<:AbstractFloat
    return igrf(
        date,
        r,
        λ,
        Ω,
        Val(:geocentric),
        P,
        dP;
        max_degree = max_degree,
        show_warns = show_warns
    )
end

function igrf(
    date::Number,
    r::S1,
    λ::S2,
    Ω::S3,
    R::T;
    max_degree::Int = 13,
    show_warns::Bool = true
) where {
    S1<:Number,
    S2<:Number,
    S3<:Number,
    T<:Union{Val{:geocentric}, Val{:geodetic}}
}
    S = promote_type(float(S1), float(S2), float(S3))

    # Currently, the maximum degree allowed for IGRF is 13.
    max_degree = clamp(max_degree, 1, 13)

    P  = Matrix{S}(undef, max_degree + 1, max_degree + 1)
    dP = Matrix{S}(undef, max_degree + 1, max_degree + 1)

    return igrf(date,
        r,
        λ,
        Ω,
        R,
        P,
        dP;
        max_degree = max_degree,
        show_warns = show_warns
    )
end

function igrf(
    date::Number,
    r::Number,
    λ::Number,
    Ω::Number,
    ::Val{:geocentric},
    P::AbstractMatrix{T},
    dP::AbstractMatrix{T};
    max_degree::Int = 13,
    show_warns::Bool = true
) where T<:AbstractFloat

    # Model data
    # ==========

    max_year = 2030    # Maximum year in which we can compute.
    rel_year = 2025    # Maximum year with reliable data.
    dat_year = 2020    # Last year in which real measurements are available.

    # Input verification
    # ==================

    # Check the data, since this model is valid for years between 1900 and
    # `max_year`.
    if ((date < 1900) || (date > max_year))
        error("This IGRF version will not work for years outside the interval [1900, $max_year).")
    end

    # Check if the latitude and longitude are valid.
    if ((λ < -π/2) || (λ > π/2))
        error("The latitude must be between -π/2 and +π/2 rad.")
    end

    if ((Ω < -π) || (Ω > π))
        error("The longitude must be between -π and +π rad.")
    end

    # Warn the user that for dates after the year `rel_year` the accuracy maybe
    # reduced.
    if show_warns && (date > rel_year)
        @warn("The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than $rel_year.")
    end

    # If the `max_degree` is equal or lower than 0, we must clamp it to 1.
    max_degree = max(max_degree, 1)

    # Input variables conversion
    # ==========================

    # Convert latitude / longitude to colatitude and east-longitude.
    θ = T(π / 2) - λ
    ϕ = (Ω >= 0) ? Ω : T(2π) + Ω

    # The input variable `r` is in [m], but all the algorithm requires it to be
    # in [km].
    r /= 1000

    # Preliminary setup
    # =================

    # Compute the epoch that will be used to obtain the coefficients. This is
    # necessary because the IGRF provides coefficients every 5 years. Between
    # two epochs, those coefficients must be interpolated.
    idx   = floor(Int, clamp((date - 1900) * 0.2 + 1, 0, (rel_year - 1900) / 5))
    epoch = 1900 + (idx - 1) * 5

    # Compute the fraction of time from the epoch of the coefficient selected by
    # `idx`.
    Δt = date - epoch

    # Compute the maximum spherical harmonic degree for the selected date.
    n_max = (epoch < 1995) ? 10 : 13

    # Check if the user wants a lower degree.
    n_max = min(max_degree, n_max)

    # Make sure we have the required amount of space in the matrices.
    rows, cols = size(P)

    if (rows < n_max + 1) || (cols < n_max + 1)
        error("Matrix `P` must have at least $(n_max + 1) rows and columns.")
    end

    rows, cols = size(dP)
    if ((rows < n_max + 1) || (cols < n_max + 1))
        error("Matrix `dP` must have at least $(n_max+1) rows and columns.")
    end

    # Compute the Schmidt quasi-normalized associated Legendre functions and
    # their first order derivative, neglecting the phase term.
    legendre!(Val(:schmidt), P, θ, false, n_max, n_max)
    dlegendre!(Val(:schmidt), dP, θ, P, false, n_max, n_max)

    # Parameters and auxiliary variables
    # ==================================

    # Auxiliary variables to select the IGRF coefficients.
    H = H_igrf
    G = G_igrf

    # Get the last index of the matrices G and H, which must be the same.
    ~, gend = size(G)
    ~, hend = size(H)

    (gend != hend) &&
    error("Internal error: the matrices G and H must have the same number of columns.")

    # Reference radius [km].
    a = T(6371.2)

    # Auxiliary variables to decrease the computational burden.
    sin_ϕ,  cos_ϕ  = sincos(1ϕ)
    sin_2ϕ, cos_2ϕ = sincos(2ϕ)

    ratio = a/r
    fact  = ratio

    # Initialization of variables
    # ===========================

    dVr = T(0)  # Derivative of the Geomagnetic potential w.r.t. r.
    dVθ = T(0)  # Derivative of the Geomagnetic potential w.r.t. θ.
    dVϕ = T(0)  # Derivative of the Geomagnetic potential w.r.t. ϕ.
    ΔG  = T(0)  # Auxiliary variable to interpolate the G coefficients.
    ΔH  = T(0)  # Auxiliary variable to interpolate the H coefficients.
    kg  = 1     # Index to obtain the values of the matrix `G`.
    kh  = 1     # Index to obtain the values of the matrix `H`.

    # Geomagnetic potential
    # =====================

    @inbounds for n in 1:n_max
        aux_dVr = T(0)
        aux_dVθ = T(0)
        aux_dVϕ = T(0)

        # Compute the contributions when `m = 0`
        # ======================================

        # Get the coefficients in the epoch and interpolate to the desired
        # time.
        Gnm_e0 = T(G[kg,idx+2])

        if date < dat_year
            Gnm_e1 = T(G[kg,idx+3])
            ΔG     = (Gnm_e1 - Gnm_e0) / 5
        else
            ΔG     = T(G[kg,gend])
        end

        Gnm  = Gnm_e0 + ΔG * Δt
        kg  += 1

        aux_dVr += -(n + 1) / r * Gnm * P[n+1,1]
        aux_dVθ += Gnm * dP[n+1,1]

        # Sine and cosine with m = 1
        # ==========================
        #
        # This values will be used to update recursively `sin(m*ϕ)` and
        # `cos(m*ϕ)`, reducing the computational burden.
        sin_mϕ   = +sin_ϕ    # sin( 1*λ_gc)
        sin_m_1ϕ = T(0)      # sin( 0*λ_gc)
        sin_m_2ϕ = -sin_ϕ    # sin(-1*λ_gc)
        cos_mϕ   = +cos_ϕ    # cos( 1*λ_gc)
        cos_m_1ϕ = T(1)      # cos( 0*λ_gc)
        cos_m_2ϕ = +cos_ϕ    # cos(-2*λ_gc)

        # Other auxiliary variables that depend only on `n`
        # =================================================

        fact_dVr = T(n + 1) / T(r)

        # Compute the contributions when `m ∈ [1,n]`
        # ==========================================

        for m = 1:n
            # Compute recursively `sin(m*ϕ)` and `cos(m*ϕ)`.
            sin_mϕ = 2cos_ϕ * sin_m_1ϕ - sin_m_2ϕ
            cos_mϕ = 2cos_ϕ * cos_m_1ϕ - cos_m_2ϕ

            # Compute the coefficients `G_nm` and `H_nm`
            # ==========================================

            # Get the coefficients in the epoch and interpolate to the
            # desired time.
            Gnm_e0 = T(G[kg,idx+2])
            Hnm_e0 = T(H[kh,idx+2])

            if date < dat_year
                Gnm_e1 = T(G[kg,idx+3])
                Hnm_e1 = T(H[kh,idx+3])
                ΔG     = (Gnm_e1 - Gnm_e0) / 5
                ΔH     = (Hnm_e1 - Hnm_e0) / 5
            else
                ΔG = T(G[kg,gend])
                ΔH = T(H[kh,hend])
            end

            Gnm  = Gnm_e0 + ΔG * Δt
            Hnm  = Hnm_e0 + ΔH * Δt
            kg  += 1
            kh  += 1

            GcHs_nm = Gnm * cos_mϕ + Hnm * sin_mϕ
            GsHc_nm = Gnm * sin_mϕ - Hnm * cos_mϕ

            # Compute the contributions for `m`
            # =================================

            aux_dVr += -fact_dVr * GcHs_nm * P[n+1,m+1]
            aux_dVθ += GcHs_nm * dP[n+1,m+1]
            aux_dVϕ += (θ == 0) ?
                -m * GsHc_nm * dP[n+1,m+1] :
                -m * GsHc_nm * P[n+1,m+1]

            # Update the values for the next step
            # ===================================

            sin_m_2ϕ = sin_m_1ϕ
            sin_m_1ϕ = sin_mϕ
            cos_m_2ϕ = cos_m_1ϕ
            cos_m_1ϕ = cos_mϕ
        end

        # Perform final computations related to the summation in `n`
        # ==========================================================

        # fact = (a/r)^(n+1)
        fact *= ratio

        # aux_<> *= (a/r)^(n+1)
        aux_dVr *= fact
        aux_dVϕ *= fact
        aux_dVθ *= fact

        dVr += aux_dVr
        dVϕ += aux_dVϕ
        dVθ += aux_dVθ
    end

    dVr *= a
    dVϕ *= a
    dVθ *= a

    # Compute the Geomagnetic field vector in the geocentric reference frame
    # ======================================================================

    x = +dVθ / r
    y = (θ == 0) ? -dVϕ / r : -dVϕ / (r * sin(θ))
    z = dVr

    B_gc = SVector{3, T}(x, y, z)

    return B_gc
end

function igrf(
    date::Number,
    h::Number,
    λ::Number,
    Ω::Number,
    ::Val{:geodetic},
    P::AbstractMatrix{T},
    dP::AbstractMatrix{T};
    max_degree::Int = 13,
    show_warns::Bool = true
) where T<:AbstractFloat

    # TODO: This method has a small error (≈ 0.01 nT) compared with the
    # `igrf12syn`. However, the result is exactly the same as the MATLAB
    # function in [3]. Hence, this does not seem to be an error in the
    # conversion from geodetic to geocentric coordinates. This is probably
    # caused by a numerical error. Further verification is necessary.

    # Convert the geodetic coordinates to geocentric coordinates.
    λ_gc, r = geodetic_to_geocentric(λ, h)

    # Compute the geomagnetic field in geocentric coordinates.
    B_gc = igrf(
        date,
        r,
        λ_gc,
        Ω,
        Val(:geocentric),
        P,
        dP;
        max_degree = max_degree,
        show_warns = show_warns
    )

    # Convert to geodetic coordinates.
    D_gd_gc = angle_to_dcm(λ_gc - λ, :Y)
    B_gd    = D_gd_gc * B_gc

    return B_gd
end
