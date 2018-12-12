#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   International Geomagnetic Field Model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
#   [2] https://www.ngdc.noaa.gov/IAGA/vmod/igrf12.f
#   [3] https://www.mathworks.com/matlabcentral/fileexchange/34388-international-geomagnetic-reference-field--igrf--model
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export igrf12
export igrf12syn

################################################################################
#                                  Functions
################################################################################

"""
    function igrf12(date::Number, r::Number, λ::Number, Ω::Number, T; show_warns = true)

**IGRF v12 Model**

Compute the geomagnetic field vector [nT] at the date `date` [Year A.D.] and
position (`r`, `λ`, `Ω`).

The position representation is defined by `T`. If `T` is `Val{:geocentric}`,
then the input must be **geocentric** coordinates:

1. Distance from the Earth center `r` [m];
1. Geocentric latitude `λ` (-π/2, +π/2) \\[rad]; and
2. Geocentric longitude `Ω` (-π, +π) \\[rad].

If `T` is `Val{:geodetic}`, then the input must be **geodetic**
coordinates:

1. Altitude above the reference ellipsoid `h` (WGS-84) \\[m];
2. Geodetic latitude `λ` (-π/2, +π/2) \\[rad]; and
3. Geodetic longitude `Ω` (-π, +π) \\[rad].

If `T` is omitted, then it defaults to `Val{:geocentric}`.

Notice that the output vector will be represented in the same reference system
selected by the parameter `T` (geocentric or geodetic). The Y-axis of the output
reference system always points East. In case of **geocentric coordinates**, the
Z-axis points toward the center of Earth and the X-axis completes a right-handed
coordinate system. In case of **geodetic coordinates**, the X-axis is tangent to
the ellipsoid at the selected location and points toward North, whereas the
Z-axis completes a right-hand coordinate system.

# Keywords

* `show_warns`: Show warnings about the data (**Default** = `true`).

# Remarks

The `date` must be greater or equal to 1900 and less than or equal 2025. Notice
that a warning message is printed for dates greater than 2020.

# Disclaimer

This function is an independent implementation of the IGRF model. It contains a
more readable code than the original one in FORTRAN, because it uses features
available in Julia language.

"""
igrf12(date::Number, r::Number, λ::Number, Ω::Number; show_warns = true) =
    igrf12(date, r, λ, Ω, Val{:geocentric}; show_warns = show_warns)

function igrf12(date::Number,
                r::Number,
                λ::Number,
                Ω::Number,
                ::Type{Val{:geocentric}};
                show_warns::Bool = true)
    # Input verification
    # ==================

    # Check the data, since this model is valid for years between 1900 and 2025.
    ( (date < 1900) || (date > 2025) ) &&
    error("This IGRF version will not work for years outside the interval [1900, 2025).")

    # Check if the latitude and longitude are valid.
    ( (λ < -pi/2) || (λ > pi/2) ) &&
    error("The latitude must be between -π/2 and +π/2 rad.")

    ( (Ω < -pi) || (Ω > pi) ) &&
    error("The longitude must be between -π and +π rad.")

    # Warn the user that for dates after the year 2020 the accuracy maybe
    # reduced.
    show_warns && (date > 2020) &&
    @warn("The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2020.")

    # Input variables conversion
    # ==========================

    # Convert latitude / longitude to colatitude and east-longitude.
    θ = pi/2 - λ
    ϕ = (Ω >= 0) ? Ω : 2pi + Ω

    # The input variable `r` is in [m], but all the algorithm requires it to be
    # in [km].
    r /= 1000

    # Preliminary setup
    # =================

    # Compute the epoch that will be used to obtain the coefficients. This is
    # necessary because the IGRF provides coefficients every 5 years. Between
    # two epochs, those coefficients must be interpolated.
    idx   = (date < 2020) ? floor(Int, (date-1900)*0.2+1) : 24
    epoch = 1900 + (idx-1)*5

    # Compute the fraction of time from the epoch of the coefficient selected by
    # `idx`.
    Δt = date - epoch

    # Compute the maximum spherical harmonic degree for the selected date.
    n_max = (epoch < 1995) ? 10 : 13

    # Compute the Schmidt quasi-normalized associated Legendre functions and
    # their first order derivative, neglecting the phase term.
    P  =  legendre(Val{:schmidt}, θ, n_max, false)
    dP = dlegendre(Val{:schmidt}, θ, P,     false)

    # Parameters and auxiliary variables
    # ==================================

    # Auxiliary variables to select the IGRF coefficients.
    H = H_igrf12
    G = G_igrf12

    # Reference radius [km].
    a = 6371.2

    # Auxiliary variables to decrease the computational burden.
    sin_ϕ,  cos_ϕ  = sincos(1ϕ)
    sin_2ϕ, cos_2ϕ = sincos(2ϕ)
    ratio   = a/r
    fact    = ratio

    # Initialization of variables
    # ===========================

    dVr = 0.0   # Derivative of the Geomagnetic potential w.r.t. r.
    dVθ = 0.0   # Derivative of the Geomagnetic potential w.r.t. θ.
    dVϕ = 0.0   # Derivative of the Geomagnetic potential w.r.t. ϕ.
    ΔG  = 0.0   # Auxiliary variable to interpolate the G coefficients.
    ΔH  = 0.0   # Auxiliary variable to interpolate the H coefficients.
    kg  = 1     # Index to obtain the values of the matrix `G`.
    kh  = 1     # Index to obtain the values of the matrix `H`.

    # Geomagnetic potential
    # =====================

    @inbounds for n = 1:n_max
        aux_dVr = 0.0
        aux_dVθ = 0.0
        aux_dVϕ = 0.0

        # Compute the contributions when `m = 0`.
        # =======================================

        # Get the coefficients in the epoch and interpolate to the desired
        # time.
        Gnm_e0 = G[kg,idx+2]

        if date < 2015
            Gnm_e1 = G[kg,idx+3]
            ΔG     = (Gnm_e1-Gnm_e0)/5
        else
            ΔG     = G[kg,27]
        end

        Gnm  = Gnm_e0 + ΔG*Δt
        kg  += 1

        aux_dVr += -(n+1)/r*Gnm*P[n+1,1]
        aux_dVθ += Gnm*dP[n+1,1]

        # Sine and cosine with m = 1.
        # ===========================
        #
        # This values will be used to update recursively `sin(m*ϕ)` and
        # `cos(m*ϕ)`, reducing the computational burden.
        sin_mϕ   = +sin_ϕ    # sin( 1*λ_gc)
        sin_m_1ϕ = 0.0       # sin( 0*λ_gc)
        sin_m_2ϕ = -sin_ϕ    # sin(-1*λ_gc)
        cos_mϕ   = +cos_ϕ    # cos( 1*λ_gc)
        cos_m_1ϕ = 1.0       # cos( 0*λ_gc)
        cos_m_2ϕ = +cos_ϕ    # cos(-2*λ_gc)

        # Compute the contributions when `m ∈ [1,n]`.
        # ===========================================

        for m = 1:n
            # Compute recursively `sin(m*ϕ)` and `cos(m*ϕ)`.
            sin_mϕ = 2cos_ϕ*sin_m_1ϕ-sin_m_2ϕ
            cos_mϕ = 2cos_ϕ*cos_m_1ϕ-cos_m_2ϕ

            # Compute the coefficients `G_nm` and `H_nm`.
            # ===========================================

            # Get the coefficients in the epoch and interpolate to the
            # desired time.
            Gnm_e0 = G[kg,idx+2]
            Hnm_e0 = H[kh,idx+2]

            if date < 2015
                Gnm_e1 = G[kg,idx+3]
                Hnm_e1 = H[kh,idx+3]
                ΔG     = (Gnm_e1-Gnm_e0)/5
                ΔH     = (Hnm_e1-Hnm_e0)/5
            else
                ΔG     = G[kg,27]
                ΔH     = H[kh,27]
            end

            Gnm    = Gnm_e0 + ΔG*Δt
            Hnm    = Hnm_e0 + ΔH*Δt
            kg    += 1
            kh    += 1

            GcHs_nm = Gnm*cos_mϕ + Hnm*sin_mϕ
            GsHc_nm = Gnm*sin_mϕ - Hnm*cos_mϕ

            # Compute the contributions for `m`.
            # ==================================
            aux_dVr += -(n+1)/r*GcHs_nm*P[n+1,m+1]
            aux_dVθ += GcHs_nm*dP[n+1,m+1]
            aux_dVϕ += (θ == 0) ? -m*GsHc_nm*dP[n+1,m+1] : -m*GsHc_nm*P[n+1,m+1]

            # Update the values for the next step.
            # ====================================

            sin_m_2ϕ = sin_m_1ϕ
            sin_m_1ϕ = sin_mϕ
            cos_m_2ϕ = cos_m_1ϕ
            cos_m_1ϕ = cos_mϕ
        end

        # Perform final computations related to the summation in `n`.
        # ===========================================================

        # fact = (a/r)^(n+1)
        fact    *= ratio

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

    # Compute the Geomagnetic field vector in the geocentric reference frame.
    # =======================================================================
    x = +1/r*dVθ
    y = (θ == 0) ? -1/r*dVϕ : -1/(r*sin(θ))*dVϕ
    z = dVr

    B_gc = SVector{3,Float64}(x,y,z)
end

function igrf12(date::Number,
                h::Number,
                λ::Number,
                Ω::Number,
                ::Type{Val{:geodetic}};
                show_warns = true)

    # TODO: This method has a small error (≈ 0.01 nT) compared with the
    # `igrf12syn`.  However, the result is exactly the same as the MATLAB
    # function in [3]. Hence, this does not seem to be an error in the
    # conversion from geodetic to geocentric coordinates. This is probably
    # caused by a numerical error. Further verification is necessary.

    # Convert the geodetic coordinates to geocentric coordinates.
    (λ_gc, r) = GeodetictoGeocentric(λ, h)

    # Compute the geomagnetic field in geocentric coordinates.
    B_gc = igrf12(date, r, λ_gc, Ω, Val{:geocentric}; show_warns = show_warns)

    # Convert to geodetic coordinates.
    D_gd_gc = angle_to_dcm(λ_gc - λ, 0., 0., :YXZ)
    B_gd    = D_gd_gc*B_gc
end

"""
    function igrf12syn(isv::Int, date::Number, itype::Int, alt::Number, colat::Number, elong::Number; show_warns = true)

This is a Julia implementation of the official IGRF source code, which was
written in Fortran [2]. The input and output variables are exactly the same as
the ones described in the function `igrf12syn` in [2].

# Args

* `isv`: `0` if main-field values are required, `1` if secular variation values
         are required.
* `date`: Year A.D.
* `itype`: `1` if geodetic (spheroid), `2` if geocentric (sphere).
* `alt`: Height above sea level [km] if `itype = 1`, or distance from the center of
         Earth [km] if `itype = 2` (must be > 3485 km).
* `colat`: Colatitude (0 - 180) [˚].
* `elong`: East-Longitude (0 - 360) [˚].

# Keywords

* `show_warns`: Show warnings about the data (**Default** = `true`).

# Returns

* The north component [nT] if `isv = 0`, or [nT/year] if `isv = 1`.
* The east component [nT] if `isv = 0`, or [nT/year] if `isv = 1`.
* The vertical component [nT] if `isv = 0`, or [nT/year] if `isv = 1`.
* The total intensity if `isv = 0`, or rubbish if `isv = 1`.

# Remarks

* The `date` must be greater or equal to 1900 and less than or equal 2025.
Notice that a warning message is printed for dates grated than 2020.

"""
function igrf12syn(isv::Int,
                   date::Number,
                   itype::Int,
                   alt::Number,
                   colat::Number,
                   elong::Number;
                   show_warns = true)

    # Check the data, since this model is valid for years between 1900 and 2025.
    ( (date < 1900) || (date > 2025) ) &&
    error("This IGRF version will not work for years outside the interval [1900, 2025).")

    # Warn the user that for dates after the year 2020 the accuracy maybe reduced.
    show_warns && (date > 2020) &&
    @warn("The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2020.")

    # Declaration of variables
    gh          = gh_igrf12
    fn::Int     = 0
    gn::Int     = 0
    kmx::Int    = 0
    ll::Int     = 0
    nc::Int     = 0
    nmx::Int    = 0
    x::Float64  = 0.0
    y::Float64  = 0.0
    z::Float64  = 0.0
    t::Float64  = 0.0
    tc::Float64 = 0.0
    cl          = zeros(13)
    sl          = zeros(13)
    p           = zeros(105)
    q           = zeros(105)

    if date < 2015
        t  = 0.2*(date - 1900)
        ll = floor(Int,t)
        t  = t - ll

        # SH models before 1995.0 are only to degree 10.
        if date < 1995
            nmx = 10

            # nc = nmx*(nmx+2)
            nc  = 120

            ll  = nc*ll

            # kmx = (nmx+1)*(nmx+2)/2
            kmx = 66
        else
            nmx = 13

            # nc  = nmx*(nmx+2)
            nc  = 195

            ll  = floor(Int,0.2*(date - 1995))

            # 19 is the number of SH models that extend to degree 10.
            ll  = 120*19 + nc*ll

            # kmx = (nmx+1)*(nmx+2)/2
            kmx = 105
        end

        tc = 1 - t

        if isv == 1
            t  = +0.2
            tc = -0.2
        end
    else
        t  = date - 2015
        tc = 1.0

        if isv == 1
            t  = 1.0
            tc = 0.0
        end

        ll  = 3060
        nmx = 13

        # nc  = nmx*(nmx+2)
        nc = 195

        # kmx = (nmx+1)*(nmx+2)/2
        kmx = 105
    end

    r::Float64  = alt
    ct          = cos(colat*pi/180)
    st          = sin(colat*pi/180)
    cl[1]       = cos(elong*pi/180)
    sl[1]       = sin(elong*pi/180)
    cd::Float64 = 1.0
    sd::Float64 = 0.0
    l::Int      = 1
    m::Int      = 1
    n::Int      = 0

    if itype != 2
        # Conversion between geodetic and geocentric coordinates using WGS84
        # spheroid.
        a2    = 40680631.6
        b2    = 40408296.0
        one   = a2*st^2
        two   = b2*ct^2
        three = one + two
        rho   = sqrt(three)
        r     = sqrt(alt*(alt + 2*rho) + (a2*one + b2*two)/three)
        cd    = (alt + rho)/r
        sd    = (a2 - b2)/rho*ct*st/r
        one   = ct
        ct    = ct*cd -  st*sd
        st    = st*cd + one*sd
    end

    ratio = 6371.2/r
    rr    = ratio^2

    # Computation of Schmidt quasi-normal coefficients p and x(=q).
    # =============================================================

    p[1] = 1.0
    p[3] = st
    q[1] = 0.0
    q[3] = ct

    for k = 2:kmx
        # There is no need to check bounds here. The code guarantees that
        # everything is inside the right bounds. This increased the source code
        # performance by 13%.
        @inbounds begin
            if n < m
                m  = 0
                n  = n+1
                rr = rr*ratio
                fn = n
                gn = n-1
            end

            fm = m

            if (m == n)
                if k != 3
                    one   = sqrt(1 - 0.5/fm)
                    j     = k - n - 1
                    p[k]  = one*st*p[j]
                    q[k]  = one*(st*q[j] + ct*p[j])
                    cl[m] = cl[m-1]*cl[1] - sl[m-1]*sl[1]
                    sl[m] = sl[m-1]*cl[1] + cl[m-1]*sl[1]
                end
            else
                gmm   = m^2
                one   = sqrt(fn^2 - gmm)
                two   = sqrt(gn^2 - gmm)/one
                three = (fn + gn)/one
                i     = k - n
                j     = i - n + 1
                p[k]  = three*ct*p[i] - two*p[j]
                q[k]  = three*(ct*q[i] - st*p[i]) - two*q[j]
            end

            # Synthesis of x, y, and z in geocentric coordinates.
            lm  = ll + l
            one = (tc*gh[lm] + t*gh[lm+nc])*rr

            if m != 0
                two   = (tc*gh[lm+1] + t*gh[lm+nc+1])*rr
                three = one*cl[m] + two*sl[m]
                x     = x + three*q[k]
                z     = z - (fn + 1)*three*p[k]

                if st != 0
                    y = y + (one*sl[m] - two*cl[m])*fm*p[k]/st
                else
                    y = y + (one*sl[m] - two*cl[m])*q[k]*ct
                end

                l = l + 2
            else
                x = x + one*q[k]
                z = z - (fn + 1)*one*p[k]
                l = l + 1
            end

            m = m + 1
        end
    end

    # Conversion to coordinate system specified by itype.
    # ===================================================
    one   = x
    x     = x*cd +   z*sd
    z     = z*cd - one*sd
    f     = sqrt(x^2 + y^2 + z^2)

    (x,y,z,f)
end
