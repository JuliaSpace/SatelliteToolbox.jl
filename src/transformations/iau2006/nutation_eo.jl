# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Compute the nutation and the Equation of the Origins (EO) as in
#   equinox-based IAU-2006 theory.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] ftp://tai.bipm.org/iers/conv2010/chapter5/tab5.2e.txt
#
#   [3] Wallace, P. T., Capitaine, N (2006). Precession-nutation procedures
#       consistent with IAU 2006 resolutions. Astronomy & Astrophysics.
#
#   [4] Capitaine, N., Wallace, P. T (2006). High precision methods for locating
#       the celestial intermediate pole and origin. Astronomy & Astrophysics.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export mean_obliquity_iau2006, nutation_eo_iau2006

"""
    mean_obliquity_iau2006(jd_tt::Number)

Compute the mean obliquity of the ecliptic [rad] using the equinox-based
IAU-2006 theory in the Julian day `jd_tt` [Terrestiral Time].

The algorithm was obtained in **[1]**.

# Reference

- **[1]**: Wallace, P. T., Capitaine, N (2006). Precession-nutation procedures
    consistent with IAU 2006 resolutions. Astronomy & Astrophysics.
"""
function mean_obliquity_iau2006(jd_tt::Number)
    # Compute the Julian Centuries from `jd_tt`.
    t_tt = (jd_tt - JD_J2000) / 36525

    # Auxiliary variables
    # ===================

    a2d = 1 / 3600
    d2r = π / 180
    a2r = a2d * d2r

    # Mean obliquity of the ecliptic
    # ==============================

    # Compute the mean obliquity of the ecliptic [s].

    # NOTE: This equation is wrong in [1](p. 216, eq. 3-68)!
    # The one used here was obtained in [3].
    mϵ_2000 = @evalpoly(
        t_tt,
        +84381.406,
        -46.836769,
        -0.0001831,
        +0.00200340,
        -0.000000576,
        -0.0000000434
    )

    # Reduce to the interval [0, 2π].
    mϵ_2000 = mod(mϵ_2000 * a2r, 2π)

    return mϵ_2000
end

"""
    nutation_eo_iau2006(jd_tt::Number)

Compute the nutation parameters and the Equation of Origins (EO) at the Julian
Day `jd_tt` [TT] using the equinox-based 2006 IAU Theory of Nutation. Notice
that one can provide corrections for the nutation in obliquity (`δΔϵ_2000`)
[rad] and in longitude (`δΔψ_2000`) [rad] that are usually obtained from IERS
EOP Data (see [`get_iers_eop`](@ref)).

# Returns

- The mean obliquity of the ecliptic [rad].
- The nutation in obliquity of the ecliptic [rad].
- The nutation in longitude [rad].
- The Equation of Origins (EO) [rad].
"""
function nutation_eo_iau2006(
    jd_tt::Number,
    δΔϵ_2000::Number = 0,
    δΔΨ_2000::Number = 0
)
    # Compute the Julian Centuries from `jd_tt`.
    t_tt = (jd_tt - JD_J2000) / 36525

    # Auxiliary variables
    # ===================

    a2d = 1 / 3600
    d2r = π / 180
    a2r = a2d * d2r
    r2d = 180 / π
    d2a = 3600
    r2a = r2d * d2a

    # Fundamental arguments
    # =====================

    # Luni-solar part.
    M_s, M_m, u_Mm, D_s, Ω_m = luni_solar_args_iau2006(jd_tt)

    # Planetary part.
    λ_M☿, λ_M♀, λ_Me, λ_M♂, λ_M♃, λ_M♄, λ_M⛢, λ_M♆, p_λ = planetary_args_iau2006(jd_tt)

    # Mean obliquity of the ecliptic
    # ==============================

    mϵ_2000 = mean_obliquity_iau2006(jd_tt)

    # Nutation in the obliquity
    # ==========================================================================

    Δϵ_2000 = _iau2006_sum(
        (
            COEFS_IAU2006_NUT_OBL0,
            COEFS_IAU2006_NUT_OBL1
        ),
        t_tt,
        M_s,
        M_m,
        u_Mm,
        D_s,
        Ω_m,
        λ_M☿,
        λ_M♀,
        λ_Me,
        λ_M♂,
        λ_M♃,
        λ_M♄,
        λ_M⛢,
        λ_M♆,
        p_λ
    )

    # Apply the correction given EOP data.
    Δϵ_2000 += δΔϵ_2000 * r2a

    # Nutation in longitude
    # ==========================================================================

    ΔΨ_2000 = _iau2006_sum(
        (
            COEFS_IAU2006_NUT_LON0,
            COEFS_IAU2006_NUT_LON1
        ),
        t_tt,
        M_s,
        M_m,
        u_Mm,
        D_s,
        Ω_m,
        λ_M☿,
        λ_M♀,
        λ_Me,
        λ_M♂,
        λ_M♃,
        λ_M♄,
        λ_M⛢,
        λ_M♆,
        p_λ
    )

    # Apply the correction given EOP data.
    ΔΨ_2000 += δΔΨ_2000 * r2a

    # Equation of Origins (EO)
    # ==========================================================================

    ΔEO = _iau2006_sum(
        (
            COEFS_IAU2006_EO_0,
            COEFS_IAU2006_EO_1
        ),
        t_tt,
        M_s,
        M_m,
        u_Mm,
        D_s,
        Ω_m,
        λ_M☿,
        λ_M♀,
        λ_Me,
        λ_M♂,
        λ_M♃,
        λ_M♄,
        λ_M⛢,
        λ_M♆,
        p_λ
    )

    # NOTE: This equation is wrong in [1](p. 216, eq. 3-66)! The coefficients
    # were obtained from [2] and [4](eq. 69).
    EO = @evalpoly(
        t_tt,
        -0.014_506,
        -4612.156_534,
        -1.391_581_7,
        +0.000_000_44,
        +0.000_029_956,
        +3.68e-8
    ) - ΔΨ_2000 * cos(mϵ_2000) - ΔEO

    # Convert to [rad].
    EO *= a2r
    Δϵ_2000 *= a2r
    ΔΨ_2000 *= a2r

    return mϵ_2000, Δϵ_2000, ΔΨ_2000, EO
end
