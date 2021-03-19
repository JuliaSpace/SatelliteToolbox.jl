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
#   [3] Wallace, P. T.; Capitaine, N (2006). Precession-nutation procedures
#       consistent with IAU 2006 resolutions. Astronomy & Astrophysics.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export nutation_eo_iau2006

"""
    nutation_eo_iau2006(JD_TT::Number)

Compute the nutation parameters and the Equation of Origins (EO) at the Julian
Day `JD_TT` [TT] using the equinox-based 2006 IAU Theory of Nutation.

# Returns

* The mean obliquity of the ecliptic [rad].
* The nutation in obliquity of the ecliptic [rad].
* The nutation in longitude [rad].
* The Equation of Origins (EO) [rad].

"""
function nutation_eo_iau2006(JD_TT::Number)
    # Compute the Julian Centuries from `JD_TT`.
    T_TT = (JD_TT - JD_J2000)/36525

    # Auxiliary variables
    # ===================

    a2d = 1/3600
    d2r = π/180
    a2r = a2d*d2r
    J2d = -2.7774e-6

    # Fundamental arguments
    # =====================

    # Luni-solar part.
    M_s, M_m, u_Mm, D_s, Ω_m = luni_solar_args_iau2006(JD_TT)

    # Planetary part.
    λ_M☿, λ_M♀, λ_Me, λ_M♂, λ_M♃, λ_M♄, λ_M⛢, λ_M♆, p_λ = planetary_args_iau2006(JD_TT)

    # Mean obliquity of the ecliptic
    # ==============================

    # Compute the mean obliquity of the ecliptic [s].

    # TODO: This equation is wrong in [1]! The one used here was obtained in
    # [3].
    mϵ_2000 = @evalpoly(T_TT, +23.439_279,
                               -0.013_010_2,
                               -5.086e-8,
                               +5.565e-7,
                               -1.6e-10,
                               -1.21e-11)

    # Reduce to the interval [0, 2π]°.
    mϵ_2000 = mod(mϵ_2000, 360)*d2r

    # Nutation in the obliquity
    # ==========================================================================

    Δϵ_2000 = _iau2006_sum((coefs_iau2006_nut_obl0, coefs_iau2006_nut_obl1),
                           T_TT, M_s, M_m, u_Mm, D_s, Ω_m, λ_M☿, λ_M♀, λ_Me,
                           λ_M♂, λ_M♃, λ_M♄, λ_M⛢, λ_M♆, p_λ)

    # Nutation in longitude
    # ==========================================================================

    ΔΨ_2000 = _iau2006_sum((coefs_iau2006_nut_lon0, coefs_iau2006_nut_lon1),
                           T_TT, M_s, M_m, u_Mm, D_s, Ω_m, λ_M☿, λ_M♀, λ_Me,
                           λ_M♂, λ_M♃, λ_M♄, λ_M⛢, λ_M♆, p_λ)

    # Equation of Origins (EO)
    # ==========================================================================

    ΔEO = _iau2006_sum((coefs_iau2006_EO_0, coefs_iau2006_EO_1),
                       T_TT, M_s, M_m, u_Mm, D_s, Ω_m, λ_M☿, λ_M♀, λ_Me, λ_M♂,
                       λ_M♃, λ_M♄, λ_M⛢, λ_M♆, p_λ)

    # TODO: This equation is wrong in [1]! The coefficients were obtained from
    # [2].
    EO = @evalpoly(T_TT,   -0.014_506,
                        -4612.156_534,
                           -1.391_581_7,
                           +0.000_000_44,
                           +0.000_029_956,
                           +3.68e-8) - ΔΨ_2000 * cos(mϵ_2000) - ΔEO

    # Convert to [rad].
    EO *= a2r
    Δϵ_2000 *= a2r
    ΔΨ_2000 *= a2r

    return mϵ_2000, Δϵ_2000, ΔΨ_2000, EO
end
