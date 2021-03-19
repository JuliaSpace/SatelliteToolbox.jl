# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions to compute the Celestial Intermediate Origin (CIO).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] Vallado, D. A (06-Feb-2018). Consolidated Errata of Fundamentals of
#       Astrodynamics and Applications 4th Ed.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export cio_iau2006

"""
    cio_iau2006(JD_TT::Number)

Compute the coordinates `X` and `Y` of the Celestial Intermediate Pole (CIP)
with respect to the Geocentric Celestial Reference Frame (GCRF), and the CIO
locator `s`. The algorithm is based on the IAU-2006 theory.

The CIO locator `s` provides the position of the CIO on the Equator of the CIP
corresponding to the kinematical definition of the non-rotation origin in the
GCRS when the CIP is moving with respect to the GCRS between the reference epoch
and the epoch due to precession and nutation [1, p. 214].

# Returns

* The coordinate `X` of the CIP w.r.t. the GCRF.
* The coordinate `Y` of the CIP w.r.t. the GCRF.
* The CIO locator `s`.

"""
function cio_iau2006(JD_TT::Number)
    # Compute the Julian Centuries from `JD_TT`.
    T_TT = (JD_TT - JD_J2000)/36525

    # Auxiliary variables
    # ===================

    a2d = 1/3600
    d2r = π/180
    a2r = a2d*d2r

    # Fundamental arguments
    # =====================

    # Luni-solar part.
    M_s, M_m, u_Mm, D_s, Ω_m = luni_solar_args_iau2006(JD_TT)

    # Planetary part.
    λ_M☿, λ_M♀, λ_Me, λ_M♂, λ_M♃, λ_M♄, λ_M⛢, λ_M♆, p_λ = planetary_args_iau2006(JD_TT)

    # X position of the CIP
    # ==========================================================================

    ΔX = _iau2006_sum((coefs_iau2006_CIP_x0, coefs_iau2006_CIP_x1,
                       coefs_iau2006_CIP_x2, coefs_iau2006_CIP_x3,
                       coefs_iau2006_CIP_x4),
                       T_TT, M_s, M_m, u_Mm, D_s, Ω_m, λ_M☿, λ_M♀, λ_Me, λ_M♂,
                       λ_M♃, λ_M♄, λ_M⛢, λ_M♆, p_λ)

    X = @evalpoly(T_TT,   -0.016_617,
                       +2004.191_898,
                          -0.429_782_9,
                          -0.198_618_34,
                          +0.000_007_578,
                          +0.000_005_928_5) + ΔX

    # Convert to [rad].
    X *= a2r

    # Y position of the CIP
    # ==========================================================================

    ΔY = _iau2006_sum((coefs_iau2006_CIP_y0, coefs_iau2006_CIP_y1,
                       coefs_iau2006_CIP_y2, coefs_iau2006_CIP_y3,
                       coefs_iau2006_CIP_y4),
                       T_TT, M_s, M_m, u_Mm, D_s, Ω_m, λ_M☿, λ_M♀, λ_Me, λ_M♂,
                       λ_M♃, λ_M♄, λ_M⛢, λ_M♆, p_λ)

    Y = @evalpoly(T_TT, -0.006_951,
                        -0.025_896,
                       -22.407_274_7,
                        +0.001_900_59,
                        +0.001_112_526,
                        +0.000_000_135_8) + ΔY

    # Convert to [rad].
    Y *= a2r

    # Parameter `s` (CIO locator)
    # ==========================================================================
    #
    # The value `s` provides the position of the CIO on the Equator of the CIP
    # corresponding to the kinematical definition of the non-rotation origin in
    # the GCRS when the CIP is moving with respect to the GCRS between the
    # reference epoch and the epoch due to precession and nutation [1, p. 214].

    Δs = _iau2006_sum((coefs_iau2006_CIO_s0, coefs_iau2006_CIO_s1,
                       coefs_iau2006_CIO_s2, coefs_iau2006_CIO_s3,
                       coefs_iau2006_CIO_s4),
                       T_TT, M_s, M_m, u_Mm, D_s, Ω_m, λ_M☿, λ_M♀, λ_Me, λ_M♂,
                       λ_M♃, λ_M♄, λ_M⛢, λ_M♆, p_λ)

    # We must convert this term to [arcsec] to match the units.
    #                      ||
    #                   |------|
    s = @evalpoly(T_TT, (-X*Y/2)/a2r + 0.000_094,
                                     +0.003_808_65,
                                     -0.000_122_68,
                                     -0.072_574_11,
                                     +0.000_027_98,
                                     +0.000_015_65) + Δs

    # Convert to [rad].
    s *= a2r

    return X, Y, s
end
