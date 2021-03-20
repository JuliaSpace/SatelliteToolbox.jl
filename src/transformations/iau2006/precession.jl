# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Compute the precession as in equinox-based IAU-2006 theory.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] Wallace, P. T., Capitaine, N (2006). Precession-nutation procedures
#       consistent with IAU 2006 resolutions. Astronomy & Astrophysics.
#
#   [3] IERS (2010). Transformation between the International Terrestrial
#       Reference System and the Geocentric Celestial Reference System. IERS
#       Technical Note No. 36, Chapter 5.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""
    precession_iau2006(JD_TT::Number)

Compute the precession angles [rad] according to equinox-based IAU-2006 theory
in the Julia day `JD_TT` [Terrestrial Time].

This algorithm was obtained from [3, p. 49].

"""
function precession_iau2006(JD_TT::Number)
    # Compute the Julian Centuries from `JD_TT`.
    T_TT = (JD_TT - JD_J2000)/36525

    # Auxiliary variables
    # ===================

    a2d = 1/3600
    d2r = π/180
    a2r = a2d*d2r

    Ψ_a = @evalpoly(T_TT,     0,
                          +5038.481_507,
                             -1.079_006_9,
                             -0.001_140_45,
                             +1.328_51e-4,
                             -9.51e-8)
    Ψ_a = mod(Ψ_a * a2r, 2π)

    ω_a = @evalpoly(T_TT, +84381.406,
                              -0.025_754,
                              +0.051_262_3,
                              -0.007_725_03,
                              -4.67e-7,
                              +3.337e-7)
    ω_a = mod(ω_a * a2r, 2π)

    χ_a = @evalpoly(T_TT,   0,
                          +10.556_403,
                           -2.381_429_2,
                           -0.001_211_97,
                           +1.706_63e-4,
                           -5.60e-8)
    χ_a = mod(χ_a * a2r, 2π)

    return Ψ_a, ω_a, χ_a
end
