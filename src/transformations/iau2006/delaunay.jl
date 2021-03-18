# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Compute the Delaunay arguments as in IAU 2006 theory.
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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""
    delaunay_arguments_iau2006(JD_TT::Number)

Compute the Delaunay arguments and they account for luni-solar nutation in
the Julian Day `JD_TT` [TT] [1, p. 210].

The returned values are in [rad].

"""
function delaunay_arguments_iau2006(JD_TT::Number)
    # Compute the Julian Centuries from `JD_TT`.
    T_TT = (JD_TT - JD_J2000)/36525

    # Auxiliary variables
    # ===================

    a2d = 1/3600
    d2r = π/180

    # Delaunay arguments of the Sun and Moon
    # ======================================
    #
    # Evaluate the Delaunay arguments associated with the Moon and the Sun in
    # [arcsec] [1, p. 210].

    M_m  = @evalpoly(T_TT, +485868.249036, +1717915923.2178, +31.8792, +0.051635, -0.00024470)
    M_s  = @evalpoly(T_TT, +1287104.79305, +129596581.0481,   -0.5532, +0.000136, -0.00001149)
    u_Mm = @evalpoly(T_TT, +335779.526232, +1739527262.8478, -12.7512, -0.001037, +0.00000417)
    D_s  = @evalpoly(T_TT, +1072260.70369, +1602961601.2090,  -6.3706, +0.006593, -0.00003169)
    Ω_m  = @evalpoly(T_TT, +450160.398036,    -6962890.5431,  +7.4722, +0.007702, -0.00005939)

    # Convert to the interval [0,2π].
    M_s  = mod(M_s*a2d,  360)*d2r
    M_m  = mod(M_m*a2d,  360)*d2r
    u_Mm = mod(u_Mm*a2d, 360)*d2r
    D_s  = mod(D_s*a2d,  360)*d2r
    Ω_m  = mod(Ω_m*a2d,  360)*d2r

    return M_s, M_m, u_Mm, D_s, Ω_m
end
