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

    # Delaunay arguments of the Sun and Moon
    # ======================================

    M_s, M_m, u_Mm, D_s, Ω_m = delaunay_arguments_iau2006(JD_TT)

    # Planetary effects of the nutation and obliquity of the ecliptic
    # ===============================================================

    λ_M☿, λ_M♀, λ_Me, λ_M♂, λ_M♃, λ_M♄, λ_M⛢, λ_M♆, p_λ =
        planetary_effects_args_iau2006(JD_TT)

    # X position of the CIP
    # ==========================================================================

    X = @evalpoly(T_TT,   -0.016_617,
                       +2004.191_898,
                          -0.429_782_9,
                          -0.198_618_34,
                          +0.000_007_578,
                          +0.000_005_928_5)

    Xc = zeros(typeof(X), 5)
    coefs_v = (nut_coefs_iau2006_X0, nut_coefs_iau2006_X1, nut_coefs_iau2006_X2,
               nut_coefs_iau2006_X3, nut_coefs_iau2006_X4)

    @inbounds for i = 1:length(coefs_v)

        coefs = coefs_v[i]
        Xp    = zero(typeof(X))

        for j = 1:size(coefs,1)
            As     = coefs[j, 2]
            Ac     = coefs[j, 3]
            ap     = coefs[j, 4]*M_m  + coefs[j, 5]*M_s  + coefs[j, 6]*u_Mm +
                     coefs[j, 7]*D_s  + coefs[j, 8]*Ω_m  + coefs[j, 9]*λ_M☿ +
                     coefs[j,10]*λ_M♀ + coefs[j,11]*λ_Me + coefs[j,12]*λ_M♂ +
                     coefs[j,13]*λ_M♃ + coefs[j,14]*λ_M♄ + coefs[j,15]*λ_M⛢ +
                     coefs[j,16]*λ_M♆ + coefs[j,17]*p_λ
            sj,cj  = sincos(ap)
            Xp    += (As*sj + Ac*cj)/1e6
        end

        Xc[i] = Xp
    end

    X += @evalpoly(T_TT, Xc[1], Xc[2], Xc[3], Xc[4], Xc[5])

    # Convert to [rad].
    X *= a2r

    # Y position of the CIP
    # ==========================================================================

    Y = @evalpoly(T_TT, -0.006_951,
                        -0.025_896,
                       -22.407_274_7,
                        +0.001_900_59,
                        +0.001_112_526,
                        +0.000_000_135_8)

    Yc = zeros(typeof(Y), 5)
    coefs_v = (nut_coefs_iau2006_Y0, nut_coefs_iau2006_Y1, nut_coefs_iau2006_Y2,
               nut_coefs_iau2006_Y3, nut_coefs_iau2006_Y4)

    @inbounds for i = 1:length(coefs_v)

        coefs = coefs_v[i]
        Yp    = zero(typeof(Y))

        for j = 1:size(coefs,1)
            As     = coefs[j, 2]
            Ac     = coefs[j, 3]
            ap     = coefs[j, 4]*M_m  + coefs[j, 5]*M_s  + coefs[j, 6]*u_Mm +
                     coefs[j, 7]*D_s  + coefs[j, 8]*Ω_m  + coefs[j, 9]*λ_M☿ +
                     coefs[j,10]*λ_M♀ + coefs[j,11]*λ_Me + coefs[j,12]*λ_M♂ +
                     coefs[j,13]*λ_M♃ + coefs[j,14]*λ_M♄ + coefs[j,15]*λ_M⛢ +
                     coefs[j,16]*λ_M♆ + coefs[j,17]*p_λ
            sj,cj  = sincos(ap)
            Yp    += (As*sj + Ac*cj)/1e6
        end

        Yc[i] = Yp
    end

    Y += @evalpoly(T_TT, Yc[1], Yc[2], Yc[3], Yc[4], Yc[5])

    # Convert to [rad].
    Y *= a2r

    # Parameter `s` (CIO locator)
    # ==========================================================================
    #
    # The value `s` provides the position of the CIO on the Equator of the CIP
    # corresponding to the kinematical definition of the non-rotation origin in
    # the GCRS when the CIP is moving with respect to the GCRS between the
    # reference epoch and the epoch due to precession and nutation [1, p. 214].

    # We must convert this term to [arcsec] to match the units.
    #                      ||
    #                   |------|
    s = @evalpoly(T_TT, (-X*Y/2)/a2r + 0.000_094,
                                     +0.003_808_65,
                                     -0.000_122_68,
                                     -0.072_574_11,
                                     +0.000_027_98,
                                     +0.000_015_65)

    sc = zeros(typeof(s), 5)
    coefs_v = (nut_coefs_iau2006_s0, nut_coefs_iau2006_s1, nut_coefs_iau2006_s2,
               nut_coefs_iau2006_s3, nut_coefs_iau2006_s4)

    @inbounds for i = 1:length(coefs_v)

        coefs = coefs_v[i]
        sp    = zero(typeof(s))

        for j = 1:size(coefs,1)
            As     = coefs[j, 2]
            Ac     = coefs[j, 3]
            ap     = coefs[j, 4]*M_m  + coefs[j, 5]*M_s  + coefs[j, 6]*u_Mm +
                     coefs[j, 7]*D_s  + coefs[j, 8]*Ω_m  + coefs[j, 9]*λ_M☿ +
                     coefs[j,10]*λ_M♀ + coefs[j,11]*λ_Me + coefs[j,12]*λ_M♂ +
                     coefs[j,13]*λ_M♃ + coefs[j,14]*λ_M♄ + coefs[j,15]*λ_M⛢ +
                     coefs[j,16]*λ_M♆ + coefs[j,17]*p_λ
            sj,cj  = sincos(ap)
            sp    += (As*sj + Ac*cj)/1e6
        end

        sc[i] = sp
    end

    s += @evalpoly(T_TT, sc[1], sc[2], sc[3], sc[4], sc[5])

    # Convert to [rad].
    s *= a2r

    return X, Y, s
end
