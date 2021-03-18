#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions to compute the precession-nutation rotation according to IAU-2006
#   theory.
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

export precession_nutation_iau2006

"""
    precession_nutation_iau2006(JD_TT::Number)

Compute the coordinates `X`, `Y`, and `s` related to the Celestial Intermediate
Pole (CIP) with respect to the Geocentric Celestial Reference Frame (GCRF). This
accounts for the effects of both precession and nutation of the CIP.

# Returns

* The coordinate `X` of the CIP w.r.t. the GCRF;
* The coordinate `Y` of the CIP w.r.t. the GCRF;
* The CIO locator `s` that provides the position of the CIO on the Equator of
  the CIP corresponding to the kinematical definition of the non-rotation origin
  in the GCRS when the CIP is moving with respect to the GCRS between the
  reference epoch and the epoch due to precession and nutation [1, p. 214].

"""
function precession_nutation_iau2006(JD_TT::Number)
    # Compute the Julian Centuries from `JD_TT`.
    T_TT = (JD_TT - JD_J2000)/36525

    # Auxiliary variables
    # ===================

    a2d = 1/3600
    d2r = π/180
    a2r = a2d*d2r

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

    # Planetary effects of the nutation and obliquity of the ecliptic
    # ===============================================================

    # Mean Heliocentric longitudes of the planets.
    #
    # TODO: In the example in [1, p. 221], the value related to Uranus is slight
    # different. The value used here is shown in [1, p. 211].

    λ_M☿ = @evalpoly(T_TT, 4.402_608_842, 2608.790_314_157_4)
    λ_M♀ = @evalpoly(T_TT, 3.176_146_697, 1021.328_554_621_1)
    λ_Me = @evalpoly(T_TT, 1.753_470_314,  628.307_584_999_1)
    λ_M♂ = @evalpoly(T_TT, 6.203_480_913,  334.061_242_670_0)
    λ_M♃ = @evalpoly(T_TT, 0.599_546_497,   52.969_096_264_1)
    λ_M♄ = @evalpoly(T_TT, 0.874_016_757,   21.329_910_496_0)
    λ_M⛢ = @evalpoly(T_TT, 5.481_293_872,    7.478_159_856_7)
    λ_M♆ = @evalpoly(T_TT, 5.311_886_287,    3.813_303_563_8)

    # General precession in longitude.
    p_λ = @evalpoly(T_TT, 0, 0.024_381_75, 0.000_005_386_91)

    # Convert to the interval [0,2π].
    λ_M☿ = mod(λ_M☿, 2π)
    λ_M♀ = mod(λ_M♀, 2π)
    λ_Me = mod(λ_Me, 2π)
    λ_M♂ = mod(λ_M♂, 2π)
    λ_M♃ = mod(λ_M♃, 2π)
    λ_M♄ = mod(λ_M♄, 2π)
    λ_M⛢ = mod(λ_M⛢, 2π)
    λ_M♆ = mod(λ_M♆, 2π)
    p_λ  = mod(p_λ,  2π)

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
