# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Compute the fundamental arguments related to the IAU-2006 theory.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export luni_solar_args_iau2006, planetary_args_iau2006

"""
    luni_solar_args_iau2006(JD_TT::Number)

Compute the fundamental arguments related to the luni-solar effect for the
IAU-2006 theory [^1](p. 211).

The returned values are in [rad].

# References

[^1]: Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
function luni_solar_args_iau2006(JD_TT::Number)
    # Compute the Julian Centuries from `JD_TT`.
    T_TT = (JD_TT - JD_J2000) / 36525

    # Auxiliary variables
    # ===================

    a2d = 1 / 3600
    d2r = π / 180

    # Delaunay arguments of the Sun and Moon
    # ======================================
    #
    # Evaluate the Delaunay arguments associated with the Moon and the Sun in
    # [arcsec] [1, p. 210].

    M_m  = @evalpoly(
        T_TT,
        +485868.249036,
        +1717915923.2178,
        +31.8792,
        +0.051635,
        -0.00024470
    )

    M_s  = @evalpoly(
        T_TT,
        +1287104.79305,
        +129596581.0481,
        -0.5532,
        +0.000136,
        -0.00001149
    )

    u_Mm = @evalpoly(
        T_TT,
        +335779.526232,
        +1739527262.8478,
        -12.7512,
        -0.001037,
        +0.00000417
    )

    D_s  = @evalpoly(
        T_TT,
        +1072260.70369,
        +1602961601.2090,
        -6.3706,
        +0.006593,
        -0.00003169
    )

    Ω_m  = @evalpoly(
        T_TT,
        +450160.398036,
        -6962890.5431,
        +7.4722,
        +0.007702,
        -0.00005939
    )

    # Convert to the interval [0,2π].
    M_s  = mod( M_s * a2d, 360) * d2r
    M_m  = mod( M_m * a2d, 360) * d2r
    u_Mm = mod(u_Mm * a2d, 360) * d2r
    D_s  = mod( D_s * a2d, 360) * d2r
    Ω_m  = mod( Ω_m * a2d, 360) * d2r

    return M_s, M_m, u_Mm, D_s, Ω_m
end

"""
    planetary_args_iau2006(JD_TT::Number)

Compute the fundamental arguments related to the planetary effects for the
IAU-2006 theory [^1](p. 211).

The returned values are in [rad].

# References

[^1]: Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
function planetary_args_iau2006(JD_TT::Number)
    # Compute the Julian Centuries from `JD_TT`.
    T_TT = (JD_TT - JD_J2000) / 36525

    # Mean Heliocentric longitudes of the planets.
    #
    # TODO: In the example in [1, p. 221], the value related to Uranus is slight
    # different. The value used here is shown in [1](p. 211).

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

    return λ_M☿, λ_M♀, λ_Me, λ_M♂, λ_M♃, λ_M♄, λ_M⛢, λ_M♆, p_λ
end
