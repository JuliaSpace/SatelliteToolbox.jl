# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Compute the moon position.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorne, CA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export moon_position_i

"""
    moon_position_i(JD_TDB::Number)

Compute the Moon position represented in the IAU-76/FK5 (mean-equator,
mean-equinox), also called as J2000, at the Julian Day `JD`. The algorithm was
adapted from [1, p. 288].

"""
function moon_position_i(JD_TDB::Number)
    # Constants.
    deg2rad = π/180

    # Number of Julian centuries from J2000 epoch.
    T_TDB = (JD_TDB-JD_J2000)/36525.0

    # Auxiliary computation to improve performance.
    sin1, cos1 = sincos( (134.9 + 477_198.85T_TDB)*deg2rad )
    sin2, cos2 = sincos( (259.2 - 413_335.38T_TDB)*deg2rad )
    sin3, cos3 = sincos( (235.7 + 890_534.23T_TDB)*deg2rad )
    sin4, cos4 = sincos( (269.9 + 954_397.70T_TDB)*deg2rad )
    sin5       = sind(357.5 +  35_999.05T_TDB)
    sin6       = sind(186.6 + 966_404.05T_TDB)

    # Ecliptic latitude of the Moon [deg].
    λₑ = 218.32 + 481_267.8813T_TDB +
         6.29sin1 - 1.27sin2 + 0.66sin3 + 0.21sin4 - 0.19sin5 - 0.11sin6

    # Ecliptic longitude of the Moon [deg].
    ϕₑ = 5.13sind( 93.3 + 483_202.03T_TDB) + 0.28sind(228.2 + 960_400.87T_TDB) -
         0.28sind(318.3 +   6_003.18T_TDB) - 0.17sind(217.6 - 407_332.20T_TDB)

    # Parallax [deg].
    P = 0.9508 + 0.0518cos1 + 0.0095cos2 + 0.0078cos3 + 0.0028cos4

    # Obliquity of the ecliptic [deg].
    ϵ = @evalpoly(T_TDB, 23.439_291, -0.013_004_2, -1.64e-7, +5.04e-7)

    # Convert to radians and limit to the interval [0,2π].
    λₑ = mod2pi(λₑ*deg2rad)
    ϕₑ = mod2pi(ϕₑ*deg2rad)
    P  = mod2pi( P*deg2rad)
    ϵ  = mod2pi( ϵ*deg2rad)

    # Compute the distance from Earth to the Moon [m].
    r = R0/sin(P)

    # Auxiliary variables.
    sinλ, cosλ = sincos(λₑ)
    sinϕ, cosϕ = sincos(ϕₑ)
    sinϵ, cosϵ = sincos(ϵ)

    # Compute the Moon vector represented in J2000 (IAU-76/KF5 mean-equator,
    # mean-equinox).
    r_J2000 = SVector{3,Float64}(r*(cosϕ*cosλ),
                                 r*(cosϵ*cosϕ*sinλ - sinϵ*sinϕ),
                                 r*(sinϵ*cosϕ*sinλ + cosϵ*sinϕ))

    return r_J2000
end
