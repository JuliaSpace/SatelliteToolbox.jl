# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Compute the sun position.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] The Astronomical Almanac for the year 2000 (p. C24).
#
#   [2] http://aa.usno.navy.mil/faq/docs/SunApprox.php
#
#   [3] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorne, CA.
#
#   [4] The Astronomical Almanac for the year 2006.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export sun_position_i

"""
    function sun_position_i(JD::Number)

Compute the Sun position represented in the Mean Equinox of Date (MOD) at the
Julian Day `JD`. The algorithm was adapted from [3, p. 277-279].

"""
function sun_position_i(JD::Number)
    # Constants.
    deg2rad = π/180

    # Number of Julian centuries from J2000 epoch.
    T_UT1 = (JD-JD_J2000)/36525.0

    # Here, we will assume that T_TBD ≈ T_UT1. This can be done because we are
    # using a low-precision algorithm [3].
    T_TBD = T_UT1

    # Mean anomaly of the Sun [deg].
    Ms = 357.529_109_2 + 35_999.050_34T_TBD

    # Convert Ms to [rad] and limit to the interval [0,2π].
    Ms = mod2pi(Ms*deg2rad)

    # Compute auxiliary variables.
    sinMs,  cosMs  = sincos(Ms)
    sin2Ms, cos2Ms = sincos(2Ms)

    # Mean longitude of the Sun [deg].
    λ_m = 280.460 + 36_000.771T_UT1

    # Ecliptic latitude of the Sun [deg].
    λ_e = λ_m + 1.914_666_471sinMs + 0.019_994_643sin2Ms

    # Obliquity of the ecliptic [deg].
    ϵ = 23.439_291 - 0.013_004_2T_TBD

    # Convert λ_e and ϵ to [rad] and limit to the interval [0,2π].
    λ_e = mod2pi(λ_e*deg2rad)
    ϵ   = mod2pi(  ϵ*deg2rad)

    # Auxiliary variables.
    sinϵ   , cosϵ   = sincos(ϵ)
    sinλ_e , cosλ_e = sincos(λ_e)

    # Distance of the Sun from Earth [m].
    r = ( 1.000_140_612 - 0.016_708_617cosMs - 0.000_139_589cos2Ms )*au2m

    # Compute the Sun vector represented in the Mean Equinox of Date (MOD).
    S_MOD = SVector{3,Float64}(r*cosλ_e, r*cosϵ*sinλ_e, r*sinϵ*sinλ_e)

    # TODO: This vector must be transformed to the inertial reference frame that
    # will be used in this toolbox (J2000, True of Date, etc.).
end
