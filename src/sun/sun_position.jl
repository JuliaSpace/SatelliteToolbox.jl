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
    # Number of Julian centuries from J2000 epoch.
    T_UT1 = (JD-JD_J2000)/36525.0

    # Mean longitude of the Sun [deg].
    λ_m = mod(280.460 + 36000.771*T_UT1, 360.0)

    # Mean anomaly of the Sun [deg].
    #
    # Here, we should use T_TBD (Barycentric Dynamical Time). However, it is
    # sufficient to use T_UT1 because this is a low precision computation [1].
    Ms = mod(357.5291092 + 35999.05034*T_UT1, 360.0)

    # Convert Ms to [rad].
    Ms = deg2rad(Ms)

    # Ecliptic latitude of the Sun [deg].
    λ_e = mod(λ_m + 1.914666471*sin(Ms) + 0.019994643*sin(2Ms), 360.0)

    # Obliquity of the ecliptic [deg].
    #
    # Here, we are also approximating T_TBD ≈ T_UT1.
    ϵ = mod(23.439291 - 0.0130042*T_UT1, 360.0)

    # Auxiliary variables.
    sin_Ms  , cos_Ms  = sincos(Ms)
    sin_2Ms , cos_2Ms = sincos(2Ms)
    sin_ϵ   , cos_ϵ   = sincos(deg2rad(ϵ))
    sin_λ_e , cos_λ_e = sincos(deg2rad(λ_e))

    # Distance of the Sun from Earth [m].
    r = ( 1.000140612 - 0.016708617cos_Ms - 0.000139589cos_2Ms )*au2m

    # Compute the Sun vector represented in the Mean Equinox of Date (MOD).
    S_MOD = [      r*cos_λ_e;
             r*cos_ϵ*sin_λ_e;
             r*sin_ϵ*sin_λ_e;]

    # TODO: This vector must be transformed to the inertial reference frame that
    # will be used in this toolbox (J2000, True of Date, etc.).
end
