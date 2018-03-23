# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    Compute the sun position.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References:
#    [1] The Astronomical Almanac for the year 2000 (p. C24).
#    [2] http://aa.usno.navy.mil/faq/docs/SunApprox.php
#    [3] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#        Microcosm Press, Hawthorne, CA.
#    [4] The Astronomical Almanac for the year 2006.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2017-08-04: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#    All DEPRECATED functions were removed.
#
# 2017-08-03: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#    - Update the function `sun_position_i` to use the algorithm in [3], which
#      was obtained from [4]. This is more accurate than the previous version.
#      However, it was noted that the output is represented in the Mean Equinox
#      of Date (MOD). Then it must be converted to the inertial system that will
#      be used (J2000, True of Data, etc.).
#
# 2016-01-26: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#    - The function sun_position_i was changed to receive the Julian Day. The
#      old function is now marked DEPRECATED.
#
# 2014-07-28: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    - The sun vector is returned in meters.
#    - The Celestial Coordinate Frame is now called Inertial coordinate frame
#      (i).
#
# 2014-07-25: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    - Change the algorithm to the one presented in [2].
#    - The function sun_position_i now returns the vector from the origin of the
#      Celestial Coordinate Frame to the center of the Sun represented in the
#      Celestial Coordinate Frame.
#
# 2014-07-16: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export sun_position_i

"""
### function sun_position_i(JD::Number)

Compute the Sun position represented in the Mean Equinox of Date (MOD) at the
Julian Day `JD`. The algorithm was adapted from [3, p. 277-279].

##### Args

* JD: Julian day.

##### Returns

The vector between the origin of the MOD and the Sun represented in the MOD.

"""

function sun_position_i(JD::Number)
    # Constants
    const deg2rad = pi/180

    # Number of Julian centuries from J2000 epoch.
    T_UT1 = (JD-JD_J2000)/36525.0

    # Mean longitude of the Sun [deg].
    λ_m = mod(280.460 + 36000.771*T_UT1, 360.0)

    # Mean anomaly of the Sun [deg].
    #
    # Here, we should use T_TBD (Barycentric Dynamical Time). However, it is
    # sufficient to use T_UT1 because this is a low precision computation [1].
    Ms = mod(357.5291092 + 35999.05034*T_UT1, 360.0)

    # Ecliptic latitude of the Sun [deg].
    λ_ecliptic = mod(λ_m + 1.914666471*sin(Ms*deg2rad) +
                           0.019994643*sin(2*Ms*deg2rad), 360.0)

    # Obliquity of the ecliptic [deg].
    #
    # Here, we are also approximating T_TBD ≈ T_UT1.
    ϵ = mod(23.439291 - 0.0130042*T_UT1, 360.0)

    # Distance of the Sun from Earth [m].
    r = (1.000140612 - 0.016708617*cos(Ms) - 0.000139589*cos(2*Ms))*au2m

    # Compute the Sun vector represented in the Mean Equinox of Date (MOD).
    S_MOD = [r*cos(λ_ecliptic*deg2rad);
             r*cos(ϵ*deg2rad)*sin(λ_ecliptic*deg2rad);
             r*sin(ϵ*deg2rad)*sin(λ_ecliptic*deg2rad);]

    # TODO: This vector must be transformed to the inertial reference frame that
    # will be used in this toolbox (J2000, True of Date, etc.).
end
