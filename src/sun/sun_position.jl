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
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
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
### function sun_position_i(day::Int, time::Int)

Compute the Sun position represented in the Inertial coordinate frame (J2000).

##### Args

* day: Number of days since 01/01/2000 (J2000.0 epoch).
* time: Time of day [s].

##### Returns

* The vector between the origin of the Inertial coordinate frame (J2000) and the
Sun represented in the Inertial coordinate frame (J2000).

##### Remarks

This function was marked **DEPRECATED** at 26 January, 2016. Use the function
sun_position_i(JD::Real) instead.

"""

function sun_position_i(day::Int, time::Int)
    println("WARNING: The function sun_position_i(day::Int, time::Int) is deprecated!")
    println("Use the function sun_position_i(JD::Real) instead.\n")

    hour = time/(86400.0)
    JD = JD_J2000 + day + hour - 0.5

    sun_position_i(JD)
end

"""
### function sun_position_i(day::Int, time::Int)

Compute the Sun position represented in the Inertial coordinate frame (J2000).

##### Args

* JD: Julian day.

##### Returns

* The vector between the origin of the Inertial coordinate frame (J2000) and the
Sun represented in the Inertial coordinate frame (J2000).

"""

function sun_position_i(JD::Real)
    # Constants
    const deg2rad = pi/180

    # Compute the days from the noon of 1 January, 2000, which is the input to
    # the algorithm.
    D = JD - JD_J2000

    # # # # # # # # # # # # # # # #
    # Ecliptic coordinates
    # # # # # # # # # # # # # # # #

    # Mean anomaly of the Sun [rad].
    g = mod(357.529 + 0.98560028*D, 360)
    g = g*deg2rad

    # Mean longitude of the Sun [deg].
    q = mod(280.459 + 0.98564736*D, 360)

    # Geocentric apparent ecliptic longitude of the Sun (adjusted for
    # aberration) [rad].
    L = mod(q + 1.915*sin(g) + 0.020*sin(2g), 360)
    L = L*deg2rad

    # Distance of the Sun from Earth [m].
    R = (1.00014 - 0.01671*cos(g) - 0.00014*cos(2g))*au2m

    # Obliquity of ecliptic [rad].
    e = (23.439 - 0.00000036*D)*deg2rad

    # # # # # # # # # # # # # # # #
    # Inertial (Celestial) coordinates
    # # # # # # # # # # # # # # # #

    # Right ascension and declination.
    cos_L = cos(L)
    sin_L = sin(L)

    cos_e = cos(e)
    sin_e = sin(e)

    ra_i = atan2( cos_e*sin_L, cos_L )
    dec_i = asin( sin_e*sin_L )

    # Compute the vector from the origin of the Celestial Coordinate Frame to
    # the center of the Sun represented in the Celestial Coordinate Frame.
    cos_ra_i = cos(ra_i)
    sin_ra_i = sin(ra_i)

    cos_dec_i = cos(dec_i)
    sin_dec_i = sin(dec_i)

    S_i = [cos_ra_i*cos_dec_i;
           sin_ra_i*cos_dec_i;
           sin_dec_i]*R
end
