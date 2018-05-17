# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Compute the equation of time.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorne, CA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2017-08-03: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export equation_of_time

"""
    function equation_of_time(JD::Number)

Compute the difference between the Sun apparent local time and the Sun mean
local time, which is called Equation of Time, at the Julian Day `JD`. The
algorithm was adapted from [1, p. 178, 277-279].

##### Args

* `JD`: Julian day.

##### Returns

The equation of time [rad].

"""
function equation_of_time(JD::Number)
    # Constants
    deg2rad = pi/180

    # Number of Julian centuries from J2000 epoch.
    T_UT1 = (JD-JD_J2000)/36525.0

    # Mean longitude of the Sun.
    λ_m = mod(280.460 + 36000.771*T_UT1, 360.0)

    # Mean anomaly of the Sun.
    #
    # Here, we should use T_TBD (Barycentric Dynamical Time). However, it is
    # sufficient to use T_UT1 because this is a low precision computation [1].
    Ms = mod(357.5291092 + 35999.05034*T_UT1, 360.0)

    # Ecliptic latitude of the Sun.
    λ_ecliptic = mod(λ_m + 1.914666471*sin(Ms*deg2rad) +
                           0.019994643*sin(2*Ms*deg2rad), 360.0)

    # Convert the angles to radians.
    Ms         *= deg2rad
    λ_ecliptic *= deg2rad

    # Compute the equation of time [rad].
    eot = -1.914666471*sin(Ms) - 0.019994643*sin(2*Ms) +
           2.466*sin(2*λ_ecliptic) - 0.0053*sin(4*λ_ecliptic)

    eot*deg2rad
end
