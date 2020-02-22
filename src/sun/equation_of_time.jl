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

export equation_of_time

"""
    equation_of_time(JD::Number)

Compute the difference between the Sun apparent local time and the Sun mean
local time [rad], which is called Equation of Time, at the Julian Day `JD`. The
algorithm was adapted from [1, p. 178, 277-279].

"""
function equation_of_time(JD::Number)
    # Constants
    deg2rad = π/180

    # Number of Julian centuries from J2000 epoch.
    T_UT1 = (JD-JD_J2000)/36525.0

    # Mean longitude of the Sun.
    λ_m = mod(280.460 + 36000.771T_UT1, 360)

    # Mean anomaly of the Sun.
    #
    # Here, we should use T_TBD (Barycentric Dynamical Time). However, it is
    # sufficient to use T_UT1 because this is a low precision computation [1].
    Ms = mod(357.5291092 + 35999.05034T_UT1, 360)*deg2rad

    # Auxiliary variables.
    sin_Ms  = sin(Ms)
    sin_2Ms = sin(2Ms)

    # Ecliptic latitude of the Sun.
    λ_ecliptic = mod(λ_m + 1.914666471sin_Ms + 0.019994643sin_2Ms, 360)*deg2rad

    # Compute the equation of time [rad].
    eot = -1.914666471sin_Ms     - 0.019994643sin_2Ms +
           2.466sin(2λ_ecliptic) - 0.0053sin(4λ_ecliptic)

    eot*deg2rad
end
