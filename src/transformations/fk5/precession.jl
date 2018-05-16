#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
#   Functions to compute the precession according to IAU-76/FK5.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-04-18: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export precession_fk5

################################################################################
#                                  Functions
################################################################################

"""
    function precession_fk5(JD_TT::Number)

Compute the angles related to the precession movement in the Julian Day
(Terrestrial Time) `JD_TT` using the theory IAU-76/FK5.

##### Args

* `JD_TT`: Julian day [TT].

##### Returns

The angles (ζ, Θ, z) as described in [1, p. 226-228].

"""
function precession_fk5(JD_TT::Number)
    # Compute the Julian Centuries from `JD_TT`.
    T_TT = (JD_TT - JD_J2000)/36525

    # Compute the angles [arcsec].
    ζ = 2306.2181*T_TT + 0.30188*T_TT^2 + 0.017998*T_TT^3
    Θ = 2004.3109*T_TT - 0.42665*T_TT^2 - 0.041833*T_TT^3
    z = 2306.2181*T_TT + 1.09468*T_TT^2 + 0.018203*T_TT^3

    # Normalize the angles in the interval [0, 86400]s and convert to rad.
    ζ = mod(ζ*pi/648000, 2*pi)
    Θ = mod(Θ*pi/648000, 2*pi)
    z = mod(z*pi/648000, 2*pi)

    # Return the date.
    (ζ, Θ, z)
end
