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
#    Compute the satellite velocity.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2016-09-14: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

import Rotations: angle2dcm, angle2dcm!

export satellite_velocity_J0_i

"""
### function satellite_velocity_J0_i(a::Number, e::Number, i::Number, RAAN::Number, w::Number, f::Number)

Compute the satellite velocity in the Inertial coordinate frame given the
orbital elements `a`, `e`, `i`, `RAAN`, `w`, and `f` neglecting all the
perturbations.

##### Args

* JD: Julian day.
* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].
* RAAN: Right ascension of the ascending node [rad].
* w: Argument of perigee [rad].
* f: True anomaly [rad].

##### Returns

* The satellite velocity vector represented in the Inertial coordinate frame
  [m/s].

"""

################################################################################
#                                 TEST RESULTS
################################################################################
#
# This function was tested using the example in [1, pp. 119-120].
#
# Scenario
# ========
#
#    p    = 11067.790 km
#    e    = 0.83285
#    i    = 87.87°
#    RAAN = 227.89°
#    w    = 53.38°
#    f    = 92.335°
#    a    = p/(1-e^2)
#
#
#                    +----------+----------+-----------+
#                    | X [km/s] | Y [km/s] | Z [km/s]  |
# +------------------+----------+----------+-----------+
# | Algorithm in [1] | 4.902276 | 5.533124 | -1.975709 |
# +------------------+----------+----------+-----------+
# | This function    | 4.902279 | 5.533140 | -1.975710 |
# +------------------+----------+----------+-----------+
# | Difference       | 0.000003 | 0.000016 |  0.000001 |
# +------------------+----------+----------+-----------+
#
# References
# ==========
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th
# ed. Microcosm Press. Hawthorn, CA, USA.
#

function satellite_velocity_J0_i(a::Number, e::Number, i::Number, RAAN::Number,
                                 w::Number, f::Number)
    # Satellite angular velocity without perturbations.
    n = period(a, e, i, :J0)

    # Compute the velocity in the q-frame, defined as:
    #   - The X axis points towards the perigee;
    #   - The Z axis is perpendicular to the orbital plane (right-hand);
    #   - The Y axis completes a right-hand coordinate system.
    v_q = n*a/sqrt(1-e^2)*[-sin(f); e+cos(f); 0]

    # Compute the matrix that rotates from the q-frame to the Inertial
    # coordinate Frame.
    Dsi = Array{Float64}(3,3)
    angle2dcm!(Dsi, RAAN, i, w, "ZXZ")

    # Compute the satellite velocity vector represented in the Inertial
    # coordinate frame.
    v_i = Dsi'*v_q
end
