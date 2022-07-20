# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions to compute RAAN to match a local time in ascending or descending
#   node.
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

export ltan_to_raan, ltdn_to_raan

"""
    ltan_to_raan(JD::Number, asc_node_lt::Number)
    ltan_to_raan(JD::Number, asc_node_lt::Time)

Compute the RAAN [0, 2π] [rad] so that the orbit plane local time of the
ascending node is `ltan` at Julian Day `jd` [UT1].

`ltan` can be represented as a number, indicating the hour, or by a `Time`
object.
"""
ltan_to_raan(jd::Number, ltan::Time) = ltan_to_raan(jd, Dates.value(ltan) / 3.6e12)

function ltan_to_raan(jd::Number, ltan::Number)
    # Get the sun position at noon (UT) represented in the Inertial ref. frame.
    s_i = sun_position_i(jd)

    # Get the desired angle between the Sun and the ascending node [deg].
    α = (ltan - 12) * π / 12

    # Get the right ascension of the Sun in the Inertial ref. frame. This is the
    # Sun apparent local time.
    salt = atan(s_i[2], s_i[1])

    # Get the equation of time to compute the Sun mean local time [rad].
    eot = equation_of_time(jd)

    # Compute the Sun mean local time.
    smlt = salt + eot

    # Compute the desired RAAN in the interval [0, 2π].
    raan = mod(smlt + α, 2π)

    return raan
end

"""
    ltdn_to_raan(JD::Number, asc_node_lt::Number)
    ltdn_to_raan(JD::Number, asc_node_lt::Time)

Compute the RAAN [0, 2π] [rad] so that the orbit plane local time of the
descending node is `ltdn` at Julian Day `jd` [UT1].

`ltdn` can be represented as a number, indicating the hour, or by a `Time`
object.
"""
ltdn_to_raan(jd::Number, ltdn::Time) = ltdn_to_raan(jd, Dates.value(ltdn) / 3.6e12)
ltdn_to_raan(jd::Number, ltdn::Number) = ltan_to_raan(jd::Number, mod(ltdn + 12, 24))

