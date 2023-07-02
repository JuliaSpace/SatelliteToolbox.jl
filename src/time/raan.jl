# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Functions to compute RAAN to match a local time in ascending or descending
#   node.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==========================================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export ltan_to_raan, ltdn_to_raan
export raan_to_ltan, raan_to_ltdn

"""
    ltan_to_raan(ltan::Union{Number, Time}, t::Union{Number, DateTime}) -> Float64

Compute the RAAN [0, 2π] [rad] so that the local time of ascending node (LTAN) is `ltan` at
instant `t` [UT1].

`ltan` can be represented as a `Number`, indicating the hour, or by a `Time` object.

`t` can be represented as a Julian Day [UT1] or `DateTime` [UT1].
"""
function ltan_to_raan(ltan::Time, t::Union{Number, DateTime})
    return ltan_to_raan(Dates.value(ltan) / 3.6e12, t)
end

function ltan_to_raan(ltan::Number, t::Union{Number, DateTime})
    # Get the sun position at noon (UT) represented in the MOD reference frame.
    s_i = sun_position_mod(t)

    # Get the desired angle between the Sun and the ascending node [deg].
    α = (ltan - 12) * π / 12

    # Get the right ascension of the Sun in the MOD reference frame. This is the Sun
    # apparent local time.
    salt = atan(s_i[2], s_i[1])

    # Get the equation of time to compute the Sun mean local time [rad].
    eot = equation_of_time(t)

    # Compute the Sun mean local time.
    smlt = salt + eot

    # Compute the desired RAAN in the interval [0, 2π].
    raan = mod(smlt + α, 2π)

    return raan
end

"""
    ltdn_to_raan(ltdn::Union{Number, Time}, t::Union{Number, DateTime}) -> Float64

Compute the RAAN [0, 2π] [rad] so that the local time of descending node (LTAN) is `ltdn` at
instant `t` [UT1].

`ltdn` can be represented as a `Number`, indicating the hour, or by a `Time` object.

`t` can be represented as a Julian Day [UT1] or `DateTime` [UT1].
"""
function ltdn_to_raan(ltdn::Time, t::Union{Number, DateTime})
    return ltdn_to_raan(Dates.value(ltdn) / 3.6e12, t)
end

function ltdn_to_raan(ltdn::Number, t::Union{Number, DateTime})
    return ltan_to_raan(mod(ltdn + 12, 24), t)
end

"""
    raan_to_ltan(raan::Number, t::Union{Number, DateTime}) -> Float64

Compute the local time of the ascending node (LTAN) given the `raan` at instant `t` [UT1].

`t` can be represented as a Julian Day [UT1] or `DateTime` [UT1].
"""
function raan_to_ltan(raan::Number, t::Union{Number, DateTime})
    # Get the sun position at noon (UT) represented in the MOD refernece frame.
    s_i = sun_position_mod(t)

    # Get the right ascension of the Sun in the MOD reference frame. This is the Sun
    # apparent local time.
    salt = atan(s_i[2], s_i[1])

    # Get the equation of time to compute the Sun mean local time [rad].
    eot = equation_of_time(t)

    # Compute the Sun mean local time.
    smlt = salt + eot

    # Get the angle between the Sun and the ascending node [rad].
    α = mod(raan - smlt, 2π)

    # Get the LTAN [hour].
    ltan = mod(α * 12 / π + 12, 24)

    return ltan
end

"""
    raan_to_ltdn(raan::Number, t::Union{Number, DateTime}) -> Float64

Compute the local time of the ascending node (LTDN) given the `raan` at instant `t` [UT1].

`t` can be represented as a Julian Day [UT1] or `DateTime` [UT1].
"""
function raan_to_ltdn(raan::Number, t::Union{Number, DateTime})
    return mod(raan_to_ltan(raan, t) + 12, 24)
end
