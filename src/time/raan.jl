## Description #############################################################################
#
# Functions to convert between the RAAN and the local time of the ascending or descending
# node. The RAAN is referenced to the MOD (mean equator and equinox of date) frame, since
# the Sun position is obtained in this frame.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorne, CA.
#
############################################################################################

export ltan_to_raan, ltdn_to_raan
export raan_to_ltan, raan_to_ltdn

"""
    ltan_to_raan(ltan::Union{Number, Time}, t::Union{Number, DateTime}) -> Float64

Compute the RAAN [rad] in the interval [0, 2π], referenced to the MOD (mean equator and
equinox of date) frame, so that the local time of ascending node (LTAN) is `ltan` at
instant `t` [UT1].

`ltan` can be represented as a `Number`, indicating the hour, or by a `Time` object.

`t` can be represented as a Julian Day [UT1] or `DateTime` [UT1]. Strictly speaking, the
Sun position model expects `t` in the TDB time scale. However, using UT1 leads to a
negligible error in this low-precision computation.
"""
function ltan_to_raan(ltan::Time, t::Union{Number, DateTime})
    return ltan_to_raan(Dates.value(ltan) / 3.6e12, t)
end

function ltan_to_raan(ltan::Number, t::Union{Number, DateTime})
    # Get the desired angle between the mean Sun and the ascending node [rad].
    α = (ltan - 12) * π / 12

    # Compute the desired RAAN in the interval [0, 2π].
    raan = mod(_mean_sun_right_ascension(t) + α, 2π)

    return raan
end

"""
    ltdn_to_raan(ltdn::Union{Number, Time}, t::Union{Number, DateTime}) -> Float64

Compute the RAAN [rad] in the interval [0, 2π], referenced to the MOD (mean equator and
equinox of date) frame, so that the local time of descending node (LTDN) is `ltdn` at
instant `t` [UT1].

`ltdn` can be represented as a `Number`, indicating the hour, or by a `Time` object.

`t` can be represented as a Julian Day [UT1] or `DateTime` [UT1]. Strictly speaking, the
Sun position model expects `t` in the TDB time scale. However, using UT1 leads to a
negligible error in this low-precision computation.
"""
function ltdn_to_raan(ltdn::Time, t::Union{Number, DateTime})
    return ltdn_to_raan(Dates.value(ltdn) / 3.6e12, t)
end

function ltdn_to_raan(ltdn::Number, t::Union{Number, DateTime})
    return ltan_to_raan(mod(ltdn + 12, 24), t)
end

"""
    raan_to_ltan(raan::Number, t::Union{Number, DateTime}) -> Float64

Compute the local time of the ascending node (LTAN) [hour] given the `raan` [rad] at
instant `t` [UT1]. The `raan` must be referenced to the MOD (mean equator and equinox of
date) frame.

`t` can be represented as a Julian Day [UT1] or `DateTime` [UT1]. Strictly speaking, the
Sun position model expects `t` in the TDB time scale. However, using UT1 leads to a
negligible error in this low-precision computation.
"""
function raan_to_ltan(raan::Number, t::Union{Number, DateTime})
    # Get the angle between the mean Sun and the ascending node [rad].
    α = mod(raan - _mean_sun_right_ascension(t), 2π)

    # Get the LTAN [hour].
    ltan = mod(α * 12 / π + 12, 24)

    return ltan
end

"""
    raan_to_ltdn(raan::Number, t::Union{Number, DateTime}) -> Float64

Compute the local time of the descending node (LTDN) [hour] given the `raan` [rad] at
instant `t` [UT1]. The `raan` must be referenced to the MOD (mean equator and equinox of
date) frame.

`t` can be represented as a Julian Day [UT1] or `DateTime` [UT1]. Strictly speaking, the
Sun position model expects `t` in the TDB time scale. However, using UT1 leads to a
negligible error in this low-precision computation.
"""
function raan_to_ltdn(raan::Number, t::Union{Number, DateTime})
    return mod(raan_to_ltan(raan, t) + 12, 24)
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

"""
    _mean_sun_right_ascension(t::Union{Number, DateTime}) -> Float64

Compute the right ascension [rad] of the mean Sun in the MOD (mean equator and equinox of
date) frame at the instant `t` [UT1], which can be represented as a Julian Day or
`DateTime`.

The mean Sun is a fictitious body that moves uniformly along the celestial equator, and its
right ascension defines the mean solar time. It is obtained by adding the equation of time
to the right ascension of the apparent Sun.
"""
function _mean_sun_right_ascension(t::Union{Number, DateTime})
    # Get the apparent Sun position at the instant `t` represented in the MOD frame.
    s_mod = sun_position_mod(t)

    # Get the right ascension of the apparent Sun in the MOD frame [rad].
    α_sun = atan(s_mod[2], s_mod[1])

    # The equation of time is the apparent solar time minus the mean solar time. Hence, the
    # right ascension of the mean Sun is the right ascension of the apparent Sun plus the
    # equation of time.
    return α_sun + equation_of_time(t)
end
