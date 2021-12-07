# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#    Compute the Greenwich Mean Sideral Time (GMST).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] http://www.navipedia.net/index.php/CEP_to_ITRF, accessed 2015-12-01.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export j2000_to_gmst, jd_to_gmst

"""
    j2000_to_gmst(j2000_ut1::Number)

Compute the Greenwich Mean Sideral Time (GMST) \\[rad] given the instant
`j2000_ut1` in J2000.0 reference [UT1].

!!! info
    The algorithm is based in **[1]**.

# References

- **[1]** http://www.navipedia.net/index.php/CEP_to_ITRF, accessed 2015-12-01.
"""
function j2000_to_gmst(j2000_ut1::Number)
	# Julian centuries elapsed from the epoch J2000.0.
	t_ut1 = j2000_ut1 / 36525

	# Greenwich Mean Sideral Time at t_ut1 [s].
    θ_gmst = @evalpoly(
        t_ut1,
        + 67310.54841,
        + (876600.0 * 3600 + 8640184.812866),
        + 0.093104,
        - 6.2e-6
    )

    # Reduce to the interval [0, 86400]s.
    θ_gmst = mod(θ_gmst, 86400)

    # Convert to radian and return.
    return θ_gmst * π / 43200
end

"""
    jd_to_gmst(jd_ut1::Number)

Compute the Greenwich Mean Sideral Time (GMST) \\[rad] for the Julian Day
`jd_ut1` [UT1].

!!! info
    The algorithm is based in **[1]**(p. 188).

# References

- **[1]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
    Microcosm Press, Hawthorn, CA, USA.
"""
jd_to_gmst(jd_ut1::Number) = j2000_to_gmst(jd_ut1 - JD_J2000)
