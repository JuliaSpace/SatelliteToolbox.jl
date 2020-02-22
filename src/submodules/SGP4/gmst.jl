# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Compute the Greenwich Mean Sideral Time (GMST).
#
# Remarks
#
#   This files was copied from SatelliteToolbox.jl because the functionality is
#   required by SGP4 model. The functions here **must not** be exported to avoid
#   interferences.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""
    J2000toGMST(J2000_UT1::Number)

Compute the Greenwich Mean Sideral Time (GMST) \\[rad] given the instant
`J2000_UT1` in J2000.0 reference [UT1].

# Remarks

Based on algorithm in [2] (http://www.navipedia.net/index.php/CEP_to_ITRF),
accessed at 2015-12-01.

"""
function J2000toGMST(J2000_UT1::Number)
	# Julian centuries elapsed from the epoch J2000.0.
	T_UT1 = J2000_UT1/36525

	# Greenwich Mean Sideral Time at T_UT1 [s].
    θ_GMST = @evalpoly(T_UT1, + 67310.54841,
                              + (876600.0*3600 + 8640184.812866),
                              + 0.093104,
                              - 6.2e-6)

    # Reduce to the interval [0, 86400]s.
    θ_GMST = mod(θ_GMST, 86400)

    # Convert to radian and return.
    θ_GMST*π/43200
end

"""
    JDtoGMST(JD_UT1::Number)

Compute the Greenwich Mean Sideral Time (GMST) \\[rad] for the Julian Day
`JD_UT1` [UT1].

# Remarks

Based on algorithm in [1, pp. 188].

"""
function JDtoGMST(JD_UT1::Number)
	J2000toGMST(JD_UT1 - JD_J2000)
end
