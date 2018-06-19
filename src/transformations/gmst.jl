#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    Compute the Greenwich Mean Sideral Time (GMST).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export J2000toGMST, JDtoGMST

"""
    function J2000toGMST(J2000_UT1::Number)

Compute the Greenwich Mean Sideral Time (GMST) given the instant `J2000` in
J2000.0 reference [UT1].

# Args

* `J2000_UT1`: Instant in J2000.0 reference [UT1].

# Returns

* Greenwich mean sideral time [rad].

# Remarks

Based on algorithm in [2] (http://www.navipedia.net/index.php/CEP_to_ITRF),
accessed at 2015-12-01.

"""
function J2000toGMST(J2000_UT1::Number)
	# Julian centuries elapsed from the epoch J2000.0.
	T_UT1 = J2000_UT1/36525

	# Greenwich Mean Sideral Time at T_UT1 [s].
    θ_GMST = 67310.54841 + (876600.0*3600 + 8640184.812866)*T_UT1 +
                            0.093104*T_UT1^2 -
                            6.2e-6*T_UT1^3

    # Reduce to the interval [0, 86400]s.
    θ_GMST = mod(θ_GMST, 86400)

    # Convert to radian and return.
    θ_GMST*pi/43200
end

"""
    function JDtoGMST(JD_UT1::Number)

Compute the Greenwich Mean Sideral Time (GMST) for a Julian Day `JD` [UT1].

# Args

* `JD_UT1`: Julian day [UT1].

# Returns

* Greenwich mean sideral time [rad].

# Remarks

Based on algorithm in [1, pp. 188].

"""
function JDtoGMST(JD_UT1::Number)
	J2000toGMST(JD_UT1 - JD_J2000)
end
