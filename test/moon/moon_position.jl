#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to the position of the Moon.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/sun/moon_position.jl
# ================================

# Function moon_position_i
# ------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Example 5-3: Finding the Moon position vector [1, p. 289]
#
#   In April 28, 1994, 00:00 UTC, one gets:
#
#             | -134_240.626 |
#   r_J2000 = | -311_571.590 | km
#             | -126_693.785 |
#
#   Note: it assumes that UTC ≈ UT1 ≈ TDB.
#
################################################################################

@testset "Function moon_position_i" begin
    JD_UTC = DatetoJD(1994,4,28,0,0,0)
    JD_TBD = JD_UTC

    r_J2000 = moon_position_i(JD_TBD)

    @test r_J2000[1]/1000 ≈ -134_240.626 atol = 1e-3
    @test r_J2000[2]/1000 ≈ -311_571.590 atol = 1e-3
    @test r_J2000[3]/1000 ≈ -126_693.785 atol = 1e-3
end
