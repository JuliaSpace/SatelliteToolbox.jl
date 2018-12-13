#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to the position of the Sun.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/sun/sun_position.jl
# ===============================

# Function sun_position_i
# -----------------------

################################################################################
#                                 Test Results
################################################################################
#
# Example 5-1: Finding the Sun position vector [1, p. 280]
#
#   Using JD_UTC ≈ JD_UT1 = 2453827.5, one gets
#
#           | 0.9771945 |
#   r_MOD = | 0.1924424 | AU
#           | 0.0834308 |
#
################################################################################

@testset "Function sun_position_i" begin
    JD_UTC = 2453827.5
    JD_UT1 = JD_UTC

    s_MOD = sun_position_i(JD_UT1)/au2m

    @test s_MOD[1] ≈ 0.9771945 atol = 2e-6
    @test s_MOD[2] ≈ 0.1924424 atol = 2e-6
    @test s_MOD[3] ≈ 0.0834308 atol = 2e-6
end
