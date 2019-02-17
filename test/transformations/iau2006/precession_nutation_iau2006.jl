#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to IAU-2006 precession-nutation algorithm.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/transformations/precession_nutation_iau2006.jl
# ==========================================================

# Function precession_nutation_iau2006
# ------------------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-14: Performing an IAU-2000 reduction [1, p. 220]
#
# According to this example, using JD_TT = 2453101.828154745, one gets:
#
#   X = 80.531880"
#   Y =  7.273921"
#   s = -0.003027"
#
################################################################################

@testset "Function precession_nutation_iau2006" begin
    JD_TT = 2453101.828_154_745

    X,Y,s = precession_nutation_iau2006(JD_TT)

    @test X*180/pi*3600 ≈ 80.531880 atol = 5e-5
    @test Y*180/pi*3600 ≈  7.273921 atol = 5e-5
    @test s*180/pi*3600 ≈ -0.003027 atol = 1e-6
end
