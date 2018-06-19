#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to IAU-76/FK5 precession algorithm.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/transformations/fk5/precession.jl
# =============================================

# Function precession_fk5
# -----------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# In this example, using JD_TT = 2453101.828154745, one gets the following data
# related to the precession:
#
#   ζ = 0.0273055˚
#   Θ = 0.0237306˚
#   z = 0.0273059˚
#
# SatelliteToolbox provides the following results:
#
#   julia> (ζ,Θ,z) = precession_fk5(2453101.828154745)
#   julia> ζ*180/pi
#   0.027305539219877804
#
#   julia> Θ*180/pi
#   0.023730619896279084
#
#   julia> z*180/pi
#   0.02730593931829398
#
################################################################################

@testset "Precession" begin
    (ζ,Θ,z) = precession_fk5(2453101.828154745)

    @test ζ*180/pi ≈ 0.0273055 atol=1e-7
    @test Θ*180/pi ≈ 0.0237306 atol=1e-7
    @test z*180/pi ≈ 0.0273059 atol=1e-7
end
