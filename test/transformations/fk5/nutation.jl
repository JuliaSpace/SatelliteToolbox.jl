#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to IAU-76/FK5 nutation algorithm.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/transformations/fk5/nutation.jl
# ===========================================

# Function nutation_fk5
# ---------------------

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
# related to the nutation:
#
#   mϵ_1980 =  23.4387368˚
#   Δϵ_1980 =   0.0020316˚
#   Δψ_1980 =  -0.0034108˚
#
# SatelliteToolbox provides the following results:
#
#   julia> (mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(2453101.828154745)
#   (0.40908312815877673, 3.5458763448549555e-5, -5.953027070867465e-5)
#
#   julia> mϵ_1980*180/pi
#   23.438736713507268
#
#   julia> Δϵ_1980*180/pi
#   0.0020316374923546377
#
#   julia> Δψ_1980*180/pi
#   -0.0034108332648783257
#
# Scenario 02
# ===========
#
# Example in [1, pp. 233 and 234].
#
# In the day 182.78495062 of year 2000, if one uses only 4 terms in nutation
# computation (IAU-76/FK5 model), then one gets:
#
#   Δψ_1980 = -0.004250260˚
#   Δϵ_1980 = -0.001260854˚
#   mϵ_1980 = +23.43922657˚
#
# If the same algorithm is applied to the day 179.78495062 of year 2000, the one
# gets:
#
#   Δψ_1980 = -0.004337544˚
#   Δϵ_1980 = -0.001247061˚
#   mϵ_1980 = +23.43922657˚
#
################################################################################

@testset "Nutation" begin
    # Scenario 01
    # ===========

    (mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(2453101.828154745)

    @test mϵ_1980*180/pi ≈ 23.4387368 atol=1e-7
    @test Δϵ_1980*180/pi ≈  0.0020316 atol=1e-7
    @test Δψ_1980*180/pi ≈ -0.0034108 atol=1e-7

    # Scenario 02
    # ===========

    JD_UTC = DatetoJD(2000,1,1,0,0,0) + 182.78495062 - 1
    JD_TT  = JD_UTCtoTT(JD_UTC)

    (mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(JD_TT, 4)

    @test mϵ_1980*180/pi ≈ 23.43922657  atol=1e-6
    @test Δϵ_1980*180/pi ≈ -0.001260854 atol=1e-8
    @test Δψ_1980*180/pi ≈ -0.004250260 atol=1e-8

    JD_UTC = DatetoJD(2000,1,1,0,0,0) + 179.78495062 - 1
    JD_TT  = JD_UTCtoTT(JD_UTC)

    (mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(JD_TT, 4)

    @test mϵ_1980*180/pi ≈ 23.43922657  atol=1e-6
    @test Δϵ_1980*180/pi ≈ -0.001247061 atol=1e-8
    @test Δψ_1980*180/pi ≈ -0.004337544 atol=1e-8
end
