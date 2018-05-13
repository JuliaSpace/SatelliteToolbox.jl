#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to IAU-76/FK5 nutation algorithm.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S (2006). Revisiting
#       Spacetrack Report #3: Rev1. AIAA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-04-18: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
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
################################################################################

@testset "Nutation" begin
    (mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(2453101.828154745)

    @test mϵ_1980*180/pi ≈ 23.4387368 atol=1e-7
    @test Δϵ_1980*180/pi ≈  0.0020316 atol=1e-7
    @test Δψ_1980*180/pi ≈ -0.0034108 atol=1e-7
end
