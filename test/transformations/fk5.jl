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
#   Tests related to IAU-76/FK5 transformations.
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

println("$(c)Testing functions in file: ./src/transformations/fk5/nutation.jl$d")
println("$(c)----------------------------------------------------------------$d")
println("")

# Function nutation_fk5
# ---------------------

println("    Testing function nutation_fk5...")

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
# SatToolbox provides the following results:
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

(mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(2453101.828154745)

@test mϵ_1980*180/pi ≈ 23.4387368 atol=1e-7
@test Δϵ_1980*180/pi ≈  0.0020316 atol=1e-7
@test Δψ_1980*180/pi ≈ -0.0034108 atol=1e-7

println("        $(b)Test passed!$d")
println("")

# File: ./src/transformations/fk5/precession.jl
# =============================================

println("$(c)Testing functions in file: ./src/transformations/fk5/precession.jl$d")
println("$(c)------------------------------------------------------------------$d")
println("")

# Function precession_fk5
# -----------------------

println("    Testing function precession_fk5...")

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
# SatToolbox provides the following results:
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

(ζ,Θ,z) = precession_fk5(2453101.828154745)

@test ζ*180/pi ≈ 0.0273055 atol=1e-7
@test Θ*180/pi ≈ 0.0237306 atol=1e-7
@test z*180/pi ≈ 0.0273059 atol=1e-7

println("        $(b)Test passed!$d")
println("")
