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
#   General tests.
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
# 2018-04-04: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/orbit/anomalies.jl
# ==============================

println("$(c)Testing functions in file: ./src/orbit/anomalies.jl$d")
println("$(c)---------------------------------------------------$d")
println("")

# Function M_to_E
# ---------------

println("    Testing function M_to_E...")

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 2-1: Using Kepler's Equation [1, p. 66-67].
#
#   M = 235.4°, e = 0.4 ===> E = 220.512074767522º
#
# Using this function, the following result was obtained:
#
#   julia> M_to_E(0.4,235.4*pi/180)*180/pi
#   220.51207476752163
#
################################################################################

E = M_to_E(0.4, 235.4*pi/180)*180/pi
@test E ≈ 220.512074767522 atol=1e-12

println("        $(b)Test passed!$d")
println("")
