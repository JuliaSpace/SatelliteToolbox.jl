#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to orbit functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/orbit/anomalies.jl
# ==============================

# Function M_to_E
# ---------------

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

@testset "Function M_to_E" begin
    # Float64
    # ==========================================================================

    E = M_to_E(0.4, 235.4 * π / 180) * 180 / π
    @test E ≈ 220.512074767522 atol=1e-12

    # Float32
    # ==========================================================================

    E = M_to_E(0.4f0, 235.4f0 * π / 180) * 180 / π
    @test E ≈ 220.512074767522 atol=1e-5
end

# File: ./src/transformation/gmst.jl
# ==================================

# Function jd_to_gmst
# -----------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-5: Finding GMST and LST (Method 1) [1, p. 188].
#
# Considering the Julian Day [UT1] 2448855.009722, the Greenwich Mean Sideral
# Time was computed as 152.578787810°.
#
# Using SatelliteToolbox, the following was obtained:
#
#   julia> jd_to_gmst(2448855.009722)*180/pi
#   152.57870762832462
#
# NOTE: The difference was also found by replicating the algorithm in MATLAB.
#
################################################################################

@testset "Function jd_to_gmst" begin
    θ_GMST = jd_to_gmst(2448855.009722)*180/pi
    @test θ_GMST ≈ 152.578787810 atol = 1e-4
end
