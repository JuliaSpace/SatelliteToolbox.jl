#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to the exponential atmospheric model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/earth/atmospheric_models/expatmosphere
# ==================================================

# Functions: expatmosphere
# -------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
#   Example 8-4 [1, p. 567]
#
#   Using the exponential atmospheric model, one gets:
#
#       ρ( 747.2119 km ) = 2.1219854 ⋅ 10⁻⁴ kg / m³
#
# Scenario 02
# ===========
#
#   By definition, one gets `ρ(h₀) = ρ₀(h₀)` for all tabulated values.
#
# Scenario 03
# ===========
#
#   Inside every interval, the atmospheric density must be monotonically
#   decreasing.
#
################################################################################

@testset "Function expatmosphere" begin
    # Scenario 01
    # ===========

    @test expatmosphere(747211.9) ≈ 2.1219854e-14 rtol = 1e-8

    # Scenario 02
    # ===========

    for i = 1:length(SatelliteToolbox._expatmosphere_h₀)
        h = SatelliteToolbox._expatmosphere_h₀[i]*1000

        @test expatmosphere(h) == SatelliteToolbox._expatmosphere_ρ₀[i]
    end

    # Scenario 03
    # ===========

    for i = 2:length(SatelliteToolbox._expatmosphere_h₀)
        h₀i = SatelliteToolbox._expatmosphere_h₀[i-1]*1000
        h₀f = SatelliteToolbox._expatmosphere_h₀[i]*1000
        Δ   = 100
        Δh₀ = (h₀f - h₀i)/Δ

        ρk = expatmosphere(h₀i)

        for k = 2:Δ
            ρk₋₁ = ρk
            h    = h₀i + Δh₀*(k-1)
            ρk   = expatmosphere(h)

            @test ρk < ρk₋₁
        end
    end
end
