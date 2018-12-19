#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to the transformations between Geodetic and Geocentric frames.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/transformations/geodetic_geocentric.jl
# ==================================================

# Function: ECEFtoGeodetic
# ------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-3: Converting ECEF to Lat Lon [1, p. 173].
#
# According to this example, using:
#
#   r_ecef = 6524.834 i + 6862.875 j + 6448.296 k [km]
#
# one gets:
#
#   Geodetic Latitude  = 34.352496°
#            Longitude = 46.4464°
#            Altitude  = 5085.22 km
#
# Scenario 02
# ===========
#
# At the poles, we have a singularity. We know that, if:
#
#   r_ecef = 0 i + 0 j + Z k
#
# then
#
#   Geodetic Latitude  = 90° for Z > 0 and -90° for Z < 0
#            Longitude =  0°
#            Altitude  = Z - b_wgs84
#
################################################################################

@testset "Function ECEFtoGeodetic" begin
    # Scenario 01
    # ===========

    r = [6524.834e3, 6862.875e3, 6448.296e3]

    ϕ_gd, λ_gd, h = ECEFtoGeodetic(r)

    @test rad2deg(ϕ_gd) ≈ 34.352496 atol = 1e-6
    @test rad2deg(λ_gd) ≈ 46.4464   atol = 1e-4
    @test h/1000        ≈ 5085.22   atol = 1e-2

    # Scenario 02
    # ===========

    aux = rand(0:1000)

    Z = R0 + aux
    ϕ_gd, λ_gd, h = ECEFtoGeodetic([0;0;Z])

    @test rad2deg(ϕ_gd) ≈ 90
    @test rad2deg(λ_gd) ≈ 0
    @test h             ≈ Z - SatelliteToolbox.b_wgs84

    Z = -R0 + aux
    ϕ_gd, λ_gd, h = ECEFtoGeodetic([0;0;Z])

    @test rad2deg(ϕ_gd) ≈ -90
    @test rad2deg(λ_gd) ≈ 0
    @test h             ≈ -Z - SatelliteToolbox.b_wgs84
end
