# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to ground repeating orbits.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/orbit/ground_repeating_orbits.jl
# ============================================

# Function ground_repeating_orbit_adjacent_track_angle
# ----------------------------------------------------

################################################################################
#                                 Test results
################################################################################
#
# Amazonia-1 orbit
# ==============================================================================
#
#   Altitude above Equator : 752.845 km
#   Orbit period           : 100 min.
#   Inclination            : 98.405 °
#   Orbit period           : 5 days
#   Minimum FOV            : 19.9361 °
#
################################################################################

@testset "Function ground_repeating_orbit_adjacent_track_angle" begin
    h = 752.845e3
    p = 6000
    i = 98.405 |> deg2rad
    c = 5
    min_fov = ground_repeating_orbit_adjacent_track_angle(h, p, i, c)

    @test (min_fov |> rad2deg) ≈ 19.9361 (atol = 1e-4)

    # Test return types.
    min_fov = ground_repeating_orbit_adjacent_track_angle(h, p, i, c)
    @test min_fov isa Float64

    T = Float32
    min_fov = ground_repeating_orbit_adjacent_track_angle(T(h), T(p), T(i), c)
    @test min_fov isa T
end

# Function ground_repeating_orbit_adjacent_track_distance
# -------------------------------------------------------

################################################################################
#                                 Test results
################################################################################
#
# Amazonia-1 orbit
# ==============================================================================
#
#   Altitude above Equator : 752.845 km
#   Orbit period           : 100 min.
#   Inclination            : 98.405 °
#   Orbit period           : 5 days
#   Track distance         : 550.604 m
#
################################################################################

@testset "Function ground_repeating_orbit_adjacent_track_distance" begin
    p = 6000
    i = 98.405 |> deg2rad
    c = 5
    d = ground_repeating_orbit_adjacent_track_distance(p, i, c)

    @test d / 1000 ≈ 550.604 (atol = 1e-3)

    # Test return types.
    d = ground_repeating_orbit_adjacent_track_distance(p, i, c)
    @test d isa Float64

    T = Float32
    d = ground_repeating_orbit_adjacent_track_distance(T(p), T(i), c)
    @test d isa T
end
