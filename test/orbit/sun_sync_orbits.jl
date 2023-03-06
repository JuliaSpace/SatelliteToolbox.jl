# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to Sun synchronous orbit functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/orbit/syn_sync_orbits.jl
# ====================================

# Function syn_sync_orbit_from_ang_vel
# ------------------------------------

################################################################################
#                                 Test results
################################################################################
#
# Amazonia-1 orbit
# ==============================================================================
#
#   Semi-major axis : 7130.982 km
#   Eccentricity    : 0.001111
#   Inclination     : 98.405°
#   Period          : 100 min.
#
################################################################################

@testset "Function sun_sync_orbit_from_ang_vel" begin
    # Angular velocity of Amazonia-1.
    n = 2π / 6000

    a, i = sun_sync_orbit_from_ang_vel(n, 0.001111)

    @test a / 1000 ≈ 7130.982 (atol = 1e-3)
    @test i * 180 / π ≈ 98.405 (atol = 1e-3)

    # Test return types.
    a, i = sun_sync_orbit_from_ang_vel(n)
    @test a isa Float64
    @test i isa Float64

    a, i = sun_sync_orbit_from_ang_vel(n, 0.0)
    @test a isa Float64
    @test i isa Float64

    a, i = sun_sync_orbit_from_ang_vel(Float32(n))
    @test a isa Float32
    @test i isa Float32

    a, i = sun_sync_orbit_from_ang_vel(Float32(n), 0.0f0)
    @test a isa Float32
    @test i isa Float32

    # Test errors.
    @test_throws ArgumentError sun_sync_orbit_from_ang_vel(-1.0)
    @test_throws ArgumentError sun_sync_orbit_from_ang_vel(n, -0.02)
    @test_throws ArgumentError sun_sync_orbit_from_ang_vel(n, 1.0)

    # Test warnings.
    @test_logs (:warn,) sun_sync_orbit_from_ang_vel(n; max_iterations = 2)
end

# Function sun_sync_orbit_semi_major_axis
# ---------------------------------------

################################################################################
#                                 Test results
################################################################################
#
# Amazonia-1 orbit
# ==============================================================================
#
#   Semi-major axis : 7130.982 km
#   Eccentricity    : 0.001111
#   Inclination     : 98.405°
#
################################################################################

@testset "Function sun_sync_orbit_semi_major_axis" begin
    # Amazonia-1 inclination.
    i = 98.40547 |> deg2rad

    a = sun_sync_orbit_semi_major_axis(i, 0.001111)
    @test a / 1000 ≈ 7130.982 (atol = 1e-3)

    # Test return types.
    a = sun_sync_orbit_semi_major_axis(i)
    @test a isa Float64

    a = sun_sync_orbit_semi_major_axis(i, 0.001111)
    @test a isa Float64

    a = sun_sync_orbit_semi_major_axis(Float32(i))
    @test a isa Float32

    a = sun_sync_orbit_semi_major_axis(Float32(i), 0.001111f0)
    @test a isa Float32

    # Test errors.
    @test_throws ArgumentError sun_sync_orbit_semi_major_axis(i, -0.02)
    @test_throws ArgumentError sun_sync_orbit_semi_major_axis(i, 1.0)
    @test_throws ErrorException sun_sync_orbit_semi_major_axis(91 |> deg2rad)
end

# Function sun_sync_orbit_inclination
# -----------------------------------

@testset "Function sun_sync_orbit_inclination" begin
    # Amazonia-1 semi-major axis.
    a = 7130.982e3

    i = sun_sync_orbit_inclination(a, 0.001111)
    @test i * 180 / π ≈ 98.405 (atol = 1e-3)

    # Test return types.
    i = sun_sync_orbit_inclination(a)
    @test i isa Float64

    i = sun_sync_orbit_inclination(a, 0.001111)
    @test i isa Float64

    i = sun_sync_orbit_inclination(Float32(a))
    @test i isa Float32

    i = sun_sync_orbit_inclination(Float32(a), 0.001111f0)
    @test i isa Float32

    # Test errors.
    @test_throws ArgumentError sun_sync_orbit_inclination(5000e3)
    @test_throws ArgumentError sun_sync_orbit_inclination(a, -0.02)
    @test_throws ArgumentError sun_sync_orbit_inclination(a, 1.0)
    @test_throws ErrorException sun_sync_orbit_inclination(16_500e3)
end
