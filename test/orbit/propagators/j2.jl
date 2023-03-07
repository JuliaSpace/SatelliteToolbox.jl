# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests of the J2 orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Source: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/91
#
# The STK provided the following result for the J2Perturbation orbit propagator.
#
# Initial mean orbital elements:
#
#            Orbit epoch : 2023-01-01T00:00:00.000
#        Semi-major axis : 8000.000 km
#           Eccentricity : 0.015
#            Inclination : 28.5°
#                   RAAN : 100°
#    Argument of perigee : 200°
#   Initial mean anomaly : 45°
#
#
# State vector on 2023-01-05T00:00:00.000
#
#   r = [-6849.654348, -2253.059809, 3574.529667] km
#   v = [    1.656142,    -6.699518,   -1.266334] km/s
#
################################################################################

@testset "J2 orbit propagator" begin
    jd₀ = date_to_jd(2023, 1, 1, 0, 0, 0)
    jd₁ = date_to_jd(2023, 1, 5, 0, 0, 0)

    # Float64
    # ==========================================================================

    let T = Float64
        orb = KeplerianElements(
            jd₀,
            T(8000e3),
            T(0.015),
            T(28.5) |>deg2rad,
            T(100) |>deg2rad,
            T(200) |>deg2rad,
            T(45) |> deg2rad
        )

        orbp = init_orbit_propagator(Val(:J2), orb; j2_gc = j2_gc_egm08)
        r, v = propagate_to_epoch!(orbp, jd₁)

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-3
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-3
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-3
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-5
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-5
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-5
        @test eltype(v) == T
    end

    let T = Float32
        orb = KeplerianElements(
            jd₀,
            T(8000e3),
            T(0.015),
            T(28.5) |>deg2rad,
            T(100) |>deg2rad,
            T(200) |>deg2rad,
            T(45) |> deg2rad
        )

        orbp = init_orbit_propagator(Val(:J2), orb; j2_gc = j2_gc_egm08_f32)
        r, v = propagate_to_epoch!(orbp, jd₁)

        @test r[1] / 1000 ≈ -6849.654348 atol = 5e-1
        @test r[2] / 1000 ≈ -2253.059809 atol = 5e-1
        @test r[3] / 1000 ≈ +3574.529667 atol = 5e-1
        @test eltype(r) == T

        @test v[1] / 1000 ≈ +1.656142 atol = 1e-3
        @test v[2] / 1000 ≈ -6.699518 atol = 1e-3
        @test v[3] / 1000 ≈ -1.266334 atol = 1e-3
        @test eltype(v) == T
    end
end
