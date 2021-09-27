# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Test conversion from osculating to mean elements using SGP4.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/orbit/mean_elements/sgp4.jl
# =======================================

# Function rv_to_mean_elements_sgp4
# ---------------------------------

@testset "Function rv_to_mean_elements_sgp4" begin
    # We just need to run the SGP4, obtain the osculating elements, convert them
    # to mean elements, and compare with the original TLE.

    tle = tle"""
    AMAZONIA 1
    1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
    2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436"""[1]

    # Generate the osculating elements (TEME).
    orbp   = init_orbit_propagator(Val(:sgp4), tle)
    vr, vv = propagate!(orbp, 0:10:12000)
    vjd    = get_epoch(orbp) .+ (0:10:12000) ./ 86400

    # Obtain the mean elements.
    epoch, xo, ~ = rv_to_mean_elements_sgp4(
        vjd,
        vr,
        vv;
        mean_elements_epoch = :begin,
        print_debug         = false,
        max_iterations      = 1000,
        atol                = 1e-10,
        rtol                = 1e-10
    )

    # Convert some elements to compare with the TLE.
    a₀ = xo[1] / (1000sgp4_gc_wgs84.R0)
    e₀ = xo[2]
    i₀ = rad2deg(xo[3])
    Ω₀ = rad2deg(xo[4])
    ω₀ = rad2deg(xo[5])
    M₀ = rad2deg(f_to_M(e₀, xo[6]))

    # Obtain the mean motion [rev / day].
    n₀ = 720 * sgp4_gc_wgs84.XKE / (sqrt(a₀ * a₀ * a₀) * π)

    # Compare.
    @test epoch == tle.epoch
    @test n₀    ≈  tle.n     atol = 1e-7
    @test e₀    ≈  tle.e     atol = 1e-7
    @test i₀    ≈  tle.i     atol = 1e-5
    @test Ω₀    ≈  tle.Ω     atol = 1e-5
    @test ω₀    ≈  tle.ω     atol = 1e-5
    @test M₀    ≈  tle.M     atol = 1e-5
    @test xo[7] ≈  tle.bstar atol = 1e-6
end

# Function rv_to_tle
# ------------------

@testset "Function rv_to_tle" begin
    # We just need to run the SGP4, obtain the osculating elements, convert them
    # to TLE, and compare with the original TLE.

    tle_input = tle"""
    AMAZONIA 1
    1 47699U 21015A   21270.48626105 -.00000044  00000-0  19860-2 0  9993
    2 47699  98.4889 344.6059 0001597  74.4244 285.7135 14.40801240 30436"""[1]

    # Generate the osculating elements (TEME).
    orbp   = init_orbit_propagator(Val(:sgp4), tle_input)
    vr, vv = propagate!(orbp, 0:10:12000)
    vjd    = get_epoch(orbp) .+ (0:10:12000) ./ 86400

    # Obtain the mean elements.
    tle, ~ = rv_to_tle(
        vjd,
        vr,
        vv;
        elem_set_number     = 999,
        int_designator      = "21015A  ",
        mean_elements_epoch = :begin,
        name                = "AMAZONIA 1",
        rev_num             = 3043,
        sat_num             = 47699,
        print_debug         = false,
        max_iterations      = 1000,
        atol                = 1e-10,
        rtol                = 1e-10
    )

    # Compare.
    @test tle.classification  == tle_input.classification
    @test tle.elem_set_number == tle_input.elem_set_number
    @test tle.epoch           == tle_input.epoch
    @test tle.epoch_year      == tle_input.epoch_year
    @test tle.int_designator  == tle_input.int_designator
    @test tle.name            == tle_input.name
    @test tle.rev_num         == tle_input.rev_num
    @test tle.sat_num         == tle_input.sat_num

    @test tle.bstar     ≈  tle_input.bstar     atol = 1e-6
    @test tle.e         ≈  tle_input.e         atol = 1e-7
    @test tle.epoch_day ≈  tle_input.epoch_day atol = 1e-8
    @test tle.i         ≈  tle_input.i         atol = 1e-4
    @test tle.M         ≈  tle_input.M         atol = 1e-4
    @test tle.n         ≈  tle_input.n         atol = 1e-7
    @test tle.Ω         ≈  tle_input.Ω         atol = 1e-4
    @test tle.ω         ≈  tle_input.ω         atol = 1e-4
end
