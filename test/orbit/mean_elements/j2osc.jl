# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Test conversion from osculating to mean elements using J2 osc. orbit
#   propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/orbit/mean_elements/j2osc.jl
# ========================================

# Function rv_to_mean_elements_j2osc
# ----------------------------------

@testset "Function rv_to_mean_elements_j2osc" begin
    # We just need to run the J2 osc. orbit propagator, obtain the osculating
    # elements, convert them to mean elements, and compare with the original
    # ones.
    jd₀ = date_to_jd(2021, 6, 19, 19, 35, 35)
    orb_input = KeplerianElements(
        jd₀,
        7130.982e3,
        0.001111,
        deg2rad(98.405),
        compute_RAAN_lt(jd₀, 22.5),
        deg2rad(90),
        deg2rad(10)
    )


    # Generate the osculating elements.
    orbp   = init_orbit_propagator(Val(:J2osc), orb_input)
    vr, vv = propagate!(orbp, 0:10:12000)
    vjd    = get_epoch(orbp) .+ (0:10:12000) ./ 86400

    # Obtain the mean elements.
    orb, ~ = rv_to_mean_elements_j2osc(
        vjd,
        vr,
        vv;
        mean_elements_epoch = :begin,
        print_debug         = false,
        max_iterations      = 1000,
        atol                = 1e-10,
        rtol                = 1e-10
    )

    # Compare.
    @test orb.t == orb_input.t
    @test orb.a ≈  orb_input.a atol = 1e-4
    @test orb.e ≈  orb_input.e atol = 1e-11
    @test orb.i ≈  orb_input.i atol = 1e-10
    @test orb.Ω ≈  orb_input.Ω atol = 1e-10
    @test orb.ω ≈  orb_input.ω atol = 1e-7
    @test orb.f ≈  orb_input.f atol = 1e-7
end
