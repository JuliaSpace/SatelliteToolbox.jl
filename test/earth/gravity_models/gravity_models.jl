#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to the gravity models.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] http://icgem.gfz-potsdam.de/home
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/earth/gravity_models/gravity_models.jl
# ==================================================

# Function compute_g
# ------------------

################################################################################
#                                 TEST RESULTS
################################################################################
#
# For each embedded model, the norm of the gravitational field obtained from
# `compute_g` is compared to the online computation available in [1]. The inputs
# for the online computation were:
#
#   * Functional selection: `gravitation_ell`.
#   * Reference system: `WGS84`.
#   * Height over Ellipsoid: 0 m.
#   * Grid step: 10 deg.
#   * Tide System: Use unmodified model.
#
################################################################################

@testset "Function compute_g" begin
    # Read the inputs.
    tests_out = (readdlm("./test_results/gravitation/EGM96_6ddd3d7d77ec891a29f7c870494a8a2e215adae4256574e5de3a233380fbecee.gdf";
                         skipstart = 34),
                 readdlm("./test_results/gravitation/JGM2_435951880302b0efcb6c8ad16ce89d1840104941a3cfcd2ad0e25831a6c0f136.gdf";
                         skipstart = 34),
                 readdlm("./test_results/gravitation/JGM3_2d444404c0f11e585c91a810332530020c8ef0a1b51697aae6b9cbedcf3d67cc.gdf";
                         skipstart = 34))

    coefs = (load_gravity_model(EGM96()),
             load_gravity_model(JGM2()),
             load_gravity_model(JGM3()))

    # Compare the results.
    for m = 1:length(coefs)
        coef = coefs[m]
        test_out = tests_out[m]

        for k = 1:size(test_out,1)
            # Get latitude and longitude.
            lat = test_out[k,2]*pi/180
            lon = test_out[k,1]*pi/180

            # Get the expected result in m/s² (the online result is in mGal).
            norm_g_online = test_out[k,3]/100000

            # Use the model to compute the gravity using all the coefficients.
            norm_g = norm(compute_g(coef, GeodetictoECEF(lat, lon, 0)))

            # Compare the results.
            @test norm_g_online ≈ norm_g atol=5e-8 rtol=0
        end
    end
end

# Function compute_U
# ------------------

################################################################################
#                                 TEST RESULTS
################################################################################
#
# For each embedded model, the gravitational potential obtained from `compute_U`
# is compared to the online computation available in [1]. The inputs for the
# online computation were:
#
#   * Functional selection: `potential_ell`.
#   * Reference system: `WGS84`.
#   * Height over Ellipsoid: 0 m.
#   * Grid step: 10 deg.
#   * Tide System: Use unmodified model.
#
################################################################################

@testset "Function compute_U" begin
    # Read the inputs.
    tests_out = (readdlm("./test_results/potential/EGM96_08b822daf0755f2f048a0af2d1afa5dedc5afa5738e991071e7b0d12cf3337ae.gdf";
                         skipstart = 35),
                 readdlm("./test_results/potential/JGM2_7801c122125724d1f5022463d2f87728edca77158ec4af1643168e31ee3ea384.gdf";
                         skipstart = 35),
                 readdlm("./test_results/potential/JGM3_719cfffa3f04b93dce106013c2b745451b656bdd0550afa7f7d4700e03901f88.gdf";
                         skipstart = 35))

    coefs = (load_gravity_model(EGM96()),
             load_gravity_model(JGM2()),
             load_gravity_model(JGM3()))

    # Compare the results.
    for m = 1:length(coefs)
        coef = coefs[m]
        test_out = tests_out[m]

        for k = 1:size(test_out,1)
            # Get latitude and longitude.
            lat = test_out[k,2]*pi/180
            lon = test_out[k,1]*pi/180

            # Get the expected result in m²/s².
            U_online = test_out[k,3]

            # Use the model to compute the gravity using all the coefficients.
            U = compute_U(coef, GeodetictoECEF(lat, lon, 0))

            # Compare the results.
            @test U_online ≈ U atol=1e-6 rtol=0
        end
    end
end

# Issue #22 - Gravity model maximum degree
# ----------------------------------------

@testset "Issue #22 - Gravity model maximum degree" begin
    v  = (1 .+ rand(3))*R0
    gm = load_gravity_model(EGM96())

    g_0 = compute_g(gm, v, -1)
    g_1 = compute_g(gm, v, 1)
    g_2 = compute_g(gm, v, 2)
    g_3 = compute_g(gm, v, 360)
    g_4 = compute_g(gm, v, 360, 0)

    @test g_0 == g_3
    @test g_1 != g_2
    @test g_1 != g_3
    @test g_2 != g_3
    @test g_3 != g_4

    U_0 = compute_U(gm, v, -1)
    U_1 = compute_U(gm, v, 1)
    U_2 = compute_U(gm, v, 2)
    U_3 = compute_U(gm, v, 360)
    U_4 = compute_U(gm, v, 360, 0)

    @test U_0 == U_3
    @test U_1 != U_2
    @test U_1 != U_3
    @test U_2 != U_3
    @test U_3 != U_4
end
