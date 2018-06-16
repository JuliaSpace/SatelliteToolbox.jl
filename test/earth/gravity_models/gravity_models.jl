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
#   Tests related to the gravity models.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] http://icgem.gfz-potsdam.de/home
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-06-16: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
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
