# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Algorithm test for SGP4. All tests are based on [1].
#
#   Notice that only the tests related to the SGP4 in [1] are considered,
#   because the SDP4 (orbit propagator for deep space) are not available yet in
#   SatelliteToolbox.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S (2006). Revisiting
#       Spacetrack Report #3: Rev1. AIAA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################################################################
#                                   Test 01
################################################################################

@testset "Tests from the paper AIAA 2006-6753" begin
    # Read all TLEs that will be used to test.
    tles = read_tle("./sgp4_tests/sgp4_tests.tle")

    st_sgp4_result = []
    SGP4_results = []

    for tle in tles
        filename = @sprintf(
            "./sgp4_tests/aiaa-2006-6753/sgp4_tle_%d_result.txt",
            tle.sat_num
        )
        SGP4_results = readdlm(filename; comments = true)

        # Initialize the orbit propagator.
        sgp4d = sgp4_init(tle, sgp4c_wgs72)

        t = SGP4_results[:,1]
        @inbounds for k = 1:length(t)

            # Propagate the orbit.
            r_TEME, v_TEME = sgp4!(sgp4d, t[k])

            # Assemble the result vector.
            st_sgp4_result = vcat(t[k], r_TEME, v_TEME)

            # Compare the values.
            @test st_sgp4_result[1] == SGP4_results[k,1]
            @test st_sgp4_result[2] ≈  SGP4_results[k,2] atol=1e-8
            @test st_sgp4_result[3] ≈  SGP4_results[k,3] atol=1e-8
            @test st_sgp4_result[4] ≈  SGP4_results[k,4] atol=1e-8
            @test st_sgp4_result[5] ≈  SGP4_results[k,5] atol=1e-9
            @test st_sgp4_result[6] ≈  SGP4_results[k,6] atol=1e-9
            @test st_sgp4_result[7] ≈  SGP4_results[k,7] atol=1e-9
        end
    end
end
