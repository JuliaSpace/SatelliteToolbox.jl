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
#   Algorithm test for SGP4. All tests are based on [1].
#
#   Notice that only the tests related to the SGP4 in [1] are considered,
#   because the SDP4 (orbit propagator for deep space) are not available yet in
#   SatToolbox.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S (2006). Revisiting
#       Spacetrack Report #3: Rev1. AIAA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-04-02: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

println("$(c)SGP4 Tests$d")
println("$(c)----------$d")
println("")

################################################################################
#                                   Test 01
################################################################################

println("    $(b)Test 01 - Tests from the paper AIAA 2006-6753$d")
println("")

# Read all TLEs that will be used to test.
tles = read_tle("sgp4_tests.tle")

# Array with the satellite numbers that must be tested.
sat_num_array = [5; 6251; 28057; 28350; 28872; 29141; 29238; 88888]

# Find the TLEs with the respective satellite numbers.
tles_test_01 = tles[find( (tle->tle.sat_num in sat_num_array), tles)]

st_sgp4_result = []
SGP4_results = []

for tle in tles_test_01
    println("    Testing the following TLE:")
    println("")
    println(tle)

    filename = @sprintf("./aiaa-2006-6753/sgp4_tle_%d_result.txt", tle.sat_num)
    SGP4_results = readdlm(filename)

    # Initialize the orbit propagator.
    orbp = init_orbit_propagator(Val{:sgp4}, tle, sgp4_gc_wgs72)

    # Propagate the orbit.
    t = SGP4_results[:,1]*60
    (orbm, r_TEME, v_TEME) = propagate!(orbp, t)

    # Compare the results.
    for k = 1:length(t)
        # Assemble the result vector.
        st_sgp4_result = [(orbm[k].t-orbm[1].t)/60 r_TEME[k]'/1000 v_TEME[k]'/1000]

        # Compare the values.
        @test st_sgp4_result[1] == SGP4_results[k,1]
        @test st_sgp4_result[2] ≈  SGP4_results[k,2] atol=1e-8
        @test st_sgp4_result[3] ≈  SGP4_results[k,3] atol=1e-8
        @test st_sgp4_result[4] ≈  SGP4_results[k,4] atol=1e-8
        @test st_sgp4_result[5] ≈  SGP4_results[k,5] atol=1e-9
        @test st_sgp4_result[6] ≈  SGP4_results[k,6] atol=1e-9
        @test st_sgp4_result[7] ≈  SGP4_results[k,7] atol=1e-9
    end

    println("    $(b)Test passed!$d")
    println("")
end
