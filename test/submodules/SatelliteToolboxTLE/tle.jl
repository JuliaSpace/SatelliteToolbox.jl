#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to TLE parser.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# Macros tle_str and tlenc_str
# ============================

@testset "Macros tle_str and tlenc_str" begin
    # Read the SCDs TLE from the file.
    tles_file = read_tle("./SCDs.tle")

    # Read the same TLEs from a string.
    tles_str = tle"""
        SCD 1
        1 22490U 93009B   18165.62596833  .00000225  00000-0  11410-4 0  9991
        2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338066

        SCD 2
        1 25504U 98060A   18165.15074951  .00000201  00000-0  55356-5 0  9994
        2 25504  24.9961  80.1303 0017060 224.4822 286.6438 14.44043397 37312
        """

    # Read the same TLES from a string with wrong checksums.
    # This should not output any exceptions.
    tles_str_nc = tlenc"""
        SCD 1
        1 22490U 93009B   18165.62596833  .00000225  00000-0  11410-4 0  9990
        2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338060

        SCD 2
        1 25504U 98060A   18165.15074951  .00000201  00000-0  55356-5 0  9990
        2 25504  24.9961  80.1303 0017060 224.4822 286.6438 14.44043397 37310
        """

    # Compare the TLEs.
    @test length(tles_file) == length(tles_str)
    @test length(tles_file) == length(tles_str_nc)

    for i = 1:length(tles_file)
        for sym in fieldnames(TLE)
            @test getfield(tles_file[i], sym) == getfield(tles_str[i], sym)

            # Skip the comparison of checksums for `tles_str_nc`.
            ( (sym == :checksum_l1) || (sym == :checksum_l2) ) && continue
            @test getfield(tles_file[i], sym) == getfield(tles_str_nc[i], sym)
        end
    end
end

# Function: read_tle_from_string
# ==============================

@testset "Function read_tle_from_string" begin
    # Read the SCD 1 TLE from the file.
    tle_scd1_expected = read_tle("./SCDs.tle")[1]

    # Parse the same TLE using the function `read_tle_from_string`.
    tle_scd1_result = read_tle_from_string(
        "1 22490U 93009B   18165.62596833  .00000225  00000-0  11410-4 0  9991",
        "2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338066")[1]

    # Compare the TLEs.
    for sym in fieldnames(TLE)
        # Skip name and satellite number.
        ( (sym == :name) || (sym == :sat_num) ) && continue

        @test getfield(tle_scd1_result, sym) == getfield(tle_scd1_expected, sym)
    end
end

# Date conversion
# ===============

@testset "Epoch year conversion" begin
    tles = tlenc"""
    1 22490U 93009B   75165.62596833  .00000225  00000-0  11410-4 0  9990
    2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338060

    1 22490U 93009B   76165.62596833  .00000225  00000-0  11410-4 0  9990
    2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338060
    """

    date_1 = JDtoDate(tles[1].epoch)
    date_2 = JDtoDate(tles[2].epoch)

    @test date_1[1] == 2075
    @test date_2[1] == 1976
end

# Test errors
# ===========

@testset "TLE format errors" begin

    @test_throws ErrorException read_tle_from_string("""
    SCD 1
    1 22490U 93009B   18165.62596833  .00000225  00000-0  11410-4 0  9991
    2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338066
    3 123
    """)

    # Errors in line 1
    # --------------------------------------------------------------------------

    @test_throws ErrorException read_tle_from_string("""
    SCD 1
    1 022490U 93009B   18165.62596833  .00000225  00000-0  11410-4 0  9991
    2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338066
    """)

    @test_throws ErrorException read_tle_from_string("""
    SCD 1
    1 022490U 93009B   181656.2596833  .00000225  00000-0  11410-4 0  9991
    2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338066
    """)

    # TODO: Errors in line 2.
end

@testset "Checksum errors" begin
    @test_throws ErrorException read_tle("./SCDs-wrong_checksum.tle")
    @test_throws LoadError @eval tle"""
    SCD 1
    1 22490U 93009B   18165.62596833  .00000225  00000-0  11410-4 0  9991
    2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338066

    SCD 2
    1 25504U 98060A   18165.15074951  .00000201  00000-0  55356-5 0  9994
    2 25504  24.9961  80.1303 0017060 224.4822 286.6438 14.44043397 37313
    """
    @test_throws ErrorException read_tle_from_string(
        "1 22490U 93009B   18165.62596833  .00000225  00000-0  11410-4 0  9991",
        "2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338063",
        verify_checksum = true)
end
