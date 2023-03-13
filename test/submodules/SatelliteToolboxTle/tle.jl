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

    date_1 = jd_to_date(tles[1].epoch)
    date_2 = jd_to_date(tles[2].epoch)

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
        true)
end

# Test conversion to string
# =========================

@testset "Conversion from TLE to string" begin

    # Conversion to string
    # ==========================================================================

    tles = read_tle("./tles_20200122.tle")
    f    = open("./tles_20200122.tle", "r")

    for i = 1:length(tles)
        stri_name = readline(f, keep = true)
        stri_l1   = readline(f, keep = true)
        stri_l2   = readline(f, keep = false)

        # The TLE for AMOS-4 has a `-` before the first time
        # derivative, even though it is 0. Since this is happening only here,
        # we will skip this case.
        tles[i].sat_num == 39237 && continue

        # The conversion of the exponent signal of BSTAR does not have a
        # defined pattern if BSTAR = 0.
        bstar_exp_le = false
        if tles[i].bstar == 0
            bstar_exp_le = stri_l1[60] == '-'
        end

        stri = stri_name * stri_l1 * stri_l2

        # If the OS is Windows, then we should remove `\r` to avoid testing
        # failure.
        Sys.iswindows() && (stri = replace(stri, "\r" => ""))

        strf = tle_to_str(tles[i], bstar_exp_le = bstar_exp_le)

        @test strf == stri
    end

    close(f)

    # Printing to IO
    # ==========================================================================

    str_tle_cbers4 = """
    CBERS 4                 
    1 40336U 14079A   18166.15595376 -.00000014  00000-0  10174-4 0  9993
    2 40336  98.4141 237.7928 0001694  75.7582 284.3804 14.35485112184485"""

    tle_cbers4 = read_tle_from_string(str_tle_cbers4)[1]

    tle_wcs = TLE(tle_cbers4; checksum_l1 = 0, checksum_l2 = 1)

    # Print without modifying the checksum.
    io = IOBuffer()
    print_tle(io, tle_wcs, recompute_checksum = false)
    result = String(take!(io))

    expected = """
    CBERS 4                 
    1 40336U 14079A   18166.15595376 -.00000014  00000-0  10174-4 0  9990
    2 40336  98.4141 237.7928 0001694  75.7582 284.3804 14.35485112184481"""

    @test result == expected

    # Print recomputing the checksum.
    io = IOBuffer()
    print_tle(io, tle_wcs)
    result = String(take!(io))
    expected = str_tle_cbers4

    @test result == expected
end
