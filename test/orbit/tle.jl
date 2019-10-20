#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to TLE parser.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/orbit/tle.jl
# ========================

# Macros tle_str and tlenc_str
# ----------------------------

@testset "Macros tle_str and tlenc_str" begin
    # Read the SCDs TLE from the file.
    tles_file = SGP4.read_tle("./SCDs.tle")

    # Read the same TLEs from a string.
    tles_str = SGP4.tle"""
        SCD 1
        1 22490U 93009B   18165.62596833  .00000225  00000-0  11410-4 0  9991
        2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338066

        SCD 2
        1 25504U 98060A   18165.15074951  .00000201  00000-0  55356-5 0  9994
        2 25504  24.9961  80.1303 0017060 224.4822 286.6438 14.44043397 37312
        """

    # Read the same TLES from a string with wrong checksums.
    # This should not output any exceptions.
    tles_str_nc = SGP4.tlenc"""
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
        for sym in fieldnames(SGP4.TLE)
            @test getfield(tles_file[i], sym) == getfield(tles_str[i], sym)

            # Skip the comparison of checksums for `tles_str_nc`.
            ( (sym == :checksum_l1) || (sym == :checksum_l2) ) && continue
            @test getfield(tles_file[i], sym) == getfield(tles_str_nc[i], sym)
        end
    end
end
