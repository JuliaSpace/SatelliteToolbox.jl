#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to the Jacchia-Bowman 2008 model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] http://sol.spacenvironment.net/~JB2008/
#   [2] http://sol.spacenvironment.net/~JB2006/
#   [3] https://github.com/JuliaSpace/JB2008_Test
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/earth/atmospheric_models/jb2008
# ===========================================

# Functions: jb2008
# -----------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
#   Compare the results with the tests in files:
#
#       JB_AUTO_OUTPUT_01.DAT
#       JB_AUTO_OUTPUT_02.DAT
#       JB_AUTO_OUTPUT_03.DAT
#       JB_AUTO_OUTPUT_04.DAT
#
#   which were created using JuliaSpace/JB2008_Test@1552551 .
#
################################################################################

@testset "Function jb2008" begin
    # Scenario 01
    # ===========

    # Initialize the space indices with the local files.
    init_space_indices(dtcfile_path = "./DTCFILE_DEMO.TXT",
                       solfsmy_path = "./SOLFSMY_DEMO.TXT")

    # Files with the test results.
    test_list = ["JB2008_AUTO_OUTPUT_01.DAT",
                 "JB2008_AUTO_OUTPUT_02.DAT",
                 "JB2008_AUTO_OUTPUT_03.DAT",
                 "JB2008_AUTO_OUTPUT_04.DAT"]

    # Execute the tests.
    for filename in test_list
        open(filename) do file
            line_num = 0
            year     = 0
            doy      = 0
            hour     = 0
            min      = 0
            sec      = 0.0

            for line in eachline(file)
                line_num += 1

                # Ignore 5 lines to skip the header.
                (line_num <= 5) && continue

                # Ignore others non important lines.
                (line_num in [7, 8, 9]) && continue

                if line_num == 6
                    # Read the next line to obtain the input data related to the
                    # time.
                    tokens = split(line)

                    year = parse(Int64,   tokens[1])
                    doy  = parse(Int64,   tokens[2])
                    hour = parse(Int64,   tokens[3])
                    min  = parse(Int64,   tokens[4])
                    sec  = parse(Float64, tokens[5])
                else
                    tokens = split(line)

                    # Read the information about the location.
                    h    = parse(Float64, tokens[1])*1000
                    glon = parse(Float64, tokens[2])*pi/180
                    glat = parse(Float64, tokens[3])*pi/180

                    # Read the model output.
                    T_exo = parse(Float64, tokens[4])
                    Tz    = parse(Float64, tokens[5])
                    ρ     = parse(Float64, tokens[6])

                    # Run the model.
                    JD  = DatetoJD(year, 1, 1, hour, min, sec) - 1 + doy
                    out = jb2008(JD, glat, glon, h)

                    # Compare the results.
                    @test out.T_exo ≈ T_exo atol = 0.6 rtol = 0.0
                    @test out.Tz    ≈ Tz    atol = 0.6 rtol = 0.0
                    @test out.rho   ≈ ρ     atol = 0.0 rtol = 5e-3
                end
            end
        end
    end
end
