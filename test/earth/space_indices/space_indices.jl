#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to the space indices interface.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] http://sol.spacenvironment.net/jb2008/indices/DTCFILE.TXT
#   [2] http://sol.spacenvironment.net/jb2008/indices/SOLFSMY.TXT
#   [3] https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/earth/space_indices/*
# =================================

# Function: get_space_indices
# ---------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Space indices at 2017-10-19 @ 06:30:00
#
#   F10.7 (observed)                     = 73.4 (from SOLFSMY.TXT)
#   F10.7 (adjusted)                     = 72.8 (from online NRLMSISE-00)
#   F10.7 (observed) 81-day average mean = 76.5 (from SOLFSMY.TXT)
#   F10.7 (adjusted) 90-day average mean = 77.8 (from online NRLMSISE-00)
#   S10                                  = 63.6 (from SOLFSMY.TXT)
#   S10 81-day average mean              = 64.9 (from SOLFSMY.TXT)
#   M10                                  = 72.7 (from SOLFSMY.TXT)
#   M10 81-day average mean              = 78.4 (from SOLFSMY.TXT)
#   Y10                                  = 82.4 (from SOLFSMY.TXT)
#   Y10 81-day average mean              = 83.6 (from SOLFSMY.TXT)
#   DstΔTc                               = 67.5 (from DTCFILE.TXT)
#
################################################################################

@testset "Function get_space_indices" begin
    # Init the space indices using the local files.
    init_space_indices(dtcfile_path = "./files/DTCFILE.TXT",
                       fluxtable_path = "./files/fluxtable.txt",
                       solfsmy_path = "./files/SOLFSMY.TXT",
                       wdcfiles_dir = "./files/",
                       wdcfiles_oldest_year = 2017,
                       wdcfiles_newest_year = 2018)

    # Get the indices.
    JD       = DatetoJD(2017,10,19,6,30,0)
    vF10adj  = get_space_index(F10()    , JD)
    vF10obs  = get_space_index(F10obs() , JD)
    vF10Madj = get_space_index(F10M()   , JD; window = 90)
    vF10Mobs = get_space_index(F10Mobs(), JD)
    vS10     = get_space_index(S10()    , JD)
    vS81a    = get_space_index(S81a()   , JD)
    vM10     = get_space_index(M10()    , JD)
    vM81a    = get_space_index(M81a()   , JD)
    vY10     = get_space_index(Y10()    , JD)
    vY81a    = get_space_index(Y81a()   , JD)
    vDstΔTc  = get_space_index(DstΔTc() , JD)

    # Test.
    @test vF10obs  ≈ 73.4 atol = 1e-2
    @test vF10adj  ≈ 72.8 atol = 1e-2
    @test vF10Mobs ≈ 76.5 atol = 1e-2
    # TODO: The adjusted mean is not exatcly what the online NRLMSISE00 is
    # computing.
    @test vF10Madj ≈ 77.8 atol = 5e-1
    @test vS10     ≈ 63.6 atol = 1e-2
    @test vS81a    ≈ 64.9 atol = 1e-2
    @test vM10     ≈ 72.7 atol = 1e-2
    @test vM81a    ≈ 78.4 atol = 1e-2
    @test vY10     ≈ 82.4 atol = 1e-2
    @test vY81a    ≈ 83.6 atol = 1e-2
end
