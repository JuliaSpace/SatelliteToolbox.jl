# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   This files must call all the tests related to the submodules.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

@testset "SatelliteToolboxSgp4" verbose = true begin
    cd("./SatelliteToolboxSgp4/")
    include("./SatelliteToolboxSgp4/sgp4.jl")
    cd("../")
end

@testset "SatelliteToolboxTle" verbose = true begin
    cd("./SatelliteToolboxTle")
    include("./SatelliteToolboxTle/tle.jl")
    cd("../")
end
