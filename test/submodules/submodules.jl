# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   This files must call all the tests related to the submodules.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

@testset "SGP4" begin
    cd("./SGP4")
    include("./SGP4/sgp4.jl")
    cd("../")
end

@testset "SatelliteToolboxTLE" begin
    cd("./SatelliteToolboxTLE")
    include("./SatelliteToolboxTLE/tle.jl")
    cd("../")
end
