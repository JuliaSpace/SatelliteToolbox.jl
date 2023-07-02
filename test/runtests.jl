using Test

using SatelliteToolbox

@testset "Orbit" verbose = true begin
    include("./orbit/general.jl")
end

@testset "Time" verbose = true begin
     include("./time/equation_of_time.jl")
     include("./time/raan.jl")
end
