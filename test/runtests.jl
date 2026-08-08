using Test

using Makie
using SatelliteToolbox

@testset "Makie Extension" verbose = true begin
    include("./makie.jl")
end

@testset "Orbit" verbose = true begin
    include("./orbit/general.jl")
end

@testset "Time" verbose = true begin
    include("./time/equation_of_time.jl")
    include("./time/raan.jl")
end
