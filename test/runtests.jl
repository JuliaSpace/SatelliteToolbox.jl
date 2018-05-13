using SatelliteToolbox
using ReferenceFrameRotations
using Compat.Test
using TestSetExtensions

# Colors

b = "\x1b[1m"
d = "\x1b[0m"
g = "\x1b[1m\x1b[32m"
y = "\x1b[1m\x1b[33m"
c = "\x1b[1m\x1b[36m"

@testset "General orbit functions" begin
    include("./orbit/general.jl")
end
println("")

@testset "Orbit propagators" begin
    cd("./orbit/propagators/")
    include("./orbit/propagators/twobody.jl")
    include("./orbit/propagators/sgp4.jl")
    cd("../../")
end
println("")

@testset "Coordinate transformations" begin
    cd("./transformations/fk5/")
    include("./transformations/fk5/fk5.jl")
    cd("../../")
    cd("./transformations/")
    include("./transformations/ecef_to_eci.jl")
    include("./transformations/eci_to_ecef.jl")
    cd("../")
end
println("")

@testset "Functions related with time" begin
    include("./time/time.jl")
end
println("")

