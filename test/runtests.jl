using Test

using DelimitedFiles
using LinearAlgebra
using Printf
using ReferenceFrameRotations
using SatelliteToolbox

@testset "Atmospheric Models" begin
    cd("./earth/atmospheric_models/expatmosphere/")
    include("./earth/atmospheric_models/expatmosphere/expatmosphere.jl")
    cd("../../../")
    cd("./earth/atmospheric_models/nrlmsise00/")
    include("./earth/atmospheric_models/nrlmsise00/nrlmsise00.jl")
    cd("../../../")
    cd("./earth/atmospheric_models/jb2008/")
    include("./earth/atmospheric_models/jb2008/jb2008.jl")
    cd("../../../")
    cd("./earth/atmospheric_models/jr1971/")
    include("./earth/atmospheric_models/jr1971/jr1971.jl")
    cd("../../../")
end
println("")

@testset "Geomagnetic Field Models" begin
    cd("./earth/geomagnetic_field_models/")
    include("./earth/geomagnetic_field_models/igrf.jl")
    cd("../../")
end
println("")

@testset "Gravity Field Models" begin
    cd("./earth/gravity_models/")
    include("./earth/gravity_models/gravity_models.jl")
    cd("../../")
end
println("")

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

@testset "TLE Parser" begin
    cd("./orbit/")
    include("./orbit/tle.jl")
    cd("../")
end
println("")

@testset "Coordinate transformations" begin
    cd("./transformations/fk5/")
    include("./transformations/fk5/fk5.jl")
    cd("../../")
    cd("./transformations/teme/")
    include("./transformations/teme/teme.jl")
    cd("../../")
    cd("./transformations/")
    # These tests cannot pass with julia 0.7 due to the following issue of
    # Interpolations.jl:
    #
    #   https://github.com/JuliaMath/Interpolations.jl/issues/204
    #
    include("./transformations/ecef_to_ecef.jl")
    include("./transformations/ecef_to_eci.jl")
    include("./transformations/eci_to_ecef.jl")
    include("./transformations/eci_to_eci.jl")
    include("./transformations/orbit_elements.jl")
    cd("../")
end
println("")

@testset "Functions related with time" begin
     include("./time/time.jl")
end
println("")

@testset "Functions related with the Sun" begin
    include("./sun/equation_of_time.jl")
end
println("")
