module SatToolbox

export SAT_LIGHTING_SUNLIGHT, SAT_LIGHTING_UMBRA, SAT_LIGHTING_PENUMBRA
export sat_lighting_condition, satellite_position_i, sun_position_i
export OrbitalParameters

import Base: asin, atan2, cos, mod, sin

type OrbitalParameters
    a::Float64
    e::Float64
    i::Float64
    RAAN::Float64
    w::Float64
    f::Float64
end

include("satellite_lighting_conditions.jl")
include("satellite_position.jl")
include("sun_position.jl")

end # module
