module SatToolbox

export SAT_LIGHTING_SUNLIGHT, SAT_LIGHTING_UMBRA, SAT_LIGHTING_PENUMBRA
export satellite_lighting_condition, satellite_position_i, sun_position_i
export OrbitalParameters
export R0, m0, J2, Rs

import Base: asin, atan2, cos, mod, sin

type OrbitalParameters
    a::Float64
    e::Float64
    i::Float64
    RAAN::Float64
    w::Float64
    f::Float64
end

# Earth radius [m].
const R0 = 6378136.3;

# Standard gravitational parameter for Earth [m^3/s^2]
const m0 = 3.986004415e14;

# J2 perturbation term.
const J2 = 1.0826269E-03

# Sun radius [m].
const Rs = 6.963e8

include("satellite_lighting_conditions.jl")
include("satellite_position.jl")
include("sun_position.jl")

end # module
