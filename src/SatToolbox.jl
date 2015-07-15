module SatToolbox

export OrbitalParameters
export R0, Rm, m0, J2, Rs, ne, au2m, sunRad

import Base: asin, atan2, cos, mod, sin

################################################################################
#                                  Structures
################################################################################

# Structure that defines the orbital parameters of a satellite.
type OrbitalParameters{T}
    a::T
    e::T
    i::T
    RAAN::T
    w::T
    f::T
end

################################################################################
#                                  Constants
################################################################################

# Earth Equatorial radius [m].
const R0 = 6378137.0;

# Earth mean radius [m].
const Rm = 6371009.0;

# Standard gravitational parameter for Earth [m^3/s^2]
const m0 = 3.986004415e14;

# J2 perturbation term.
const J2 = 0.0010826267

# Sun radius [m].
const Rs = 6.963e8

# Earth's orbit mean motion [rad/s]
const ne = (360.0/365.2421897)*pi/180/86400

# Conversion factor from AU to m.
const au2m = 149597870700.0

# Sun radiation emitted [J/sec].
const sunRad = 3.826e26

################################################################################
#                                  Exceptions
################################################################################

# Exception: The perigee is inside the Earth.
type OrbitInvalidPerigee <: Exception
    R_p::Real
end
Base.showerror(io::IO, e::OrbitInvalidPerigee) =
    print(io, "The orbit perigee (", e.R_p, " m) is inside Earth!")

################################################################################
#                                    Files
################################################################################

include("sun_position.jl")

include("orbit.jl")
include("orbit_sun_sync.jl")
include("orbit_sun_sync_ground_reap.jl")

include("payload_optical_analysis.jl")

include("satellite_orbit_step.jl")
include("satellite_lighting_conditions.jl")
include("satellite_position.jl")
include("satellite_beta_angle.jl")
include("satellite_eclipse_time.jl")

include("satellite_sun_angle.jl")
include("satellite_sun_radiation.jl")

end # module
