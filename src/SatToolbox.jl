VERSION >= v"0.4.0-dev+6521" && __precompile__()

module SatToolbox

export JD_J2000, R0, Rm, m0, J2, Rs, ne, au2m, sunRad
export a_wgs84, b_wgs84, f_wgs84, e_wgs84, el_wgs84

import Base: asin, atan2, cos, mod, sin

################################################################################
#                                  Structures
################################################################################

################################################################################
#                                  Constants
################################################################################

# Julian Day of J2000.0 epoch.
const JD_J2000 = 2451545.0

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

# WGS-84 Data.
const a_wgs84  = 6378137.0
const f_wgs84  = 1/298.257223563
const b_wgs84  = a_wgs84*(1-f_wgs84)
const e_wgs84  = sqrt( (a_wgs84^2-b_wgs84^2)/a_wgs84^2 )
const el_wgs84 = sqrt( (a_wgs84^2-b_wgs84^2)/b_wgs84^2 )

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

include("coordinate_transformations.jl")

include("sun_position.jl")

include("orbit.jl")
include("orbit_sun_sync.jl")
include("orbit_sun_sync_ground_reap.jl")

include("payload_optical_analysis.jl")

include("satellite_orbit_step.jl")
include("satellite_lighting_conditions.jl")
include("satellite_position.jl")
include("satellite_position_countries.jl")
include("satellite_position_groundstations.jl")
include("satellite_beta_angle.jl")
include("satellite_eclipse_time.jl")

include("satellite_sun_angle.jl")
include("satellite_sun_radiation.jl")

end # module
