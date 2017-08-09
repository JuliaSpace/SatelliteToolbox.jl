VERSION >= v"0.4.0-dev+6521" && __precompile__()

module SatToolbox

export JD_J2000, R0, Rm, m0, J2, Rs, ne, au2m, sunRad
export a_wgs84, b_wgs84, f_wgs84, e_wgs84, el_wgs84

import Base: asin, atan2, cos, getindex, mod, sin

################################################################################
#                                  Structures
################################################################################

################################################################################
#                                  Constants
################################################################################

include("constants.jl")

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

include("coordinate_transformations/coordinate_transformations.jl")
include("coordinate_transformations/orbit_elements.jl")
include("coordinate_transformations/position.jl")

include("sun/equation_of_time.jl")
include("sun/sun_position.jl")

include("orbit/orbit.jl")
include("orbit/orbit_aux.jl")
include("orbit/orbit_sun_sync.jl")
include("orbit/orbit_sun_sync_ground_reap.jl")
include("orbit/satellite_velocity.jl")
include("orbit/sgp4.jl")

include("analysis/beta_angle.jl")
include("analysis/eclipse_time.jl")
include("analysis/lighting_conditions.jl")
include("analysis/raan.jl")
include("analysis/payload_optical_analysis.jl")
include("analysis/satellite_position_countries.jl")
include("analysis/satellite_position_groundstations.jl")
include("analysis/sun_angle.jl")
include("analysis/sun_radiation.jl")

end # module
