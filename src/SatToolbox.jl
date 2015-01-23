module SatToolbox

export SAT_LIGHTING_SUNLIGHT, SAT_LIGHTING_UMBRA, SAT_LIGHTING_PENUMBRA
export compute_ss_orbit_by_ang_vel, compute_ss_orbit_by_num_rev_per_day
export compute_ss_orbit_by_semi_major_axis
export compute_ss_orbit_by_inclination
export list_ss_orbits_by_rep_period, sort_list_ss_orbits_by_height
export dRAAN_J2, dw_J2, n_J0, n_J2, t_J0, t_J2
export satellite_lighting_condition, satellite_position_i, sun_position_i
export satellite_orbit_compute_f
export satellite_beta_angle
export satellite_eclipse_time
export compute_RAAN_lt
export minimum_swath_grss, minimum_swath_orbit_grss
export OrbitalParameters
export R0, m0, J2, Rs

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

# Earth radius [m].
const R0 = 6378136.3;

# Standard gravitational parameter for Earth [m^3/s^2]
const m0 = 3.986004415e14;

# J2 perturbation term.
const J2 = 1.0826269E-03

# Sun radius [m].
const Rs = 6.963e8

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

include("satellite_orbit_step.jl")
include("satellite_lighting_conditions.jl")
include("satellite_position.jl")
include("satellite_beta_angle.jl")
include("satellite_eclipse_time.jl")

end # module
