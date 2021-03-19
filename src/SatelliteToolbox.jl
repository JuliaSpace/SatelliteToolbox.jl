module SatelliteToolbox

export JD_J2000, R0, Rm, m0, J2, Rs, ne, au2m, sunRad
export a_wgs84, b_wgs84, f_wgs84, e_wgs84, el_wgs84

import Base: asin, atan, convert, copy, cos, deepcopy, getindex, length,
             iterate, mod, setindex!, sin, show

using Crayons
using Dates
using DelimitedFiles
using Interpolations
using LinearAlgebra
using OptionalData
using Parameters
using PolynomialRoots
using PrettyTables
using Printf
using ReferenceFrameRotations
using RemoteFiles
using StaticArrays
using SparseArrays
using Statistics
using Reexport

# Re-exporting symbols from ReferenceFrameRotations.jl.
export DCM
export Quaternion

################################################################################
#                                  Submodules
################################################################################

include("submodules/SatelliteToolboxTLE/SatelliteToolboxTLE.jl")
@reexport using .SatelliteToolboxTLE

include("submodules/SGP4/SGP4.jl")
@reexport using .SGP4

################################################################################
#                             Types and Structures
################################################################################

include("./types/types.jl")

################################################################################
#                                  Constants
################################################################################

include("constants.jl")

# Pre-defined crayons.
const _reset_crayon = Crayon(reset = true)
const _crayon_bold  = crayon"bold"
const _crayon_g     = crayon"bold green"
const _crayon_u     = crayon"bold blue"
const _crayon_y     = crayon"bold yellow"

# Escape sequences related to the crayons.
const _b = _crayon_bold
const _d = _reset_crayon
const _g = _crayon_g
const _y = _crayon_y
const _u = _crayon_u

################################################################################
#                                  Exceptions
################################################################################

################################################################################
#                                    Files
################################################################################

include("analysis/beta_angle.jl")
include("analysis/eclipse_time.jl")
include("analysis/lighting_conditions.jl")
include("analysis/ground_stations.jl")
include("analysis/ground_trace.jl")
include("analysis/raan.jl")
include("analysis/payload_optical_analysis.jl")
include("analysis/satellite_position_countries.jl")
include("analysis/sun_angle.jl")
include("analysis/sun_radiation.jl")

include("earth/atmospheric_models/expatmosphere/expatmosphere.jl")
include("earth/atmospheric_models/jb2008/jb2008.jl")
include("earth/atmospheric_models/jr1971/jr1971.jl")
include("earth/atmospheric_models/nrlmsise00/nrlmsise00.jl")
include("earth/gravity_models/embedded_gravity_models.jl")
include("earth/gravity_models/gravity_model.jl")
include("earth/geomagnetic_field_models/dipole/geomag_dipole.jl")
include("earth/geomagnetic_field_models/dipole/geomag_dipole_constants.jl")
include("earth/geomagnetic_field_models/igrf/igrf.jl")
include("earth/geomagnetic_field_models/igrf/igrf_coefs.jl")
include("earth/geomagnetic_field_models/igrf/igrf12/igrf12syn.jl")
include("earth/geomagnetic_field_models/igrf/igrf12/igrf12syn_coefs.jl")
include("earth/geomagnetic_field_models/igrf/igrf13/igrf13syn.jl")
include("earth/space_indices/space_indices.jl")

include("./misc/crossing.jl")
include("./misc/legendre.jl")
include("./misc/dlegendre.jl")
include("./misc/icgem.jl")

include("./moon/moon_position.jl")

include("sun/equation_of_time.jl")
include("sun/sun_position.jl")
include("sun/sun_velocity.jl")

include("orbit/general.jl")
include("orbit/anomalies.jl")
include("orbit/mean_elements.jl")
include("orbit/orbit_sun_sync.jl")
include("orbit/orbit_sun_sync_ground_reap.jl")

include("orbit/propagators/j2.jl")
include("orbit/propagators/j2osc.jl")
include("orbit/propagators/j4.jl")
include("orbit/propagators/twobody.jl")
include("orbit/propagators/api/api.jl")
include("orbit/propagators/api/j2.jl")
include("orbit/propagators/api/j4.jl")
include("orbit/propagators/api/sgp4.jl")
include("orbit/propagators/api/twobody.jl")

include("orbit/representations/api.jl")
include("orbit/representations/conversions.jl")
include("orbit/representations/keplerian_elements.jl")
include("orbit/representations/rv.jl")
include("orbit/representations/state_vector.jl")

include("transformations/eop.jl")
include("transformations/ecef_to_ecef.jl")
include("transformations/ecef_to_eci.jl")
include("transformations/eci_to_ecef.jl")
include("transformations/eci_to_eci.jl")
include("transformations/geodetic_geocentric.jl")
include("transformations/gmst.jl")
include("transformations/misc.jl")
include("transformations/orbit_elements.jl")
include("transformations/position.jl")
include("transformations/sv_ecef_to_ecef.jl")
include("transformations/sv_ecef_to_eci.jl")
include("transformations/sv_eci_to_ecef.jl")
include("transformations/sv_eci_to_eci.jl")

include("transformations/fk5/fk5.jl")
include("transformations/fk5/nutation.jl")
include("transformations/fk5/precession.jl")

include("transformations/iau2006/cio.jl")
include("transformations/iau2006/iau2006.jl")
include("transformations/iau2006/fundamental_args.jl")
include("transformations/iau2006/misc.jl")

include("transformations/iau2006/constants/cio_s.jl")
include("transformations/iau2006/constants/cip_x.jl")
include("transformations/iau2006/constants/cip_y.jl")

include("transformations/teme/teme.jl")

include("time/julian_day.jl")
include("time/time.jl")

include("deprecations.jl")

end # module
