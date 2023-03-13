module SatelliteToolboxSgp4

using Dates, Parameters, StaticArrays
using SatelliteToolbox.SatelliteToolboxTle

################################################################################
#                                  Constants
################################################################################

# Earth Equatorial radius [m].
const R0 = 6378137.0

# Julian Day of J2000.0 epoch.
const JD_J2000 = 2451545.0

################################################################################
#                                    Types
################################################################################

include("types.jl")

################################################################################
#                                   Includes
################################################################################

include("gmst.jl")
include("sgp4_model.jl")
include("helpers.jl")

end
