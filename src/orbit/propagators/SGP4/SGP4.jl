module SGP4

using Parameters, StaticArrays, Printf, Dates

# Earth Equatorial radius [m].
const R0 = 6378137.0

# Julian Day of J2000.0 epoch.
const JD_J2000 = 2451545.0

include("sgp4_types.jl")
include("../../../transformations/gmst.jl")
include("tle.jl")
include("sgp4_model.jl")

function sgp4_init(tle::TLE,
                   sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84)

    sgp4_init(sgp4_gc,
              tle.epoch,
              tle.n*2*pi/(24*60),
              tle.e,
              tle.i*pi/180,
              tle.Ω*pi/180,
              tle.ω*pi/180,
              tle.M*pi/180,
              tle.bstar)
end

end
