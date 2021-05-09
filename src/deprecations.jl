# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Deprecation warnings.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Deprecations introduced in SatelliteToolbox v0.9
# ==============================================================================

@deprecate satsv(args...) orbsv(args...)
@deprecate dArgPer(args...) dargp(args...)
@deprecate dRAAN(args...) draan(args...)

@deprecate ECEFtoGeodetic(args...) ecef_to_geodetic(args...)
@deprecate GeodetictoECEF(args...) geodetic_to_ecef(args...)
@deprecate GeodetictoGeocentric(args...) geodetic_to_geocentric(args...)

@deprecate rECEFtoECEF(args...) r_ecef_to_ecef(args...)
@deprecate rECEFtoECI(args...) r_ecef_to_eci(args...)
@deprecate rECItoECEF(args...) r_eci_to_ecef(args...)
@deprecate rECItoECI(args...) r_eci_to_eci(args...)

@deprecate svECEFtoECEF(args...) sv_ecef_to_ecef(args...)
