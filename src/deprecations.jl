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
