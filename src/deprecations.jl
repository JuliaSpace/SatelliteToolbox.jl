# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Deprecation warnings.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Deprecations introduced in SatelliteToolbox v0.11
# ==============================================================================

@deprecate angvel(a, e, i, pert) angvel(a, e, i; pert = pert)

# Deprecations introduced in SatelliteToolbox v0.10
# ==============================================================================

@deprecate compute_RAAN_lt ltan_to_raan
