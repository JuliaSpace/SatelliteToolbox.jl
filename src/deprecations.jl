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
@deprecate angvel_to_a(a, e, i, pert; kwargs...) angvel_to_a(a, e, i; pert = pert, kwargs...)

# Deprecations introduced in SatelliteToolbox v0.10
# ==============================================================================

@deprecate compute_RAAN_lt ltan_to_raan
