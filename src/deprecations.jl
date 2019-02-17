#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Deprecation warnings.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

#                     Introduced in SatelliteToolbox v0.4
# ==============================================================================

@deprecate parse_gfc(filename) create_gravity_model_coefs(parse_icgem(filename))
@deprecate satellite_orbit_compute_f(a, e, i, M) M_to_f(e, M)
@deprecate satellite_orbit_compute_f(a, e, i, M, tol) M_to_f(e, M, tol)
