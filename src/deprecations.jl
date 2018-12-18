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

@deprecate rECEFtoECI(T, M::Type{Val{:FK5}}, ECEF, ECI, JD_UTC) rECEFtoECI(T, ECEF, ECI, JD_UTC)
@deprecate rECEFtoECI(M::Type{Val{:FK5}}, ECEF, ECI, JD_UTC)    rECEFtoECI(ECEF, ECI, JD_UTC)

@deprecate rECEFtoECI(T, M::Type{Val{:FK5}}, ECEF, ECI, JD_UTC, eop_data) rECEFtoECI(T, ECEF, ECI, JD_UTC, eop_data)
@deprecate rECEFtoECI(M::Type{Val{:FK5}}, ECEF, ECI, JD_UTC, eop_data)    rECEFtoECI(ECEF, ECI, JD_UTC, eop_data)

@deprecate rECItoECEF(T, M::Type{Val{:FK5}}, ECI, ECEF, JD_UTC) rECItoECEF(T, ECI, ECEF, JD_UTC)
@deprecate rECItoECEF(M::Type{Val{:FK5}}, ECI, ECEF, JD_UTC)    rECItoECEF(ECI, ECEF, JD_UTC)

@deprecate rECItoECEF(T, M::Type{Val{:FK5}}, ECI, ECEF, JD_UTC, eop_data) rECItoECEF(T, ECI, ECEF, JD_UTC, eop_data)
@deprecate rECItoECEF(M::Type{Val{:FK5}}, ECI, ECEF, JD_UTC, eop_data)    rECItoECEF(ECI, ECEF, JD_UTC, eop_data)

@deprecate rECEFtoECEF(T, M::Type{Val{:FK5}}, ECEFo, ECEFf, JD_UTC, eop_data) rECEFtoECEF(T, ECEFo, ECEFf, JD_UTC, eop_data)
@deprecate rECEFtoECEF(M::Type{Val{:FK5}}, ECEFo, ECEFf, JD_UTC, eop_data)    rECEFtoECEF(ECEFo, ECEFf, JD_UTC, eop_data)

@deprecate rECItoECI(T, M::Type{Val{:FK5}}, ECIo, ECIf, JD_UTC, eop_data) rECItoECI(T, ECIo, ECIf, JD_UTC, eop_data)
@deprecate rECItoECI(M::Type{Val{:FK5}}, ECIo, ECIf, JD_UTC, eop_data)    rECItoECI(ECIo, ECIf, JD_UTC, eop_data)

@deprecate rECItoECI(T, M::Type{Val{:FK5}}, ECIo, JD_UTCo, ECIf, JD_UTCf, eop_data) rECItoECI(T, ECIo, JD_UTCo, ECIf, JD_UTCf, eop_data)
@deprecate rECItoECI(M::Type{Val{:FK5}}, ECIo, JD_UTCo, ECIf, JD_UTCf, eop_data)    rECItoECI(ECIo, JD_UTCo, ECIf, JD_UTCf, eop_data)


