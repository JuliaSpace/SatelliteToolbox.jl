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
@deprecate svECEFtoECI(args...) sv_ecef_to_eci(args...)
@deprecate svECItoECEF(args...) sv_eci_to_ecef(args...)
@deprecate svECItoECI(args...) sv_eci_to_eci(args...)

@deprecate rITRFtoPEF_fk5(args...) r_itrf_to_pef_fk5(args...)
@deprecate rPEFtoITRF_fk5(args...) r_pef_to_itrf_fk5(args...)
@deprecate rPEFtoTOD_fk5(args...) r_pef_to_tod_fk5(args...)
@deprecate rTODtoPEF_fk5(args...) r_tod_to_pef_fk5(args...)
@deprecate rTODtoMOD_fk5(args...) r_tod_to_mod_fk5(args...)
@deprecate rMODtoTOD_fk5(args...) r_mod_to_tod_fk5(args...)
@deprecate rMODtoGCRF_fk5(args...) r_mod_to_gcrf_fk5(args...)
@deprecate rGCRFtoMOD_fk5(args...) r_gcrf_to_mod_fk5(args...)
@deprecate rITRFtoGCRF_fk5(args...) r_itrf_to_gcrf_fk5(args...)
@deprecate rGCRFtoITRF_fk5(args...) r_gcrf_to_itrf_fk5(args...)
@deprecate rPEFtoMOD_fk5(args...) r_pef_to_mod_fk5(args...)
@deprecate rMODtoPEF_fk5(args...) r_mod_to_pef_fk5(args...)
