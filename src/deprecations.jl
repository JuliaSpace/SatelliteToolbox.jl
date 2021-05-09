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

@deprecate rITRFtoTIRS_iau2006(args...) r_itrf_to_tirs_iau2006(args...)
@deprecate rTIRStoITRF_iau2006(args...) r_tirs_to_itrf_iau2006(args...)
@deprecate rTIRStoCIRS_iau2006(args...) r_tirs_to_cirs_iau2006(args...)
@deprecate rCIRStoTIRS_iau2006(args...) r_cirs_to_tirs_iau2006(args...)
@deprecate rCIRStoGCRF_iau2006(args...) r_cirs_to_gcrf_iau2006(args...)
@deprecate rGCRFtoCIRS_iau2006(args...) r_gcrf_to_cirs_iau2006(args...)

@deprecate rTIRStoERS_iau2006(args...) r_tirs_to_ers_iau2006(args...)
@deprecate rERStoTIRS_iau2006(args...) r_ers_to_tirs_iau2006(args...)
@deprecate rERStoMOD_iau2006(args...) r_ers_to_mod_iau2006(args...)
@deprecate rMODtoERS_iau2006(args...) r_mod_to_ers_iau2006(args...)
@deprecate rMODtoMJ2000_iau2006(args...) r_mod_to_mj2000_iau2006(args...)
@deprecate rMJ2000toMOD_iau2006(args...) r_mj2000_to_mod_iau2006(args...)
@deprecate rMJ2000toGCRF_iau2006(args...) r_mj2000_to_gcrf_iau2006(args...)
@deprecate rGCRFtoMJ2000_iau2006(args...) r_gcrf_to_mj2000_iau2006(args...)
@deprecate rTIRStoMOD_iau2006(args...) r_tirs_to_mod_iau2006(args...)
@deprecate rMODtoTIRS_iau2006(args...) r_mod_to_tirs_iau2006(args...)

@deprecate JD_UT1toUTC(args...) jd_ut1_to_utc(args...)
@deprecate JD_UTCtoUT1(args...) jd_utc_to_ut1(args...)
@deprecate JD_UTCtoTT(args...) jd_utc_to_tt(args...)
@deprecate JD_TTtoUTC(args...) jd_tt_to_utc(args...)
@deprecate get_ΔAT(args...) get_Δat(args...)

@deprecate JDtoDate jd_to_date
@deprecate DatetoJD date_to_jd

@deprecate J2000toGMST j2000_to_gmst
@deprecate JDtoGMST jd_to_gmst

