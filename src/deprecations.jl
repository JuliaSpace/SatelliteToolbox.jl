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

@deprecate satsv orbsv
@deprecate dArgPer dargp
@deprecate dRAAN draan

@deprecate ECEFtoGeodetic ecef_to_geodetic
@deprecate GeodetictoECEF geodetic_to_ecef
@deprecate GeodetictoGeocentric geodetic_to_geocentric

@deprecate rECEFtoECEF r_ecef_to_ecef
@deprecate rECEFtoECI r_ecef_to_eci
@deprecate rECItoECEF r_eci_to_ecef
@deprecate rECItoECI r_eci_to_eci

@deprecate svECEFtoECEF sv_ecef_to_ecef
@deprecate svECEFtoECI sv_ecef_to_eci
@deprecate svECItoECEF sv_eci_to_ecef
@deprecate svECItoECI sv_eci_to_eci

@deprecate rITRFtoPEF_fk5 r_itrf_to_pef_fk5
@deprecate rPEFtoITRF_fk5 r_pef_to_itrf_fk5
@deprecate rPEFtoTOD_fk5 r_pef_to_tod_fk5
@deprecate rTODtoPEF_fk5 r_tod_to_pef_fk5
@deprecate rTODtoMOD_fk5 r_tod_to_mod_fk5
@deprecate rMODtoTOD_fk5 r_mod_to_tod_fk5
@deprecate rMODtoGCRF_fk5 r_mod_to_gcrf_fk5
@deprecate rGCRFtoMOD_fk5 r_gcrf_to_mod_fk5
@deprecate rITRFtoGCRF_fk5 r_itrf_to_gcrf_fk5
@deprecate rGCRFtoITRF_fk5 r_gcrf_to_itrf_fk5
@deprecate rPEFtoMOD_fk5 r_pef_to_mod_fk5
@deprecate rMODtoPEF_fk5 r_mod_to_pef_fk5

@deprecate rITRFtoTIRS_iau2006 r_itrf_to_tirs_iau2006
@deprecate rTIRStoITRF_iau2006 r_tirs_to_itrf_iau2006
@deprecate rTIRStoCIRS_iau2006 r_tirs_to_cirs_iau2006
@deprecate rCIRStoTIRS_iau2006 r_cirs_to_tirs_iau2006
@deprecate rCIRStoGCRF_iau2006 r_cirs_to_gcrf_iau2006
@deprecate rGCRFtoCIRS_iau2006 r_gcrf_to_cirs_iau2006

@deprecate rTIRStoERS_iau2006 r_tirs_to_ers_iau2006
@deprecate rERStoTIRS_iau2006 r_ers_to_tirs_iau2006
@deprecate rERStoMOD_iau2006 r_ers_to_mod_iau2006
@deprecate rMODtoERS_iau2006 r_mod_to_ers_iau2006
@deprecate rMODtoMJ2000_iau2006 r_mod_to_mj2000_iau2006
@deprecate rMJ2000toMOD_iau2006 r_mj2000_to_mod_iau2006
@deprecate rMJ2000toGCRF_iau2006 r_mj2000_to_gcrf_iau2006
@deprecate rGCRFtoMJ2000_iau2006 r_gcrf_to_mj2000_iau2006
@deprecate rTIRStoMOD_iau2006 r_tirs_to_mod_iau2006
@deprecate rMODtoTIRS_iau2006 r_mod_to_tirs_iau2006

@deprecate rTEMEtoTOD r_teme_to_tod
@deprecate rTODtoTEME r_tod_to_teme
@deprecate rTEMEtoMOD r_teme_to_mod
@deprecate rMODtoTEME r_mod_to_teme
@deprecate rTEMEtoGCRF r_teme_to_gcrf
@deprecate rGCRFtoTEME r_gcrf_to_teme
@deprecate rTEMEtoPEF r_teme_to_pef
@deprecate rPEFtoTEME r_pef_to_teme

@deprecate JD_UT1toUTC jd_ut1_to_utc
@deprecate JD_UTCtoUT1 jd_utc_to_ut1
@deprecate JD_UTCtoTT jd_utc_to_tt
@deprecate JD_TTtoUTC jd_tt_to_utc
@deprecate get_ΔAT get_Δat

@deprecate JDtoDate jd_to_date
@deprecate DatetoJD date_to_jd

@deprecate J2000toGMST j2000_to_gmst
@deprecate JDtoGMST jd_to_gmst

@deprecate dEps_dPsi deps_dpsi
