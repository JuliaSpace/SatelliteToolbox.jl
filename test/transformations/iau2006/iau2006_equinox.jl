# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to equinox-based IAU-2006 transformations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/transformations/iau2006/iau2006_equinox.jl
# ======================================================

# Functions rTIRStoERS_iau2006 and rERStoTIRS_iau2006
# ---------------------------------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-14: Performing an IAU-2000 reduction [1, p. 220]
#
# According to this example and Table 3-6, using:
#
#   JD_UT1 = 2453101.827406783
#   JD_TT  = 2453101.828154745
#   r_tirs  = -1033.47503120   i + 7901.30558560   j + 6380.34453270   k [km]
#   v_tirs  =    -3.2256327470 i -    2.8724425110 j +    5.5319312880 k [km/s]
#
# one gets the following:
#
#   r_ers  = +5094.51462800   i + 6127.36658790   j + 6380.34453270   k [km]
#   v_ers  =    -4.7460885870 i -    0.7860771040 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Functions rTIRStoERS_iau2006 and rERStoTIRS_iau2006" begin
    JD_UT1 = 2453101.827406783
    JD_TT  = 2453101.828154745
    LOD    = 0.0015563
    w      = 7.292115146706979e-5*(1-LOD/86400)

    # rTIRStoERS_iau2006
    # ==================

    r_tirs  = [-1033.47503120; 7901.30558560; 6380.34453270]
    v_tirs  = [-3.2256327470; -2.8724425110; +5.5319312880]

    # DCM
    # ---

    D_ERS_TIRS = rTIRStoERS_iau2006(JD_UT1, JD_TT)

    r_ers = D_ERS_TIRS*r_tirs
    v_ers = D_ERS_TIRS*(v_tirs + [0;0;w] × r_tirs)

    @test r_ers[1] ≈ +5094.51462800 atol=5e-6
    @test r_ers[2] ≈ +6127.36658790 atol=5e-6
    @test r_ers[3] ≈ +6380.34453270 atol=5e-6

    @test v_ers[1] ≈ -4.7460885870  atol=1e-9
    @test v_ers[2] ≈ +0.7860771040  atol=1e-9
    @test v_ers[3] ≈ +5.5319312880  atol=1e-9

    # Quaternion
    # ----------

    q_ERS_TIRS = rTIRStoERS_iau2006(Quaternion, JD_UT1, JD_TT)

    r_ers = vect(q_ERS_TIRS\r_tirs*q_ERS_TIRS)
    v_ers = vect(q_ERS_TIRS\(v_tirs + [0;0;w] × r_tirs)*q_ERS_TIRS)

    @test r_ers[1] ≈ +5094.51462800 atol=5e-6
    @test r_ers[2] ≈ +6127.36658790 atol=5e-6
    @test r_ers[3] ≈ +6380.34453270 atol=5e-6

    @test v_ers[1] ≈ -4.7460885870  atol=1e-9
    @test v_ers[2] ≈ +0.7860771040  atol=1e-9
    @test v_ers[3] ≈ +5.5319312880  atol=1e-9

    # rERStoTIRS_iau2006
    # ==================

    r_ers = [+5094.51462800; +6127.36658790; +6380.34453270]
    v_ers = [-4.7460885870; +0.7860771040; +5.5319312880]

    # DCM
    # ---

    D_TIRS_ERS = rERStoTIRS_iau2006(JD_UT1, JD_TT)

    r_tirs = D_TIRS_ERS*r_ers
    v_tirs = D_TIRS_ERS*v_ers - [0;0;w] × r_tirs

    @test r_tirs[1] ≈ -1033.47503120 atol=5e-6
    @test r_tirs[2] ≈ +7901.30558560 atol=5e-6
    @test r_tirs[3] ≈ +6380.34453270 atol=5e-6

    @test v_tirs[1] ≈ -3.2256327470  atol=1e-9
    @test v_tirs[2] ≈ -2.8724425110  atol=1e-9
    @test v_tirs[3] ≈ +5.5319312880  atol=1e-9

    # Quaternion
    # ----------

    q_TIRS_ERS = rERStoTIRS_iau2006(Quaternion, JD_UT1, JD_TT)

    r_tirs = vect(q_TIRS_ERS\r_ers*q_TIRS_ERS)
    v_tirs = vect(q_TIRS_ERS\v_ers*q_TIRS_ERS) - [0;0;w] × r_tirs

    @test r_tirs[1] ≈ -1033.47503120 atol=5e-6
    @test r_tirs[2] ≈ +7901.30558560 atol=5e-6
    @test r_tirs[3] ≈ +6380.34453270 atol=5e-6

    @test v_tirs[1] ≈ -3.2256327470  atol=1e-9
    @test v_tirs[2] ≈ -2.8724425110  atol=1e-9
    @test v_tirs[3] ≈ +5.5319312880  atol=1e-9
end

# Functions rERStoMOD_iau2006 and rMODtoERS_iau2006
# --------------------------------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-14: Performing an IAU-2000 reduction [1, p. 220]
#
# According to this example and Table 3-6, using:
#
#   JD_TT  = 2453101.828154745
#   r_ers  = +5094.51462800   i + 6127.36658790   j + 6380.34453270   k [km]
#   v_ers  =    -4.7460885870 i -    0.7860771040 j +    5.5319312880 k [km/s]
#
# one gets the following:
#
#   r_mod  = +5094.02896110   i + 6127.87113500   j + 6380.24774200   k [km]
#   v_mod  =    -4.7462624800 i +    0.7860141930 j +    5.5317910320 k [km/s]
#
################################################################################

@testset "Functions rERStoMOD_iau2006 and rMODtoERS_iau2006" begin
    JD_TT  = 2453101.828154745

    # rERStoMOD_iau2006
    # =================

    r_ers = [+5094.51462800; +6127.36658790; +6380.34453270]
    v_ers = [-4.7460885870; +0.7860771040; +5.5319312880]

    # DCM
    # ---

    D_MOD_ERS = rERStoMOD_iau2006(JD_TT)

    r_mod = D_MOD_ERS*r_ers
    v_mod = D_MOD_ERS*v_ers

    @test r_mod[1] ≈ +5094.02896110 atol=1e-7
    @test r_mod[2] ≈ +6127.87113500 atol=1e-7
    @test r_mod[3] ≈ +6380.24774200 atol=1e-7

    @test v_mod[1] ≈ -4.7462624800  atol=1e-9
    @test v_mod[2] ≈ +0.7860141930  atol=1e-9
    @test v_mod[3] ≈ +5.5317910320  atol=1e-9

    # Quaternion
    # ----------

    q_MOD_ERS = rERStoMOD_iau2006(Quaternion, JD_TT)

    r_mod = vect(q_MOD_ERS\r_ers*q_MOD_ERS)
    v_mod = vect(q_MOD_ERS\v_ers*q_MOD_ERS)

    @test r_mod[1] ≈ +5094.02896110 atol=1e-7
    @test r_mod[2] ≈ +6127.87113500 atol=1e-7
    @test r_mod[3] ≈ +6380.24774200 atol=1e-7

    @test v_mod[1] ≈ -4.7462624800  atol=1e-9
    @test v_mod[2] ≈ +0.7860141930  atol=1e-9
    @test v_mod[3] ≈ +5.5317910320  atol=1e-9

    # rMODtoERS_iau2006
    # =================

    r_mod = [+5094.02896110; +6127.87113500; +6380.24774200]
    v_mod = [-4.7462624800; +0.7860141930; +5.5317910320]

    # DCM
    # ---

    D_ERS_MOD = rMODtoERS_iau2006(JD_TT)

    r_ers = D_ERS_MOD*r_mod
    v_ers = D_ERS_MOD*v_mod

    @test r_ers[1] ≈ +5094.51462800 atol=1e-7
    @test r_ers[2] ≈ +6127.36658790 atol=1e-7
    @test r_ers[3] ≈ +6380.34453270 atol=1e-7

    @test v_ers[1] ≈ -4.7460885870  atol=1e-9
    @test v_ers[2] ≈ +0.7860771040  atol=1e-9
    @test v_ers[3] ≈ +5.5319312880  atol=1e-9

    # Quaternion
    # ----------

    q_ERS_MOD = rMODtoERS_iau2006(Quaternion, JD_TT)

    r_ers = vect(q_ERS_MOD\r_mod*q_ERS_MOD)
    v_ers = vect(q_ERS_MOD\v_mod*q_ERS_MOD)

    @test r_ers[1] ≈ +5094.51462800 atol=1e-7
    @test r_ers[2] ≈ +6127.36658790 atol=1e-7
    @test r_ers[3] ≈ +6380.34453270 atol=1e-7

    @test v_ers[1] ≈ -4.7460885870  atol=1e-9
    @test v_ers[2] ≈ +0.7860771040  atol=1e-9
    @test v_ers[3] ≈ +5.5319312880  atol=1e-9
end

# Functions rMODtoMJ2000_iau2006, rMJ2000toMOD_iau2006, rMJ2000toGCRF_iau2006, and rGCRFtoMJ2000_iau2006
# ------------------------------------------------------------------------------------------------------

# We are testing all these functions together because we have in [1] the result
# MOD <=> GCRF using equinox-based IAU-2006 theory, but we do not have the
# intermediate result MOD <=> MJ2000.

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-14: Performing an IAU-2000 reduction [1, p. 220]
#
# According to this example and Table 3-6, using:
#
#   JD_TT  = 2453101.828154745
#   r_mod  = +5094.02896110   i + 6127.87113500   j + 6380.24774200   k [km]
#   v_mod  =    -4.7462624800 i +    0.7860141930 j +    5.5317910320 k [km/s]
#
# one gets the following:
#
#   r_gcrf = 5102.50895780  i + 6123.01140380 j + 6378.13692530 k [km]
#   v_gcrf =  -4.7432201560 i + 0.7905364970  j + 5.5337557280  k [km/s]
#
################################################################################

@testset "Functions rMODtoMJ2000_iau2006, rMJ2000toMOD_iau2006, rMJ2000toGCRF_iau2006, and rGCRFtoMJ2000_iau2006" begin
    JD_TT  = 2453101.828154745

    # MOD to GCRF
    # ===========

    r_mod = [+5094.02896110; +6127.87113500; +6380.24774200]
    v_mod = [-4.7462624800; +0.7860141930; +5.5317910320]

    # DCM
    # ---

    D_GCRF_MOD = rMJ2000toGCRF_iau2006(JD_TT)*rMODtoMJ2000_iau2006(JD_TT)

    r_gcrf = D_GCRF_MOD*r_mod
    v_gcrf = D_GCRF_MOD*v_mod

    @test r_gcrf[1] ≈ +5102.50895780 atol=1e-7
    @test r_gcrf[2] ≈ +6123.01140380 atol=1e-7
    @test r_gcrf[3] ≈ +6378.13692530 atol=1e-7

    @test v_gcrf[1] ≈ -4.7432201560  atol=1e-9
    @test v_gcrf[2] ≈ +0.7905364970  atol=1e-9
    @test v_gcrf[3] ≈ +5.5337557280  atol=1e-9

    # Quaternion
    # ----------

    q_GCRF_MOD = rMODtoMJ2000_iau2006(Quaternion, JD_TT) *
                 rMJ2000toGCRF_iau2006(Quaternion, JD_TT)

    r_gcrf = vect(q_GCRF_MOD\r_mod*q_GCRF_MOD)
    v_gcrf = vect(q_GCRF_MOD\v_mod*q_GCRF_MOD)

    @test r_gcrf[1] ≈ +5102.50895780 atol=1e-7
    @test r_gcrf[2] ≈ +6123.01140380 atol=1e-7
    @test r_gcrf[3] ≈ +6378.13692530 atol=1e-7

    @test v_gcrf[1] ≈ -4.7432201560  atol=1e-9
    @test v_gcrf[2] ≈ +0.7905364970  atol=1e-9
    @test v_gcrf[3] ≈ +5.5337557280  atol=1e-9

    # GCRF to MOD
    # ===========

    r_gcrf = [+5102.50895780; +6123.01140380; +6378.13692530]
    v_gcrf = [-4.7432201560; +0.7905364970; +5.5337557280]

    # DCM
    # ---

    D_MOD_GCRF = rMJ2000toMOD_iau2006(JD_TT) * rGCRFtoMJ2000_iau2006(JD_TT)

    r_mod = D_MOD_GCRF*r_gcrf
    v_mod = D_MOD_GCRF*v_gcrf

    @test r_mod[1] ≈ +5094.02896110 atol=1e-7
    @test r_mod[2] ≈ +6127.87113500 atol=1e-7
    @test r_mod[3] ≈ +6380.24774200 atol=1e-7

    @test v_mod[1] ≈ -4.7462624800  atol=1e-9
    @test v_mod[2] ≈ +0.7860141930  atol=1e-9
    @test v_mod[3] ≈ +5.5317910320  atol=1e-9

    # Quaternion
    # ----------

    q_MOD_GCRF = rGCRFtoMJ2000_iau2006(Quaternion, JD_TT) *
                 rMJ2000toMOD_iau2006(Quaternion, JD_TT)

    r_mod = vect(q_MOD_GCRF\r_gcrf*q_MOD_GCRF)
    v_mod = vect(q_MOD_GCRF\v_gcrf*q_MOD_GCRF)

    @test r_mod[1] ≈ +5094.02896110 atol=1e-7
    @test r_mod[2] ≈ +6127.87113500 atol=1e-7
    @test r_mod[3] ≈ +6380.24774200 atol=1e-7

    @test v_mod[1] ≈ -4.7462624800  atol=1e-9
    @test v_mod[2] ≈ +0.7860141930  atol=1e-9
    @test v_mod[3] ≈ +5.5317910320  atol=1e-9
end
