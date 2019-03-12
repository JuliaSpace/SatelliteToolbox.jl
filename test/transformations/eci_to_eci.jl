#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to ECI to ECI transformations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# Get the current EOP Data.
#
# TODO: The EOP data obtained from IERS website does not match the values in the
# examples in [1]. However, this should be enough, because 1) the individual
# functions at the low level are tested using the same values of [1], and 2) the
# difference is smaller than 30 cm.

eop_iau1980  = read_iers_eop("./eop_IAU1980.txt",  :IAU1980)
eop_iau2000a = read_iers_eop("./eop_IAU2000A.txt", :IAU2000A)

# File: ./src/transformations/eci_to_eci.jl
# =========================================

################################################################################
#                                  IAU-76/FK5
################################################################################

# Functions: rECItoECI
# --------------------

# GCRF <=> J2000
# ==============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_gcrf = 5102.50895790   i + 6123.01140070   j + 6378.13692820   k [km]
#   v_gcrf =   -4.7432201570 i +    0.7905364970 j +    5.5337557270 k [km/s]
#
# one gets:
#
#   r_j2000 = 5102.50960000   i + 6123.01152000   j + 6378.13630000   k [km]
#   v_j2000 =   -4.7432196000 i +    0.7905366000 j +    5.5337561900 k [km/s]
#
################################################################################

@testset "Function rECItoECI GCRF <=> J2000" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## GCRF => J2000
    ## =============

    r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]
    v_gcrf = [-4.7432201570; 0.7905364970; 5.5337557270]

    ## DCM
    ## ---

    D_J2000_GCRF = rECItoECI(GCRF(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = D_J2000_GCRF*r_gcrf
    v_j2000 = D_J2000_GCRF*v_gcrf

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    @test v_j2000[1] ≈ -4.7432196000  atol=1e-7
    @test v_j2000[2] ≈ +0.7905366000  atol=1e-7
    @test v_j2000[3] ≈ +5.5337561900  atol=1e-7

    ## Quaternion
    ## ----------

    q_J2000_GCRF = rECItoECI(Quaternion, GCRF(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = vect(conj(q_J2000_GCRF)*r_gcrf*q_J2000_GCRF)
    v_j2000 = vect(conj(q_J2000_GCRF)*v_gcrf*q_J2000_GCRF)

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    @test v_j2000[1] ≈ -4.7432196000  atol=1e-7
    @test v_j2000[2] ≈ +0.7905366000  atol=1e-7
    @test v_j2000[3] ≈ +5.5337561900  atol=1e-7

    ## J2000 => GCRF
    ## =============

    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]
    v_j2000 = [-4.7432196000; 0.7905366000; 5.5337561900]

    ## DCM
    ## ---

    D_GCRF_J2000 = rECItoECI(J2000(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = D_GCRF_J2000*r_j2000
    v_gcrf = D_GCRF_J2000*v_j2000

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4

    @test v_gcrf[1] ≈ -4.7432201570  atol=1e-7
    @test v_gcrf[2] ≈ +0.7905364970  atol=1e-7
    @test v_gcrf[3] ≈ +5.5337557270  atol=1e-7

    ## Quaternion
    ## ----------

    q_GCRF_J2000 = rECItoECI(Quaternion, J2000(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = vect(conj(q_GCRF_J2000)*r_j2000*q_GCRF_J2000)
    v_gcrf = vect(conj(q_GCRF_J2000)*v_j2000*q_GCRF_J2000)

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4

    @test v_gcrf[1] ≈ -4.7432201570  atol=1e-7
    @test v_gcrf[2] ≈ +0.7905364970  atol=1e-7
    @test v_gcrf[3] ≈ +5.5337557270  atol=1e-7
end

# GCRF <=> MOD
# ============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_gcrf = 5102.50895790   i + 6123.01140070   j + 6378.13692820   k [km]
#   v_gcrf =   -4.7432201570 i +    0.7905364970 j +    5.5337557270 k [km/s]
#
# one gets:
#
#   r_mod = 5094.02837450   i + 6127.87081640   j + 6380.24851640   k [km]
#   v_mod =   -4.7462630520 i +    0.7860140450 j +    5.5317905620 k [km/s]
#
################################################################################

@testset "Function rECItoECI GCRF <=> MOD" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## GCRF => MOD
    ## ===========

    r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]
    v_gcrf = [-4.7432201570; 0.7905364970; 5.5337557270]

    ## DCM
    ## ---

    D_MOD_GCRF = rECItoECI(GCRF(), MOD(), JD_UTC, eop_iau1980)

    r_mod = D_MOD_GCRF*r_gcrf
    v_mod = D_MOD_GCRF*v_gcrf

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    @test v_mod[1] ≈ -4.7462630520  atol=1e-7
    @test v_mod[2] ≈ +0.7860140450  atol=1e-7
    @test v_mod[3] ≈ +5.5317905620  atol=1e-7

    ## Quaternion
    ## ----------

    q_MOD_GCRF = rECItoECI(Quaternion, GCRF(), MOD(), JD_UTC, eop_iau1980)

    r_mod = vect(conj(q_MOD_GCRF)*r_gcrf*q_MOD_GCRF)
    v_mod = vect(conj(q_MOD_GCRF)*v_gcrf*q_MOD_GCRF)

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    @test v_mod[1] ≈ -4.7462630520  atol=1e-7
    @test v_mod[2] ≈ +0.7860140450  atol=1e-7
    @test v_mod[3] ≈ +5.5317905620  atol=1e-7

    ## MOD => GCRF
    ## ===========

    r_mod = [5094.02837450; 6127.87081640; 6380.24851640]
    v_mod = [-4.7462630520; 0.7860140450; 5.5317905620]

    ## DCM
    ## ---

    D_GCRF_MOD = rECItoECI(MOD(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = D_GCRF_MOD*r_mod
    v_gcrf = D_GCRF_MOD*v_mod

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4

    @test v_gcrf[1] ≈ -4.7432201570  atol=1e-7
    @test v_gcrf[2] ≈ +0.7905364970  atol=1e-7
    @test v_gcrf[3] ≈ +5.5337557270  atol=1e-7

    ## Quaternion
    ## ----------

    q_GCRF_MOD = rECItoECI(Quaternion, MOD(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = vect(conj(q_GCRF_MOD)*r_mod*q_GCRF_MOD)
    v_gcrf = vect(conj(q_GCRF_MOD)*v_mod*q_GCRF_MOD)

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4

    @test v_gcrf[1] ≈ -4.7432201570  atol=1e-7
    @test v_gcrf[2] ≈ +0.7905364970  atol=1e-7
    @test v_gcrf[3] ≈ +5.5337557270  atol=1e-7
end

# GCRF <=> TOD
# ============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_gcrf = 5102.50895790   i + 6123.01140070   j + 6378.13692820   k [km]
#   v_gcrf =   -4.7432201570 i +    0.7905364970 j +    5.5337557270 k [km/s]
#
# one gets:
#
#   r_tod = 5094.51620300   i + 6127.36527840   j + 6380.34453270   k [km]
#   v_tod =   -4.7460883850 i +    0.7860783240 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Function rECItoECI GCRF <=> TOD" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## GCRF => TOD
    ## ===========

    r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]
    v_gcrf = [-4.7432201570; 0.7905364970; 5.5337557270]

    ## DCM
    ## ---

    D_TOD_GCRF = rECItoECI(GCRF(), TOD(), JD_UTC, eop_iau1980)

    r_tod = D_TOD_GCRF*r_gcrf
    v_tod = D_TOD_GCRF*v_gcrf

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    @test v_tod[1] ≈ -4.7460883850  atol=1e-7
    @test v_tod[2] ≈ +0.7860783240  atol=1e-7
    @test v_tod[3] ≈ +5.5319312880  atol=1e-7

    ## Quaternion
    ## ----------

    q_TOD_GCRF = rECItoECI(Quaternion, GCRF(), TOD(), JD_UTC, eop_iau1980)

    r_tod = vect(conj(q_TOD_GCRF)*r_gcrf*q_TOD_GCRF)
    v_tod = vect(conj(q_TOD_GCRF)*v_gcrf*q_TOD_GCRF)

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    @test v_tod[1] ≈ -4.7460883850  atol=1e-7
    @test v_tod[2] ≈ +0.7860783240  atol=1e-7
    @test v_tod[3] ≈ +5.5319312880  atol=1e-7

    ## TOD => GCRF
    ## ===========

    r_tod = [5094.51620300; 6127.36527840; 6380.34453270]
    v_tod = [-4.7460883850; 0.7860783240; 5.5319312880]

    ## DCM
    ## ---

    D_GCRF_TOD = rECItoECI(TOD(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = D_GCRF_TOD*r_tod
    v_gcrf = D_GCRF_TOD*v_tod

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4

    @test v_gcrf[1] ≈ -4.7432201570  atol=1e-7
    @test v_gcrf[2] ≈ +0.7905364970  atol=1e-7
    @test v_gcrf[3] ≈ +5.5337557270  atol=1e-7

    ## Quaternion
    ## ----------

    q_GCRF_TOD = rECItoECI(Quaternion, TOD(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = vect(conj(q_GCRF_TOD)*r_tod*q_GCRF_TOD)
    v_gcrf = vect(conj(q_GCRF_TOD)*v_tod*q_GCRF_TOD)

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4

    @test v_gcrf[1] ≈ -4.7432201570  atol=1e-7
    @test v_gcrf[2] ≈ +0.7905364970  atol=1e-7
    @test v_gcrf[3] ≈ +5.5337557270  atol=1e-7
end

# GCRF <=> TEME
# =============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_gcrf = 5102.50895790   i + 6123.01140070   j + 6378.13692820   k [km]
#   v_gcrf =   -4.7432201570 i +    0.7905364970 j +    5.5337557270 k [km/s]
#
# one gets:
#
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#   v_teme =   -4.7461314870 i +    0.7858180410 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Function rECItoECI GCRF <=> TEME" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## GCRF => TEME
    ## ===========

    r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]
    v_gcrf = [-4.7432201570; 0.7905364970; 5.5337557270]

    ## DCM
    ## ---

    D_TEME_GCRF = rECItoECI(GCRF(), TEME(), JD_UTC, eop_iau1980)

    r_teme = D_TEME_GCRF*r_gcrf
    v_teme = D_TEME_GCRF*v_gcrf

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    @test v_teme[1] ≈ -4.7461314870  atol=1e-7
    @test v_teme[2] ≈ +0.7858180410  atol=1e-7
    @test v_teme[3] ≈ +5.5319312880  atol=1e-7

    ## Quaternion
    ## ----------

    q_TEME_GCRF = rECItoECI(Quaternion, GCRF(), TEME(), JD_UTC, eop_iau1980)

    r_teme = vect(conj(q_TEME_GCRF)*r_gcrf*q_TEME_GCRF)
    v_teme = vect(conj(q_TEME_GCRF)*v_gcrf*q_TEME_GCRF)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    @test v_teme[1] ≈ -4.7461314870  atol=1e-7
    @test v_teme[2] ≈ +0.7858180410  atol=1e-7
    @test v_teme[3] ≈ +5.5319312880  atol=1e-7

    ## TEME => GCRF
    ## ============

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]
    v_teme = [-4.7461314870; 0.7858180410; 5.5319312880]

    ## DCM
    ## ---

    D_GCRF_TEME = rECItoECI(TEME(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = D_GCRF_TEME*r_teme
    v_gcrf = D_GCRF_TEME*v_teme

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4

    @test v_gcrf[1] ≈ -4.7432201570  atol=1e-7
    @test v_gcrf[2] ≈ +0.7905364970  atol=1e-7
    @test v_gcrf[3] ≈ +5.5337557270  atol=1e-7

    ## Quaternion
    ## ----------

    q_GCRF_TEME = rECItoECI(Quaternion, TEME(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = vect(conj(q_GCRF_TEME)*r_teme*q_GCRF_TEME)
    v_gcrf = vect(conj(q_GCRF_TEME)*v_teme*q_GCRF_TEME)

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4

    @test v_gcrf[1] ≈ -4.7432201570  atol=1e-7
    @test v_gcrf[2] ≈ +0.7905364970  atol=1e-7
    @test v_gcrf[3] ≈ +5.5337557270  atol=1e-7
end

# J2000 <=> MOD
# =============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# According to this example and Table 3-6, using:
#
#   UTC     = April 6, 2004, 07:51:28.386009
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#   v_j2000 =  -4.7432196000 i + 0.7905366000    j + 5.5337561900    k [km/s]
#
# one gets:
#
#   r_mod = 5094.02837450   i + 6127.87081640   j + 6380.24851640   k [km]
#   v_mod =   -4.7462630520 i +    0.7860140450 j +    5.5317905620 k [km/s]
#
################################################################################

@testset "Function rECItoECI J2000 <=> MOD" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## J2000 => MOD
    ## ============

    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]
    v_j2000 = [-4.7432196000; 0.7905366000; 5.5337561900]

    ## DCM
    ## ---

    D_MOD_J2000 = rECItoECI(J2000(), MOD(), JD_UTC, eop_iau1980)

    r_mod = D_MOD_J2000*r_j2000
    v_mod = D_MOD_J2000*v_j2000

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    @test v_mod[1] ≈ -4.7462630520  atol=1e-7
    @test v_mod[2] ≈ +0.7860140450  atol=1e-7
    @test v_mod[3] ≈ +5.5317905620  atol=1e-7

    ## Quaternion
    ## ----------

    q_MOD_J2000 = rECItoECI(Quaternion, J2000(), MOD(), JD_UTC, eop_iau1980)

    r_mod = vect(conj(q_MOD_J2000)*r_j2000*q_MOD_J2000)
    v_mod = vect(conj(q_MOD_J2000)*v_j2000*q_MOD_J2000)

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    @test v_mod[1] ≈ -4.7462630520  atol=1e-7
    @test v_mod[2] ≈ +0.7860140450  atol=1e-7
    @test v_mod[3] ≈ +5.5317905620  atol=1e-7

    ## MOD => J2000
    ## ============

    r_mod = [5094.02837450; 6127.87081640; 6380.24851640]
    v_mod = [-4.7462630520; 0.7860140450; 5.5317905620]

    ## DCM
    ## ---

    D_J2000_MOD = rECItoECI(MOD(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = D_J2000_MOD*r_mod
    v_j2000 = D_J2000_MOD*v_mod

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    @test v_j2000[1] ≈ -4.7432196000  atol=1e-7
    @test v_j2000[2] ≈ +0.7905366000  atol=1e-7
    @test v_j2000[3] ≈ +5.5337561900  atol=1e-7

    ## Quaternion
    ## ----------

    q_J2000_MOD = rECItoECI(Quaternion, MOD(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = vect(conj(q_J2000_MOD)*r_mod*q_J2000_MOD)
    v_j2000 = vect(conj(q_J2000_MOD)*v_mod*q_J2000_MOD)

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    @test v_j2000[1] ≈ -4.7432196000  atol=1e-7
    @test v_j2000[2] ≈ +0.7905366000  atol=1e-7
    @test v_j2000[3] ≈ +5.5337561900  atol=1e-7
end

# J2000 <=> TOD
# =============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#   v_j2000 =  -4.7432196000 i + 0.7905366000    j + 5.5337561900    k [km/s]
#
# one gets:
#
#   r_tod = 5094.51620300   i + 6127.36527840   j + 6380.34453270   k [km]
#   v_tod =   -4.7460883850 i +    0.7860783240 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Function rECItoECI J2000 <=> TOD" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## J2000 => TOD
    ## ============

    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]
    v_j2000 = [-4.7432196000; 0.7905366000; 5.5337561900]

    ## DCM
    ## ---

    D_TOD_J2000 = rECItoECI(J2000(), TOD(), JD_UTC, eop_iau1980)

    r_tod = D_TOD_J2000*r_j2000
    v_tod = D_TOD_J2000*v_j2000

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    @test v_tod[1] ≈ -4.7460883850  atol=1e-7
    @test v_tod[2] ≈ +0.7860783240  atol=1e-7
    @test v_tod[3] ≈ +5.5319312880  atol=1e-7

    ## Quaternion
    ## ----------

    q_TOD_J2000 = rECItoECI(Quaternion, J2000(), TOD(), JD_UTC, eop_iau1980)

    r_tod = vect(conj(q_TOD_J2000)*r_j2000*q_TOD_J2000)
    v_tod = vect(conj(q_TOD_J2000)*v_j2000*q_TOD_J2000)

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    @test v_tod[1] ≈ -4.7460883850  atol=1e-7
    @test v_tod[2] ≈ +0.7860783240  atol=1e-7
    @test v_tod[3] ≈ +5.5319312880  atol=1e-7

    ## TOD => J2000
    ## ============

    r_tod = [5094.51620300; 6127.36527840; 6380.34453270]
    v_tod = [-4.7460883850; 0.7860783240; 5.5319312880]

    ## DCM
    ## ---

    D_J2000_TOD = rECItoECI(TOD(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = D_J2000_TOD*r_tod
    v_j2000 = D_J2000_TOD*v_tod

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    @test v_j2000[1] ≈ -4.7432196000  atol=1e-7
    @test v_j2000[2] ≈ +0.7905366000  atol=1e-7
    @test v_j2000[3] ≈ +5.5337561900  atol=1e-7

    ## Quaternion
    ## ----------

    q_J2000_TOD = rECItoECI(Quaternion, TOD(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = vect(conj(q_J2000_TOD)*r_tod*q_J2000_TOD)
    v_j2000 = vect(conj(q_J2000_TOD)*v_tod*q_J2000_TOD)

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    @test v_j2000[1] ≈ -4.7432196000  atol=1e-7
    @test v_j2000[2] ≈ +0.7905366000  atol=1e-7
    @test v_j2000[3] ≈ +5.5337561900  atol=1e-7
end

# J2000 <=> TEME
# ==============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#   v_j2000 =  -4.7432196000 i + 0.7905366000    j + 5.5337561900    k [km/s]
#
# one gets:
#
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#   v_teme =   -4.7461314870 i +    0.7858180410 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Function rECItoECI J2000 <=> TEME" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## J2000 => TEME
    ## =============

    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]
    v_j2000 = [-4.7432196000; 0.7905366000; 5.5337561900]

    ## DCM
    ## ---

    D_TEME_J2000 = rECItoECI(J2000(), TEME(), JD_UTC, eop_iau1980)

    r_teme = D_TEME_J2000*r_j2000
    v_teme = D_TEME_J2000*v_j2000

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    @test v_teme[1] ≈ -4.7461314870  atol=1e-7
    @test v_teme[2] ≈ +0.7858180410  atol=1e-7
    @test v_teme[3] ≈ +5.5319312880  atol=1e-7

    ## Quaternion
    ## ----------

    q_TEME_J2000 = rECItoECI(Quaternion, J2000(), TEME(), JD_UTC)

    r_teme = vect(conj(q_TEME_J2000)*r_j2000*q_TEME_J2000)
    v_teme = vect(conj(q_TEME_J2000)*v_j2000*q_TEME_J2000)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    @test v_teme[1] ≈ -4.7461314870  atol=1e-7
    @test v_teme[2] ≈ +0.7858180410  atol=1e-7
    @test v_teme[3] ≈ +5.5319312880  atol=1e-7

    ## TEME => J2000
    ## =============

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]
    v_teme = [-4.7461314870; 0.7858180410; 5.5319312880]

    ## DCM
    ## ---

    D_J2000_TEME = rECItoECI(TEME(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = D_J2000_TEME*r_teme
    v_j2000 = D_J2000_TEME*v_teme

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    @test v_j2000[1] ≈ -4.7432196000  atol=1e-7
    @test v_j2000[2] ≈ +0.7905366000  atol=1e-7
    @test v_j2000[3] ≈ +5.5337561900  atol=1e-7

    ## Quaternion
    ## ----------

    q_J2000_TEME = rECItoECI(Quaternion, TEME(), J2000(), JD_UTC)

    r_j2000 = vect(conj(q_J2000_TEME)*r_teme*q_J2000_TEME)
    v_j2000 = vect(conj(q_J2000_TEME)*v_teme*q_J2000_TEME)

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    @test v_j2000[1] ≈ -4.7432196000  atol=1e-7
    @test v_j2000[2] ≈ +0.7905366000  atol=1e-7
    @test v_j2000[3] ≈ +5.5337561900  atol=1e-7
end

# MOD <=> TOD
# ===========

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_mod = 5094.02837450   i + 6127.87081640   j + 6380.24851640   k [km]
#   v_mod =   -4.7462630520 i +    0.7860140450 j +    5.5317905620 k [km/s]
#
# one gets:
#
#   r_tod = 5094.51620300   i + 6127.36527840   j + 6380.34453270   k [km]
#   v_tod =   -4.7460883850 i +    0.7860783240 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Function rECItoECI MOD <=> TOD" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## MOD => TOD
    ## ==========

    r_mod = [5094.02837450; 6127.87081640; 6380.24851640]
    v_mod = [-4.7462630520; 0.7860140450; 5.5317905620]

    ## DCM
    ## ---

    D_TOD_MOD = rECItoECI(MOD(), JD_UTC, TOD(), JD_UTC, eop_iau1980)

    r_tod = D_TOD_MOD*r_mod
    v_tod = D_TOD_MOD*v_mod

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    @test v_tod[1] ≈ -4.7460883850  atol=1e-7
    @test v_tod[2] ≈ +0.7860783240  atol=1e-7
    @test v_tod[3] ≈ +5.5319312880  atol=1e-7

    ## Quaternion
    ## ----------

    q_TOD_MOD = rECItoECI(Quaternion, MOD(), JD_UTC, TOD(), JD_UTC, eop_iau1980)

    r_tod = vect(conj(q_TOD_MOD)*r_mod*q_TOD_MOD)
    v_tod = vect(conj(q_TOD_MOD)*v_mod*q_TOD_MOD)

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    @test v_tod[1] ≈ -4.7460883850  atol=1e-7
    @test v_tod[2] ≈ +0.7860783240  atol=1e-7
    @test v_tod[3] ≈ +5.5319312880  atol=1e-7

    ## TOD => MOD
    ## ==========

    r_tod = [5094.51620300; 6127.36527840; 6380.34453270]
    v_tod = [-4.7460883850; 0.7860783240; 5.5319312880]

    ## DCM
    ## ---

    D_MOD_TOD = rECItoECI(TOD(), JD_UTC, MOD(), JD_UTC, eop_iau1980)

    r_mod = D_MOD_TOD*r_tod
    v_mod = D_MOD_TOD*v_tod

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    @test v_mod[1] ≈ -4.7462630520  atol=1e-7
    @test v_mod[2] ≈ +0.7860140450  atol=1e-7
    @test v_mod[3] ≈ +5.5317905620  atol=1e-7

    ## Quaternion
    ## ----------

    q_MOD_TOD = rECItoECI(Quaternion, TOD(), JD_UTC, MOD(), JD_UTC, eop_iau1980)

    r_mod = vect(conj(q_MOD_TOD)*r_tod*q_MOD_TOD)
    v_mod = vect(conj(q_MOD_TOD)*v_tod*q_MOD_TOD)

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    @test v_mod[1] ≈ -4.7462630520  atol=1e-7
    @test v_mod[2] ≈ +0.7860140450  atol=1e-7
    @test v_mod[3] ≈ +5.5317905620  atol=1e-7
end

# MOD <=> TEME
# ============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_mod = 5094.02837450   i + 6127.87081640   j + 6380.24851640   k [km]
#   v_mod =   -4.7462630520 i +    0.7860140450 j +    5.5317905620 k [km/s]
#
# one gets:
#
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#   v_teme =   -4.7461314870 i +    0.7858180410 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Function rECItoECI MOD <=> TEME" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## MOD => TEME
    ## ===========

    r_mod = [5094.02837450; 6127.87081640; 6380.24851640]
    v_mod = [-4.7462630520; 0.7860140450; 5.5317905620]

    ## DCM
    ## ---

    D_TEME_MOD = rECItoECI(MOD(), JD_UTC, TEME(), JD_UTC, eop_iau1980)

    r_teme = D_TEME_MOD*r_mod
    v_teme = D_TEME_MOD*v_mod

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    @test v_teme[1] ≈ -4.7461314870  atol=1e-7
    @test v_teme[2] ≈ +0.7858180410  atol=1e-7
    @test v_teme[3] ≈ +5.5319312880  atol=1e-7

    ## Quaternion
    ## ----------

    q_TEME_MOD = rECItoECI(Quaternion, MOD(), JD_UTC, TEME(), JD_UTC, eop_iau1980)

    r_teme = vect(conj(q_TEME_MOD)*r_mod*q_TEME_MOD)
    v_teme = vect(conj(q_TEME_MOD)*v_mod*q_TEME_MOD)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    @test v_teme[1] ≈ -4.7461314870  atol=1e-7
    @test v_teme[2] ≈ +0.7858180410  atol=1e-7
    @test v_teme[3] ≈ +5.5319312880  atol=1e-7

    ## TEME => MOD
    ## ===========

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]
    v_teme = [-4.7461314870; 0.7858180410; 5.5319312880]

    ## DCM
    ## ---

    D_MOD_TEME = rECItoECI(TEME(), JD_UTC, MOD(), JD_UTC, eop_iau1980)

    r_mod = D_MOD_TEME*r_teme
    v_mod = D_MOD_TEME*v_teme

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    @test v_mod[1] ≈ -4.7462630520  atol=1e-7
    @test v_mod[2] ≈ +0.7860140450  atol=1e-7
    @test v_mod[3] ≈ +5.5317905620  atol=1e-7

    ## Quaternion
    ## ----------

    q_MOD_TEME = rECItoECI(Quaternion, TEME(), JD_UTC, MOD(), JD_UTC, eop_iau1980)

    r_mod = vect(conj(q_MOD_TEME)*r_teme*q_MOD_TEME)
    v_mod = vect(conj(q_MOD_TEME)*v_teme*q_MOD_TEME)

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    @test v_mod[1] ≈ -4.7462630520  atol=1e-7
    @test v_mod[2] ≈ +0.7860140450  atol=1e-7
    @test v_mod[3] ≈ +5.5317905620  atol=1e-7
end

# TOD <=> TEME
# ============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_tod = 5094.51620300   i + 6127.36527840   j + 6380.34453270   k [km]
#   v_tod =   -4.7460883850 i +    0.7860783240 j +    5.5319312880 k [km/s]
#
# one gets:
#
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#   v_teme =   -4.7461314870 i +    0.7858180410 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Function rECItoECI TOD <=> TEME" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## TOD => TEME
    ## ===========

    r_tod = [5094.51620300; 6127.36527840; 6380.34453270]
    v_tod = [-4.7460883850; 0.7860783240; 5.5319312880]

    ## DCM
    ## ---

    D_TEME_TOD = rECItoECI(TOD(), JD_UTC, TEME(), JD_UTC, eop_iau1980)

    r_teme = D_TEME_TOD*r_tod
    v_teme = D_TEME_TOD*v_tod

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    @test v_teme[1] ≈ -4.7461314870  atol=1e-7
    @test v_teme[2] ≈ +0.7858180410  atol=1e-7
    @test v_teme[3] ≈ +5.5319312880  atol=1e-7

    ## Quaternion
    ## ----------

    q_TEME_TOD = rECItoECI(Quaternion, TOD(), JD_UTC, TEME(), JD_UTC, eop_iau1980)

    r_teme = vect(conj(q_TEME_TOD)*r_tod*q_TEME_TOD)
    v_teme = vect(conj(q_TEME_TOD)*v_tod*q_TEME_TOD)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    @test v_teme[1] ≈ -4.7461314870  atol=1e-7
    @test v_teme[2] ≈ +0.7858180410  atol=1e-7
    @test v_teme[3] ≈ +5.5319312880  atol=1e-7

    ## TEME => TOD
    ## ===========

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]
    v_teme = [-4.7461314870; 0.7858180410; 5.5319312880]

    ## DCM
    ## ---

    D_TOD_TEME = rECItoECI(TEME(), JD_UTC, TOD(), JD_UTC, eop_iau1980)

    r_tod = D_TOD_TEME*r_teme
    v_tod = D_TOD_TEME*v_teme

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    @test v_tod[1] ≈ -4.7460883850  atol=1e-7
    @test v_tod[2] ≈ +0.7860783240  atol=1e-7
    @test v_tod[3] ≈ +5.5319312880  atol=1e-7

    ## Quaternion
    ## ----------

    q_TOD_TEME = rECItoECI(Quaternion, TEME(), JD_UTC, TOD(), JD_UTC, eop_iau1980)

    r_tod = vect(conj(q_TOD_TEME)*r_teme*q_TOD_TEME)
    v_tod = vect(conj(q_TOD_TEME)*v_teme*q_TOD_TEME)

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    @test v_tod[1] ≈ -4.7460883850  atol=1e-7
    @test v_tod[2] ≈ +0.7860783240  atol=1e-7
    @test v_tod[3] ≈ +5.5319312880  atol=1e-7
end

################################################################################
#                                IAU-2006/2010
################################################################################

## GCRF <=> CIRS
## =============

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
# NOTE: It seems that the results in the Example 3-14 is computed **without**
# the `dX` and `dY` corrections, whereas in the Table 3-6 they are computed
# **with** the corrections.
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_cirs = -5100.01840470   i + 6122.78636480   j + 6380.34453270   k [km]
#   v_cirs =    -4.7453803300 i -    0.7903414530 j +    5.5319312880 k [km/s]
#
# one gets the following (this is the result in Table 3-6):
#
#   r_gcrf = 5102.50895290  i + 6123.01139910 j + 6378.13693380 k [km]
#   v_gcrf =  -4.7432201610 i + 0.7905364950  j + 5.5337557240  k [km/s]
#
################################################################################

@testset "Function rECItoECI GCRF <=> CIRS" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## GCRF => CIRS
    ## ============

    r_gcrf = [5102.50895290; 6123.01139910; 6378.13693380]
    v_gcrf = [-4.7432201610; 0.7905364950; 5.5337557240]

    ## DCM
    ## ---

    D_CIRS_GCRF = rECItoECI(GCRF(), CIRS(), JD_UTC, eop_iau2000a)

    r_cirs = D_CIRS_GCRF*r_gcrf
    v_cirs = D_CIRS_GCRF*v_gcrf

    @test r_cirs[1] ≈ +5100.01840470 atol=1e-4
    @test r_cirs[2] ≈ +6122.78636480 atol=1e-4
    @test r_cirs[3] ≈ +6380.34453270 atol=1e-4

    @test v_cirs[1] ≈ -4.7453803300  atol=1e-7
    @test v_cirs[2] ≈ +0.7903414530  atol=1e-7
    @test v_cirs[3] ≈ +5.5319312880  atol=1e-7

    ## Quaternion
    ## ----------

    q_CIRS_GCRF = rECItoECI(Quaternion, GCRF(), CIRS(), JD_UTC, eop_iau2000a)

    r_cirs = vect(conj(q_CIRS_GCRF)*r_gcrf*q_CIRS_GCRF)
    v_cirs = vect(conj(q_CIRS_GCRF)*v_gcrf*q_CIRS_GCRF)

    @test r_cirs[1] ≈ +5100.01840470 atol=1e-4
    @test r_cirs[2] ≈ +6122.78636480 atol=1e-4
    @test r_cirs[3] ≈ +6380.34453270 atol=1e-4

    @test v_cirs[1] ≈ -4.7453803300  atol=1e-7
    @test v_cirs[2] ≈ +0.7903414530  atol=1e-7
    @test v_cirs[3] ≈ +5.5319312880  atol=1e-7

    ## CIRS => GCRF
    ## ============

    r_cirs = [+5100.01840470; +6122.78636480; +6380.34453270]
    v_cirs = [-4.7453803300; +0.7903414530; +5.5319312880]

    ## DCM
    ## ---

    D_GCRF_CIRS = rECItoECI(CIRS(), GCRF(), JD_UTC, eop_iau2000a)

    r_gcrf = D_GCRF_CIRS*r_cirs
    v_gcrf = D_GCRF_CIRS*v_cirs

    @test r_gcrf[1] ≈ +5102.50895290 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01139910 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13693380 atol=1e-4

    @test v_gcrf[1] ≈ -4.7432201610  atol=1e-7
    @test v_gcrf[2] ≈ +0.7905364950  atol=1e-7
    @test v_gcrf[3] ≈ +5.5337557240  atol=1e-7

    ## Quaternion
    ## ----------

    q_GCRF_CIRS = rECItoECI(Quaternion, CIRS(), GCRF(), JD_UTC, eop_iau2000a)

    r_gcrf = vect(conj(q_GCRF_CIRS)*r_cirs*q_GCRF_CIRS)
    v_gcrf = vect(conj(q_GCRF_CIRS)*v_cirs*q_GCRF_CIRS)

    @test r_gcrf[1] ≈ +5102.50895290 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01139910 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13693380 atol=1e-4

    @test v_gcrf[1] ≈ -4.7432201610  atol=1e-7
    @test v_gcrf[2] ≈ +0.7905364950  atol=1e-7
    @test v_gcrf[3] ≈ +5.5337557240  atol=1e-7

end
