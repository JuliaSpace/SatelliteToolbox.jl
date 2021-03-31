# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to ECI to ECI transformations.
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

# Get the current EOP Data.
#
# TODO: The EOP data obtained from IERS website does not match the values in
# the examples in [1]. However, this should be enough, because 1) the
# individual functions at the low level are tested using the same values of
# [1], and 2) the difference is smaller than 30 cm.

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
#
# one gets:
#
#   r_j2000 = 5102.50960000   i + 6123.01152000   j + 6378.13630000   k [km]
#
################################################################################

@testset "Function rECItoECI GCRF <=> J2000" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## GCRF => J2000
    ## =============

    r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]

    ## DCM
    ## ---

    D_J2000_GCRF = rECItoECI(GCRF(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = D_J2000_GCRF*r_gcrf

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    ## Quaternion
    ## ----------

    q_J2000_GCRF = rECItoECI(Quaternion, GCRF(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = vect(conj(q_J2000_GCRF)*r_gcrf*q_J2000_GCRF)

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    ## J2000 => GCRF
    ## =============

    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]

    ## DCM
    ## ---

    D_GCRF_J2000 = rECItoECI(J2000(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = D_GCRF_J2000*r_j2000

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4

    ## Quaternion
    ## ----------

    q_GCRF_J2000 = rECItoECI(Quaternion, J2000(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = vect(conj(q_GCRF_J2000)*r_j2000*q_GCRF_J2000)

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4
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
#
# one gets:
#
#   r_mod = 5094.02837450   i + 6127.87081640   j + 6380.24851640   k [km]
#
################################################################################

@testset "Function rECItoECI GCRF <=> MOD" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## GCRF => MOD
    ## ===========

    r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]

    ## DCM
    ## ---

    D_MOD_GCRF = rECItoECI(GCRF(), MOD(), JD_UTC, eop_iau1980)

    r_mod = D_MOD_GCRF*r_gcrf

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    ## Quaternion
    ## ----------

    q_MOD_GCRF = rECItoECI(Quaternion, GCRF(), MOD(), JD_UTC, eop_iau1980)

    r_mod = vect(conj(q_MOD_GCRF)*r_gcrf*q_MOD_GCRF)

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    ## MOD => GCRF
    ## ===========

    r_mod = [5094.02837450; 6127.87081640; 6380.24851640]
    v_mod = [-4.7462630520; 0.7860140450; 5.5317905620]

    ## DCM
    ## ---

    D_GCRF_MOD = rECItoECI(MOD(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = D_GCRF_MOD*r_mod

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4

    ## Quaternion
    ## ----------

    q_GCRF_MOD = rECItoECI(Quaternion, MOD(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = vect(conj(q_GCRF_MOD)*r_mod*q_GCRF_MOD)

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4
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
#
# one gets:
#
#   r_tod = 5094.51620300   i + 6127.36527840   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECI GCRF <=> TOD" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## GCRF => TOD
    ## ===========

    r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]

    ## DCM
    ## ---

    D_TOD_GCRF = rECItoECI(GCRF(), TOD(), JD_UTC, eop_iau1980)

    r_tod = D_TOD_GCRF*r_gcrf

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    ## Quaternion
    ## ----------

    q_TOD_GCRF = rECItoECI(Quaternion, GCRF(), TOD(), JD_UTC, eop_iau1980)

    r_tod = vect(conj(q_TOD_GCRF)*r_gcrf*q_TOD_GCRF)

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    ## TOD => GCRF
    ## ===========

    r_tod = [5094.51620300; 6127.36527840; 6380.34453270]

    ## DCM
    ## ---

    D_GCRF_TOD = rECItoECI(TOD(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = D_GCRF_TOD*r_tod

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4

    ## Quaternion
    ## ----------

    q_GCRF_TOD = rECItoECI(Quaternion, TOD(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = vect(conj(q_GCRF_TOD)*r_tod*q_GCRF_TOD)

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4
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
#
# one gets:
#
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECI GCRF <=> TEME" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## GCRF => TEME
    ## ===========

    r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]

    ## DCM
    ## ---

    D_TEME_GCRF = rECItoECI(GCRF(), TEME(), JD_UTC, eop_iau1980)

    r_teme = D_TEME_GCRF*r_gcrf

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    ## Quaternion
    ## ----------

    q_TEME_GCRF = rECItoECI(Quaternion, GCRF(), TEME(), JD_UTC, eop_iau1980)

    r_teme = vect(conj(q_TEME_GCRF)*r_gcrf*q_TEME_GCRF)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    ## TEME => GCRF
    ## ============

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]

    ## DCM
    ## ---

    D_GCRF_TEME = rECItoECI(TEME(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = D_GCRF_TEME*r_teme

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4

    ## Quaternion
    ## ----------

    q_GCRF_TEME = rECItoECI(Quaternion, TEME(), GCRF(), JD_UTC, eop_iau1980)

    r_gcrf = vect(conj(q_GCRF_TEME)*r_teme*q_GCRF_TEME)

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-4
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
#
# one gets:
#
#   r_mod = 5094.02837450   i + 6127.87081640   j + 6380.24851640   k [km]
#
# If not EOP correction is used, then we get:
#
#   r_mod = 5094.02901670   i + 6127.87093630   j + 6380.24788850   k [km]
#
################################################################################

@testset "Function rECItoECI J2000 <=> MOD" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## J2000 => MOD
    ## ============

    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]

    ## DCM
    ## ---

    D_MOD_J2000 = rECItoECI(J2000(), MOD(), JD_UTC, eop_iau1980)

    r_mod = D_MOD_J2000*r_j2000

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    ## Quaternion
    ## ----------

    q_MOD_J2000 = rECItoECI(Quaternion, J2000(), MOD(), JD_UTC, eop_iau1980)

    r_mod = vect(conj(q_MOD_J2000)*r_j2000*q_MOD_J2000)

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    ## MOD => J2000
    ## ============

    r_mod = [5094.02837450; 6127.87081640; 6380.24851640]

    ## DCM
    ## ---

    D_J2000_MOD = rECItoECI(MOD(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = D_J2000_MOD*r_mod

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    ## Quaternion
    ## ----------

    q_J2000_MOD = rECItoECI(Quaternion, MOD(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = vect(conj(q_J2000_MOD)*r_mod*q_J2000_MOD)

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    # No EOP Corrections
    # ==========================================================================

    ## J2000 => MOD
    ## ============

    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]

    ## DCM
    ## ---

    D_MOD_J2000 = rECItoECI(J2000(), MOD(), JD_UTC)

    r_mod = D_MOD_J2000*r_j2000

    @test r_mod[1] ≈ +5094.02901670 atol=1e-7
    @test r_mod[2] ≈ +6127.87093630 atol=1e-7
    @test r_mod[3] ≈ +6380.24788850 atol=1e-7

    ## Quaternion
    ## ----------

    q_MOD_J2000 = rECItoECI(Quaternion, J2000(), MOD(), JD_UTC)

    r_mod = vect(conj(q_MOD_J2000)*r_j2000*q_MOD_J2000)

    @test r_mod[1] ≈ +5094.02901670 atol=1e-7
    @test r_mod[2] ≈ +6127.87093630 atol=1e-7
    @test r_mod[3] ≈ +6380.24788850 atol=1e-7

    ## MOD => J2000
    ## ============

    r_mod = [5094.02901670; 6127.87093630; 6380.24788850]

    ## DCM
    ## ---

    D_J2000_MOD = rECItoECI(MOD(), J2000(), JD_UTC)

    r_j2000 = D_J2000_MOD*r_mod

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-7
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-7
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-7

    ## Quaternion
    ## ----------

    q_J2000_MOD = rECItoECI(Quaternion, MOD(), J2000(), JD_UTC)

    r_j2000 = vect(conj(q_J2000_MOD)*r_mod*q_J2000_MOD)

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-7
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-7
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-7
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
#
# one gets:
#
#   r_tod = 5094.51620300   i + 6127.36527840   j + 6380.34453270   k [km]
#
# If not EOP correction is used, then we get:
#
#   r_tod = 5094.51478040   i + 6127.36646120   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECI J2000 <=> TOD" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## J2000 => TOD
    ## ============

    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]

    ## DCM
    ## ---

    D_TOD_J2000 = rECItoECI(J2000(), TOD(), JD_UTC, eop_iau1980)

    r_tod = D_TOD_J2000*r_j2000

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    ## Quaternion
    ## ----------

    q_TOD_J2000 = rECItoECI(Quaternion, J2000(), TOD(), JD_UTC, eop_iau1980)

    r_tod = vect(conj(q_TOD_J2000)*r_j2000*q_TOD_J2000)

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    ## TOD => J2000
    ## ============

    r_tod = [5094.51620300; 6127.36527840; 6380.34453270]

    ## DCM
    ## ---

    D_J2000_TOD = rECItoECI(TOD(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = D_J2000_TOD*r_tod

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    ## Quaternion
    ## ----------

    q_J2000_TOD = rECItoECI(Quaternion, TOD(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = vect(conj(q_J2000_TOD)*r_tod*q_J2000_TOD)

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    # No EOP corrections
    # ==========================================================================

    ## J2000 => TOD
    ## ============

    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]

    ## DCM
    ## ---

    D_TOD_J2000 = rECItoECI(J2000(), TOD(), JD_UTC)

    r_tod = D_TOD_J2000*r_j2000

    @test r_tod[1] ≈ +5094.51478040 atol=1e-7
    @test r_tod[2] ≈ +6127.36646120 atol=1e-7
    @test r_tod[3] ≈ +6380.34453270 atol=1e-7

    ## Quaternion
    ## ----------

    q_TOD_J2000 = rECItoECI(Quaternion, J2000(), TOD(), JD_UTC)

    r_tod = vect(conj(q_TOD_J2000)*r_j2000*q_TOD_J2000)

    @test r_tod[1] ≈ +5094.51478040 atol=1e-7
    @test r_tod[2] ≈ +6127.36646120 atol=1e-7
    @test r_tod[3] ≈ +6380.34453270 atol=1e-7

    ## TOD => J2000
    ## ============

    r_tod = [5094.51478040; 6127.36646120; 6380.34453270]

    ## DCM
    ## ---

    D_J2000_TOD = rECItoECI(TOD(), J2000(), JD_UTC)

    r_j2000 = D_J2000_TOD*r_tod

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-7
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-7
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-7

    ## Quaternion
    ## ----------

    q_J2000_TOD = rECItoECI(Quaternion, TOD(), J2000(), JD_UTC)

    r_j2000 = vect(conj(q_J2000_TOD)*r_tod*q_J2000_TOD)

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-7
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-7
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-7
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
#
# one gets:
#
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECI J2000 <=> TEME" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## J2000 => TEME
    ## =============

    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]

    ## DCM
    ## ---

    D_TEME_J2000 = rECItoECI(J2000(), TEME(), JD_UTC, eop_iau1980)

    r_teme = D_TEME_J2000*r_j2000

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    ## Quaternion
    ## ----------

    q_TEME_J2000 = rECItoECI(Quaternion, J2000(), TEME(), JD_UTC)

    r_teme = vect(conj(q_TEME_J2000)*r_j2000*q_TEME_J2000)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    ## TEME => J2000
    ## =============

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]

    ## DCM
    ## ---

    D_J2000_TEME = rECItoECI(TEME(), J2000(), JD_UTC, eop_iau1980)

    r_j2000 = D_J2000_TEME*r_teme

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4

    ## Quaternion
    ## ----------

    q_J2000_TEME = rECItoECI(Quaternion, TEME(), J2000(), JD_UTC)

    r_j2000 = vect(conj(q_J2000_TEME)*r_teme*q_J2000_TEME)

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-4
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-4
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-4
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
#
# one gets:
#
#   r_tod = 5094.51620300   i + 6127.36527840   j + 6380.34453270   k [km]
#
# If not EOP correction is used, then we get:
#
#   r_mod = 5094.02901670   i + 6127.87093630   j + 6380.24788850   k [km]
#
#   r_tod = 5094.51478040   i + 6127.36646120   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECI MOD <=> TOD" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## MOD => TOD
    ## ==========

    r_mod = [5094.02837450; 6127.87081640; 6380.24851640]

    ## DCM
    ## ---

    D_TOD_MOD = rECItoECI(MOD(), JD_UTC, TOD(), JD_UTC, eop_iau1980)

    r_tod = D_TOD_MOD*r_mod

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    ## Quaternion
    ## ----------

    q_TOD_MOD = rECItoECI(Quaternion, MOD(), JD_UTC, TOD(), JD_UTC, eop_iau1980)

    r_tod = vect(conj(q_TOD_MOD)*r_mod*q_TOD_MOD)

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    ## TOD => MOD
    ## ==========

    r_tod = [5094.51620300; 6127.36527840; 6380.34453270]

    ## DCM
    ## ---

    D_MOD_TOD = rECItoECI(TOD(), JD_UTC, MOD(), JD_UTC, eop_iau1980)

    r_mod = D_MOD_TOD*r_tod

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    ## Quaternion
    ## ----------

    q_MOD_TOD = rECItoECI(Quaternion, TOD(), JD_UTC, MOD(), JD_UTC, eop_iau1980)

    r_mod = vect(conj(q_MOD_TOD)*r_tod*q_MOD_TOD)

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    # No EOP corrections
    # ==========================================================================

    ## MOD => TOD
    ## ==========

    r_mod = [5094.02901670; 6127.87093630; 6380.24788850]

    ## DCM
    ## ---

    D_TOD_MOD = rECItoECI(MOD(), JD_UTC, TOD(), JD_UTC)

    r_tod = D_TOD_MOD*r_mod

    @test r_tod[1] ≈ +5094.51478040 atol=1e-7
    @test r_tod[2] ≈ +6127.36646120 atol=1e-7
    @test r_tod[3] ≈ +6380.34453270 atol=1e-7

    ## Quaternion
    ## ----------

    q_TOD_MOD = rECItoECI(Quaternion, MOD(), JD_UTC, TOD(), JD_UTC)

    r_tod = vect(conj(q_TOD_MOD)*r_mod*q_TOD_MOD)

    @test r_tod[1] ≈ +5094.51478040 atol=1e-7
    @test r_tod[2] ≈ +6127.36646120 atol=1e-7
    @test r_tod[3] ≈ +6380.34453270 atol=1e-7

    ## TOD => MOD
    ## ==========

    r_tod = [5094.51478040; 6127.36646120; 6380.34453270]

    ## DCM
    ## ---

    D_MOD_TOD = rECItoECI(TOD(), JD_UTC, MOD(), JD_UTC)

    r_mod = D_MOD_TOD*r_tod

    @test r_mod[1] ≈ +5094.02901670 atol=1e-7
    @test r_mod[2] ≈ +6127.87093630 atol=1e-7
    @test r_mod[3] ≈ +6380.24788850 atol=1e-7

    ## Quaternion
    ## ----------

    q_MOD_TOD = rECItoECI(Quaternion, TOD(), JD_UTC, MOD(), JD_UTC)

    r_mod = vect(conj(q_MOD_TOD)*r_tod*q_MOD_TOD)

    @test r_mod[1] ≈ +5094.02901670 atol=1e-7
    @test r_mod[2] ≈ +6127.87093630 atol=1e-7
    @test r_mod[3] ≈ +6380.24788850 atol=1e-7
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
#
# one gets:
#
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#
# If not EOP correction is used, then we get the same result by using:
#
#   r_mod = 5094.02901670   i + 6127.87093630   j + 6380.24788850   k [km]
#
################################################################################

@testset "Function rECItoECI MOD <=> TEME" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## MOD => TEME
    ## ===========

    r_mod = [5094.02837450; 6127.87081640; 6380.24851640]

    ## DCM
    ## ---

    D_TEME_MOD = rECItoECI(MOD(), JD_UTC, TEME(), JD_UTC, eop_iau1980)

    r_teme = D_TEME_MOD*r_mod

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    ## Quaternion
    ## ----------

    q_TEME_MOD = rECItoECI(Quaternion, MOD(), JD_UTC, TEME(), JD_UTC, eop_iau1980)

    r_teme = vect(conj(q_TEME_MOD)*r_mod*q_TEME_MOD)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    ## TEME => MOD
    ## ===========

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]

    ## DCM
    ## ---

    D_MOD_TEME = rECItoECI(TEME(), JD_UTC, MOD(), JD_UTC, eop_iau1980)

    r_mod = D_MOD_TEME*r_teme

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    ## Quaternion
    ## ----------

    q_MOD_TEME = rECItoECI(Quaternion, TEME(), JD_UTC, MOD(), JD_UTC, eop_iau1980)

    r_mod = vect(conj(q_MOD_TEME)*r_teme*q_MOD_TEME)

    @test r_mod[1] ≈ +5094.02837450 atol=1e-4
    @test r_mod[2] ≈ +6127.87081640 atol=1e-4
    @test r_mod[3] ≈ +6380.24851640 atol=1e-4

    # No EOP corrections
    # ==========================================================================

    ## MOD => TEME
    ## ===========

    r_mod = [5094.02901670; 6127.87093630; 6380.24788850]

    ## DCM
    ## ---

    D_TEME_MOD = rECItoECI(MOD(), JD_UTC, TEME(), JD_UTC)

    r_teme = D_TEME_MOD*r_mod

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    ## Quaternion
    ## ----------

    q_TEME_MOD = rECItoECI(Quaternion, MOD(), JD_UTC, TEME(), JD_UTC)

    r_teme = vect(conj(q_TEME_MOD)*r_mod*q_TEME_MOD)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    ## TEME => MOD
    ## ===========

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]

    ## DCM
    ## ---

    D_MOD_TEME = rECItoECI(TEME(), JD_UTC, MOD(), JD_UTC)

    r_mod = D_MOD_TEME*r_teme

    @test r_mod[1] ≈ +5094.02901670 atol=1e-7
    @test r_mod[2] ≈ +6127.87093630 atol=1e-7
    @test r_mod[3] ≈ +6380.24788850 atol=1e-7

    ## Quaternion
    ## ----------

    q_MOD_TEME = rECItoECI(Quaternion, TEME(), JD_UTC, MOD(), JD_UTC)

    r_mod = vect(conj(q_MOD_TEME)*r_teme*q_MOD_TEME)

    @test r_mod[1] ≈ +5094.02901670 atol=1e-7
    @test r_mod[2] ≈ +6127.87093630 atol=1e-7
    @test r_mod[3] ≈ +6380.24788850 atol=1e-7
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
#
# one gets:
#
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECI TOD <=> TEME" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## TOD => TEME
    ## ===========

    r_tod = [5094.51620300; 6127.36527840; 6380.34453270]

    ## DCM
    ## ---

    D_TEME_TOD = rECItoECI(TOD(), JD_UTC, TEME(), JD_UTC, eop_iau1980)

    r_teme = D_TEME_TOD*r_tod

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    ## Quaternion
    ## ----------

    q_TEME_TOD = rECItoECI(Quaternion, TOD(), JD_UTC, TEME(), JD_UTC, eop_iau1980)

    r_teme = vect(conj(q_TEME_TOD)*r_tod*q_TEME_TOD)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    ## TEME => TOD
    ## ===========

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]

    ## DCM
    ## ---

    D_TOD_TEME = rECItoECI(TEME(), JD_UTC, TOD(), JD_UTC, eop_iau1980)

    r_tod = D_TOD_TEME*r_teme

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    ## Quaternion
    ## ----------

    q_TOD_TEME = rECItoECI(Quaternion, TEME(), JD_UTC, TOD(), JD_UTC, eop_iau1980)

    r_tod = vect(conj(q_TOD_TEME)*r_teme*q_TOD_TEME)

    @test r_tod[1] ≈ +5094.51620300 atol=1e-4
    @test r_tod[2] ≈ +6127.36527840 atol=1e-4
    @test r_tod[3] ≈ +6380.34453270 atol=1e-4

    # No EOP corrections
    # ==========================================================================

    ## TOD => TEME
    ## ===========

    r_tod = [5094.51478040; 6127.36646120; 6380.34453270]

    ## DCM
    ## ---

    D_TEME_TOD = rECItoECI(TOD(), JD_UTC, TEME(), JD_UTC)

    r_teme = D_TEME_TOD*r_tod

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    ## Quaternion
    ## ----------

    q_TEME_TOD = rECItoECI(Quaternion, TOD(), JD_UTC, TEME(), JD_UTC)

    r_teme = vect(conj(q_TEME_TOD)*r_tod*q_TEME_TOD)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-4
    @test r_teme[2] ≈ +6127.64465950 atol=1e-4
    @test r_teme[3] ≈ +6380.34453270 atol=1e-4

    ## TEME => TOD
    ## ===========

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]

    ## DCM
    ## ---

    D_TOD_TEME = rECItoECI(TEME(), JD_UTC, TOD(), JD_UTC)

    r_tod = D_TOD_TEME*r_teme

    @test r_tod[1] ≈ +5094.51478040 atol=1e-7
    @test r_tod[2] ≈ +6127.36646120 atol=1e-7
    @test r_tod[3] ≈ +6380.34453270 atol=1e-7

    ## Quaternion
    ## ----------

    q_TOD_TEME = rECItoECI(Quaternion, TEME(), JD_UTC, TOD(), JD_UTC)

    r_tod = vect(conj(q_TOD_TEME)*r_teme*q_TOD_TEME)

    @test r_tod[1] ≈ +5094.51478040 atol=1e-7
    @test r_tod[2] ≈ +6127.36646120 atol=1e-7
    @test r_tod[3] ≈ +6380.34453270 atol=1e-7
end

################################################################################
#                           IAU-2006/2010 CIO-based
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
#
# one gets the following (this is the result in Table 3-6):
#
#   r_gcrf = 5102.50895290  i + 6123.01139910 j + 6378.13693380 k [km]
#
################################################################################

@testset "Function rECItoECI GCRF <=> CIRS" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## GCRF => CIRS
    ## ============

    r_gcrf = [5102.50895290; 6123.01139910; 6378.13693380]

    ## DCM
    ## ---

    D_CIRS_GCRF = rECItoECI(GCRF(), CIRS(), JD_UTC, eop_iau2000a)

    r_cirs = D_CIRS_GCRF*r_gcrf

    @test r_cirs[1] ≈ +5100.01840470 atol=1e-4
    @test r_cirs[2] ≈ +6122.78636480 atol=1e-4
    @test r_cirs[3] ≈ +6380.34453270 atol=1e-4

    ## Quaternion
    ## ----------

    q_CIRS_GCRF = rECItoECI(Quaternion, GCRF(), CIRS(), JD_UTC, eop_iau2000a)

    r_cirs = vect(conj(q_CIRS_GCRF)*r_gcrf*q_CIRS_GCRF)

    @test r_cirs[1] ≈ +5100.01840470 atol=1e-4
    @test r_cirs[2] ≈ +6122.78636480 atol=1e-4
    @test r_cirs[3] ≈ +6380.34453270 atol=1e-4

    ## CIRS => GCRF
    ## ============

    r_cirs = [+5100.01840470; +6122.78636480; +6380.34453270]

    ## DCM
    ## ---

    D_GCRF_CIRS = rECItoECI(CIRS(), GCRF(), JD_UTC, eop_iau2000a)

    r_gcrf = D_GCRF_CIRS*r_cirs

    @test r_gcrf[1] ≈ +5102.50895290 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01139910 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13693380 atol=1e-4

    ## Quaternion
    ## ----------

    q_GCRF_CIRS = rECItoECI(Quaternion, CIRS(), GCRF(), JD_UTC, eop_iau2000a)

    r_gcrf = vect(conj(q_GCRF_CIRS)*r_cirs*q_GCRF_CIRS)

    @test r_gcrf[1] ≈ +5102.50895290 atol=1e-4
    @test r_gcrf[2] ≈ +6123.01139910 atol=1e-4
    @test r_gcrf[3] ≈ +6378.13693380 atol=1e-4
end

################################################################################
#                         IAU-2006/2010 equinox-based
################################################################################

## GCRF <=> MJ2000
## ===============

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
#   UTC     = April 6, 2004, 07:51:28.386009
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#
# one gets the following (this is the result in Table 3-6):
#
#   r_gcrf = 5102.50895780  i + 6123.01140380 j + 6378.13692530 k [km]
#
# Notice that the transformation we are testing here does not convert to the
# original J2000. However, this result is close enough for a test comparison.
#
################################################################################

@testset "Function rECItoECI GCRF <=> MJ2000" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## MJ2000 => GCRF
    ## ==============

    r_mj2000 = [5102.50960000; 6123.01152000; 6378.13630000]

    ## DCM
    ## ---

    D_GCRF_MJ2000 = rECItoECI(MJ2000(), GCRF(), JD_UTC, eop_iau2000a)

    r_gcrf = D_GCRF_MJ2000*r_mj2000

    @test r_gcrf[1] ≈ +5102.50895780 atol=8e-4
    @test r_gcrf[2] ≈ +6123.01140380 atol=8e-4
    @test r_gcrf[3] ≈ +6378.13692530 atol=8e-4

    ## Quaternion
    ## ----------

    q_GCRF_MJ2000 = rECItoECI(Quaternion, MJ2000(), GCRF(), JD_UTC, eop_iau2000a)

    r_gcrf = vect(conj(q_GCRF_MJ2000)*r_mj2000*q_GCRF_MJ2000)

    @test r_gcrf[1] ≈ +5102.50895780 atol=8e-4
    @test r_gcrf[2] ≈ +6123.01140380 atol=8e-4
    @test r_gcrf[3] ≈ +6378.13692530 atol=8e-4

    ## GCRF => MJ2000
    ## ==============

    r_gcrf = [+5102.50895780; +6123.01140380; +6378.13692530]

    ## DCM
    ## ---

    D_MJ2000_GCRF = rECItoECI(GCRF(), MJ2000(), JD_UTC, eop_iau2000a)

    r_mj2000 = D_MJ2000_GCRF*r_gcrf

    @test r_mj2000[1] ≈ +5102.50960000 atol=8e-4
    @test r_mj2000[2] ≈ +6123.01152000 atol=8e-4
    @test r_mj2000[3] ≈ +6378.13630000 atol=8e-4

    ## Quaternion
    ## ----------

    q_MJ2000_GCRF = rECItoECI(Quaternion, GCRF(), MJ2000(), JD_UTC, eop_iau2000a)

    r_mj2000 = vect(conj(q_MJ2000_GCRF)*r_gcrf*q_MJ2000_GCRF)

    @test r_mj2000[1] ≈ +5102.50960000 atol=8e-4
    @test r_mj2000[2] ≈ +6123.01152000 atol=8e-4
    @test r_mj2000[3] ≈ +6378.13630000 atol=8e-4
end

## GCRF <=> MOD
## ============

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
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_gcrf = 5102.50895780  i + 6123.01140380 j + 6378.13692530 k [km]
#
# one gets the following (this is the result in Table 3-6):
#
#   r_mod  = +5094.02896110   i + 6127.87113500   j + 6380.24774200   k [km]
#
################################################################################

@testset "Function rECItoECI GCRF <=> MOD" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## MOD => GCRF
    ## ===========

    r_mod = [+5094.02896110; +6127.87113500; +6380.24774200]

    ## DCM
    ## ---

    D_GCRF_MOD = rECItoECI(MOD06(), GCRF(), JD_UTC, eop_iau2000a)

    r_gcrf = D_GCRF_MOD*r_mod

    @test r_gcrf[1] ≈ +5102.50895780 atol=1e-7
    @test r_gcrf[2] ≈ +6123.01140380 atol=1e-7
    @test r_gcrf[3] ≈ +6378.13692530 atol=1e-7

    ## Quaternion
    ## ----------

    q_GCRF_MOD = rECItoECI(Quaternion, MOD06(), GCRF(), JD_UTC, eop_iau2000a)

    r_gcrf = vect(conj(q_GCRF_MOD)*r_mod*q_GCRF_MOD)

    @test r_gcrf[1] ≈ +5102.50895780 atol=1e-7
    @test r_gcrf[2] ≈ +6123.01140380 atol=1e-7
    @test r_gcrf[3] ≈ +6378.13692530 atol=1e-7

    ## GCRF => MOD
    ## ===========

    r_gcrf = [+5102.50895780; +6123.01140380; +6378.13692530]

    ## DCM
    ## ---

    D_MOD_GCRF = rECItoECI(GCRF(), MOD06(), JD_UTC, eop_iau2000a)

    r_mod = D_MOD_GCRF*r_gcrf

    @test r_mod[1] ≈ +5094.02896110 atol=1e-7
    @test r_mod[2] ≈ +6127.87113500 atol=1e-7
    @test r_mod[3] ≈ +6380.24774200 atol=1e-7

    ## Quaternion
    ## ----------

    q_MOD_GCRF = rECItoECI(Quaternion, GCRF(), MOD06(), JD_UTC, eop_iau2000a)

    r_mod = vect(conj(q_MOD_GCRF)*r_gcrf*q_MOD_GCRF)

    @test r_mod[1] ≈ +5094.02896110 atol=1e-7
    @test r_mod[2] ≈ +6127.87113500 atol=1e-7
    @test r_mod[3] ≈ +6380.24774200 atol=1e-7
end

## GCRF <=> ERS
## ============

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
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_gcrf = 5102.50895780  i + 6123.01140380 j + 6378.13692530 k [km]
#
# one gets the following (this is the result in Table 3-6):
#
#   r_ers  = +5094.51462800   i + 6127.36658790   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECI GCRF <=> ERS" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## ERS => GCRF
    ## ===========

    r_ers = [+5094.51462800; +6127.36658790; +6380.34453270]

    ## DCM
    ## ---

    D_GCRF_ERS = rECItoECI(ERS(), GCRF(), JD_UTC, eop_iau2000a)

    r_gcrf = D_GCRF_ERS*r_ers

    @test r_gcrf[1] ≈ +5102.50895780 atol=5e-5
    @test r_gcrf[2] ≈ +6123.01140380 atol=5e-5
    @test r_gcrf[3] ≈ +6378.13692530 atol=5e-5

    ## Quaternion
    ## ----------

    q_GCRF_ERS = rECItoECI(Quaternion, ERS(), GCRF(), JD_UTC, eop_iau2000a)

    r_gcrf = vect(conj(q_GCRF_ERS)*r_ers*q_GCRF_ERS)

    @test r_gcrf[1] ≈ +5102.50895780 atol=5e-5
    @test r_gcrf[2] ≈ +6123.01140380 atol=5e-5
    @test r_gcrf[3] ≈ +6378.13692530 atol=5e-5

    ## GCRF => ERS
    ## ===========

    r_gcrf = [+5102.50895780; +6123.01140380; +6378.13692530]

    ## DCM
    ## ---

    D_ERS_GCRF = rECItoECI(GCRF(), ERS(), JD_UTC, eop_iau2000a)

    r_ers = D_ERS_GCRF*r_gcrf

    @test r_ers[1] ≈ +5094.51462800 atol=5e-5
    @test r_ers[2] ≈ +6127.36658790 atol=5e-5
    @test r_ers[3] ≈ +6380.34453270 atol=5e-5

    ## Quaternion
    ## ----------

    q_ERS_GCRF = rECItoECI(Quaternion, GCRF(), ERS(), JD_UTC, eop_iau2000a)

    r_ers = vect(conj(q_ERS_GCRF)*r_gcrf*q_ERS_GCRF)

    @test r_ers[1] ≈ +5094.51462800 atol=5e-5
    @test r_ers[2] ≈ +6127.36658790 atol=5e-5
    @test r_ers[3] ≈ +6380.34453270 atol=5e-5
end

## MJ2000 <=> MOD
## ==============

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
#   UTC     = April 6, 2004, 07:51:28.386009
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#
# one gets the following (this is the result in Table 3-6):
#
#   r_mod  = +5094.02896110   i + 6127.87113500   j + 6380.24774200   k [km]
#
# Notice that the transformation we are testing here does not convert to the
# original J2000. However, this result is close enough for a test comparison.
#
################################################################################

@testset "Function rECItoECI MJ2000 <=> MOD" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## MOD => MJ2000
    ## =============

    r_mod = [+5094.02896110; +6127.87113500; +6380.24774200]

    ## DCM
    ## ---

    D_MJ2000_MOD = rECItoECI(MOD06(), MJ2000(), JD_UTC, eop_iau2000a)

    r_mj2000 = D_MJ2000_MOD*r_mod

    @test r_mj2000[1] ≈ +5102.50960000 atol=8e-4
    @test r_mj2000[2] ≈ +6123.01152000 atol=8e-4
    @test r_mj2000[3] ≈ +6378.13630000 atol=8e-4

    ## Quaternion
    ## ----------

    q_MJ2000_MOD = rECItoECI(Quaternion, MOD06(), MJ2000(), JD_UTC, eop_iau2000a)

    r_mj2000 = vect(conj(q_MJ2000_MOD)*r_mod*q_MJ2000_MOD)

    @test r_mj2000[1] ≈ +5102.50960000 atol=8e-4
    @test r_mj2000[2] ≈ +6123.01152000 atol=8e-4
    @test r_mj2000[3] ≈ +6378.13630000 atol=8e-4

    ## MJ2000 => MOD
    ## =============

    r_mj2000 = [5102.50960000; 6123.01152000; 6378.13630000]

    ## DCM
    ## ---

    D_MOD_MJ2000 = rECItoECI(MJ2000(), MOD06(), JD_UTC, eop_iau2000a)

    r_mod = D_MOD_MJ2000*r_mj2000

    @test r_mod[1] ≈ +5094.02896110 atol=8e-4
    @test r_mod[2] ≈ +6127.87113500 atol=8e-4
    @test r_mod[3] ≈ +6380.24774200 atol=8e-4

    ## Quaternion
    ## ----------

    q_MOD_MJ2000 = rECItoECI(Quaternion, MJ2000(), MOD06(), JD_UTC, eop_iau2000a)

    r_mod = vect(conj(q_MOD_MJ2000)*r_mj2000*q_MOD_MJ2000)

    @test r_mod[1] ≈ +5094.02896110 atol=8e-4
    @test r_mod[2] ≈ +6127.87113500 atol=8e-4
    @test r_mod[3] ≈ +6380.24774200 atol=8e-4
end

## MJ2000 <=> ERS
## ==============

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
#   UTC     = April 6, 2004, 07:51:28.386009
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#
# one gets the following (this is the result in Table 3-6):
#
#   r_ers  = +5094.51462800   i + 6127.36658790   j + 6380.34453270   k [km]
#
# Notice that the transformation we are testing here does not convert to the
# original J2000. However, this result is close enough for a test comparison.
#
################################################################################

@testset "Function rECItoECI MJ2000 <=> ERS" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## ERS => MJ2000
    ## =============

    r_ers = [+5094.51462800; +6127.36658790; +6380.34453270]

    ## DCM
    ## ---

    D_MJ2000_ERS = rECItoECI(ERS(), MJ2000(), JD_UTC, eop_iau2000a)

    r_mj2000 = D_MJ2000_ERS*r_ers

    @test r_mj2000[1] ≈ +5102.50960000 atol=8e-4
    @test r_mj2000[2] ≈ +6123.01152000 atol=8e-4
    @test r_mj2000[3] ≈ +6378.13630000 atol=8e-4

    ## Quaternion
    ## ----------

    q_MJ2000_ERS = rECItoECI(Quaternion, ERS(), MJ2000(), JD_UTC, eop_iau2000a)

    r_mj2000 = vect(conj(q_MJ2000_ERS)*r_ers*q_MJ2000_ERS)

    @test r_mj2000[1] ≈ +5102.50960000 atol=8e-4
    @test r_mj2000[2] ≈ +6123.01152000 atol=8e-4
    @test r_mj2000[3] ≈ +6378.13630000 atol=8e-4

    ## MJ2000 => ERS
    ## =============

    r_mj2000 = [5102.50960000; 6123.01152000; 6378.13630000]

    ## DCM
    ## ---

    D_ERS_MJ2000 = rECItoECI(MJ2000(), ERS(), JD_UTC, eop_iau2000a)

    r_ers = D_ERS_MJ2000*r_mj2000

    @test r_ers[1] ≈ +5094.51462800 atol=8e-4
    @test r_ers[2] ≈ +6127.36658790 atol=8e-4
    @test r_ers[3] ≈ +6380.34453270 atol=8e-4

    ## Quaternion
    ## ----------

    q_ERS_MJ2000 = rECItoECI(Quaternion, MJ2000(), ERS(), JD_UTC, eop_iau2000a)

    r_ers = vect(conj(q_ERS_MJ2000)*r_mj2000*q_ERS_MJ2000)

    @test r_ers[1] ≈ +5094.51462800 atol=8e-4
    @test r_ers[2] ≈ +6127.36658790 atol=8e-4
    @test r_ers[3] ≈ +6380.34453270 atol=8e-4
end

## MOD <=> ERS
## ==============

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
#   UTC     = April 6, 2004, 07:51:28.386009
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#
# one gets the following (this is the result in Table 3-6):
#
#   r_ers  = +5094.51462800   i + 6127.36658790   j + 6380.34453270   k [km]
#
# Notice that the transformation we are testing here does not convert to the
# original J2000. However, this result is close enough for a test comparison.
#
################################################################################

@testset "Function rECItoECI MOD <=> ERS" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## ERS => MOD
    ## ==========

    r_ers = [+5094.51462800; +6127.36658790; +6380.34453270]

    ## DCM
    ## ---

    D_MOD_ERS = rECItoECI(ERS(), JD_UTC, MOD06(), JD_UTC, eop_iau2000a)

    r_mod = D_MOD_ERS*r_ers

    @test r_mod[1] ≈ +5094.02896110 atol=8e-4
    @test r_mod[2] ≈ +6127.87113500 atol=8e-4
    @test r_mod[3] ≈ +6380.24774200 atol=8e-4

    ## Quaternion
    ## ----------

    q_MOD_ERS = rECItoECI(Quaternion, ERS(), JD_UTC, MOD06(), JD_UTC, eop_iau2000a)

    r_mod = vect(conj(q_MOD_ERS)*r_ers*q_MOD_ERS)

    @test r_mod[1] ≈ +5094.02896110 atol=8e-4
    @test r_mod[2] ≈ +6127.87113500 atol=8e-4
    @test r_mod[3] ≈ +6380.24774200 atol=8e-4

    ## MOD => ERS
    ## ==========

    r_mod = [+5094.02896110; +6127.87113500; +6380.24774200]

    ## DCM
    ## ---

    D_ERS_MOD = rECItoECI(MOD06(), JD_UTC, ERS(), JD_UTC, eop_iau2000a)

    r_ers = D_ERS_MOD*r_mod

    @test r_ers[1] ≈ +5094.51462800 atol=8e-4
    @test r_ers[2] ≈ +6127.36658790 atol=8e-4
    @test r_ers[3] ≈ +6380.34453270 atol=8e-4

    ## Quaternion
    ## ----------

    q_ERS_MOD = rECItoECI(Quaternion, MOD06(), JD_UTC, ERS(), JD_UTC, eop_iau2000a)

    r_ers = vect(conj(q_ERS_MOD)*r_mod*q_ERS_MOD)

    @test r_ers[1] ≈ +5094.51462800 atol=8e-4
    @test r_ers[2] ≈ +6127.36658790 atol=8e-4
    @test r_ers[3] ≈ +6380.34453270 atol=8e-4
end
