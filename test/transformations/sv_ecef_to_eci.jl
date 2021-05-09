# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to ECEF to ECI transformations using satellite state vector.
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
# TODO: The EOP data obtained from IERS website does not match the values in the
# examples in [1]. However, this should be enough, because 1) the individual
# functions at the low level are tested using the same values of [1], and 2) the
# difference is smaller than 30 cm.

eop_iau1980  = read_iers_eop("./eop_IAU1980.txt",  :IAU1980)
eop_iau2000a = read_iers_eop("./eop_IAU2000A.txt", :IAU2000A)

# File: ./src/transformations/sv_ecef_to_eci.jl
# =============================================

# Functions: sv_ecef_to_eci
# -------------------------

# The rotations functions were already heavily tested in `ecef_to_eci.jl`,
# `fk5/fk5.jl`, and `iau2006/iau2006.jl`. Hence, here we will do only some minor
# testing involving the following transformations:
#
#   ITRF => GCRF  (FK5)
#   PEF  => J2000 (FK5, no corrections)
#   ITRF => GCRF  (IAU-2010)
#   TIRS => GCRF  (IAU-2010, no corrections)

################################################################################
#                                  IAU-76/FK5
################################################################################

# ITRF to GCRF
# ============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction [1, p. 230].
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#   v_itrf =    -3.225636520  i -    2.872451450  j +    5.531924446  k [km/s]
#
# one gets:
#
#   r_gcrf =  5102.50895790   i + 6123.01140070   j + 6378.13692820   k [km]
#   v_gcrf =    -4.7432201570 i +    0.7905364970 j +    5.533755270  k [km/s]
#
################################################################################

@testset "Function sv_ecef_to_eci ITRF => GCRF" begin
    JD_UTC = date_to_jd(2004, 4, 6, 7, 51, 28.386009)

    r_itrf  = [-1033.4793830; 7901.2952754; 6380.3565958]
    v_itrf  = [-3.225636520; -2.872451450; +5.531924446]
    sv_itrf = orbsv(JD_UTC, r_itrf, v_itrf)
    sv_gcrf = sv_ecef_to_eci(sv_itrf, ITRF(), GCRF(), JD_UTC, eop_iau1980)

    @test sv_gcrf.t === JD_UTC

    @test sv_gcrf.r[1] ≈ +5102.50895790 atol=3e-4
    @test sv_gcrf.r[2] ≈ +6123.01140070 atol=3e-4
    @test sv_gcrf.r[3] ≈ +6378.13692820 atol=3e-4

    @test sv_gcrf.v[1] ≈ -4.7432201570  atol=8e-7
    @test sv_gcrf.v[2] ≈ +0.7905364970  atol=8e-7
    @test sv_gcrf.v[3] ≈ +5.5337552700  atol=8e-7
end

# PEF to J2000
# ============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction [1, p. 230].
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   UT1    = UTC - 0.4399619/86400
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#   v_pef  =    -3.2256327470 i -    2.8724425110 j +    5.5319312880 k [km/s]
#
# one gets:
#
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#   v_j2000 =  -4.7432196000 i + 0.7905366000    j + 5.5337561900    k [km/s]
#
################################################################################

@testset "Function sv_ecef_to_eci PEF => J2000" begin
    JD_UT1 = date_to_jd(2004,4,6,7,51,28.386009) - 0.4399619/86400

    r_pef    = [-1033.47503130; 7901.30558560; 6380.34453270]
    v_pef    = [-3.2256327470; -2.8724425110; +5.5319312880]
    sv_pef   = orbsv(JD_UT1, r_pef, v_pef)
    sv_j2000 = sv_ecef_to_eci(sv_pef, PEF(), J2000(), JD_UT1)

    @test sv_j2000.t === JD_UT1

    @test sv_j2000.r[1] ≈ +5102.50960000 atol=1e-7
    @test sv_j2000.r[2] ≈ +6123.01152000 atol=1e-7
    @test sv_j2000.r[3] ≈ +6378.13630000 atol=1e-7

    @test sv_j2000.v[1] ≈ -4.7432196000  atol=1e-7
    @test sv_j2000.v[2] ≈ +0.7905366000  atol=1e-7
    @test sv_j2000.v[3] ≈ +5.5337561900  atol=1e-7
end

################################################################################
#                                IAU-2006/2010
################################################################################

# ITRF to GCRF
# ============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-2010 reduction [1, p. 220].
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#   v_itrf =    -3.225636520  i -    2.872451450  j +    5.531924446  k [km/s]
#
# one gets:
#
#   r_gcrf = 5102.50895290   i + 6123.01139910   j + 6378.13693380   k [km]
#   v_gcrf =   -4.7432201610 i +    0.7905364950 j +    5.533755240  k [km/s]
#
################################################################################

@testset "Function sv_ecef_to_eci ITRF => GCRF" begin
    JD_UTC = date_to_jd(2004, 4, 6, 7, 51, 28.386009)

    r_itrf  = [-1033.4793830; 7901.2952754; 6380.3565958]
    v_itrf  = [-3.225636520; -2.872451450; +5.531924446]
    sv_itrf = orbsv(JD_UTC, r_itrf, v_itrf)
    sv_gcrf = sv_ecef_to_eci(sv_itrf, ITRF(), GCRF(), JD_UTC, eop_iau2000a)

    @test sv_gcrf.t === JD_UTC

    @test sv_gcrf.r[1] ≈ +5102.50895290 atol=3e-4
    @test sv_gcrf.r[2] ≈ +6123.01139910 atol=3e-4
    @test sv_gcrf.r[3] ≈ +6378.13693380 atol=3e-4

    @test sv_gcrf.v[1] ≈ -4.7432201610  atol=8e-7
    @test sv_gcrf.v[2] ≈ +0.7905364950  atol=8e-7
    @test sv_gcrf.v[3] ≈ +5.5337557240  atol=8e-7
end

# TIRS to GCRF
# ============

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-2010 reduction [1, p. 220].
#
# According to this example and Table 3-6, using:
#
#   UTC    = April 6, 2004, 07:51:28.386009
#   UT1    = UTC - 0.4399619/86400
#   r_tirs  = -1033.47503120   i + 7901.30558560   j + 6380.34453270   k [km]
#   v_tirs  =    -3.2256327470 i -    2.8724425110 j +    5.5319312880 k [km/s]
#
# one gets:
#
#   r_gcrf = 5102.50895290   i + 6123.01139910   j + 6378.13693380   k [km]
#   v_gcrf =   -4.7432201610 i +    0.7905364950 j +    5.533755240  k [km/s]
#
################################################################################

@testset "Function sv_ecef_to_eci TIRS => GCRF" begin
    JD_UT1 = date_to_jd(2004, 4, 6, 7, 51, 28.386009) - 0.4399619/86400

    r_tirs  = [-1033.47503120; 7901.30558560; 6380.34453270]
    v_tirs  = [-3.2256327470; -2.8724425110; +5.5319312880]
    sv_tirs = orbsv(JD_UT1, r_tirs, v_tirs)
    sv_gcrf = sv_ecef_to_eci(sv_tirs, TIRS(), GCRF(), JD_UT1)

    @test sv_gcrf.t === JD_UT1

    @test sv_gcrf.r[1] ≈ +5102.50895290 atol=3e-4
    @test sv_gcrf.r[2] ≈ +6123.01139910 atol=3e-4
    @test sv_gcrf.r[3] ≈ +6378.13693380 atol=3e-4

    @test sv_gcrf.v[1] ≈ -4.7432201610  atol=8e-7
    @test sv_gcrf.v[2] ≈ +0.7905364950  atol=8e-7
    @test sv_gcrf.v[3] ≈ +5.5337557240  atol=8e-7
end
