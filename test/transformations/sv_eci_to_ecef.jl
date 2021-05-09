#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to ECI to ECEF transformations using satellite state vector.
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

# File: ./src/transformations/sv_eci_to_ecef.jl
# =============================================

# Functions: sv_eci_to_ecef
# -------------------------

# The rotations functions were already heavily tested in `ecef_to_eci.jl`,
# `fk5/fk5.jl`, and `iau2006/iau2006.jl`. Hence, here we will do only some minor
# testing involving the following transformations:
#
#   GCRF  => ITRF (FK5)
#   J2000 => PEF  (FK5, no corrections)
#   GCRF  => ITRF (IAU-2010)
#   GCRF  => TIRS (IAU-2010, no corrections)

################################################################################
#                                  IAU-76/FK5
################################################################################

# GCRF to ICRF
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
#   r_gcrf =  5102.50895790   i + 6123.01140070   j + 6378.13692820   k [km]
#   v_gcrf =    -4.7432201570 i +    0.7905364970 j +    5.533755270  k [km/s]
#
# one gets:
#
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#   v_itrf =    -3.225636520  i -    2.872451450  j +    5.531924446  k [km/s]
#
################################################################################

@testset "Function sv_eci_to_ecef GCRF => ITRF" begin
    JD_UTC = date_to_jd(2004, 4, 6, 7, 51, 28.386009)

    r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]
    v_gcrf = [-4.7432201570; 0.7905364970; 5.5337557270]
    sv_gcrf = orbsv(JD_UTC, r_gcrf, v_gcrf)
    sv_itrf = sv_eci_to_ecef(sv_gcrf, GCRF(), ITRF(), JD_UTC, eop_iau1980)

    @test sv_itrf.t === JD_UTC

    @test sv_itrf.r[1] ≈ -1033.4793830 atol=3e-4
    @test sv_itrf.r[2] ≈ +7901.2952754 atol=3e-4
    @test sv_itrf.r[3] ≈ +6380.3565958 atol=3e-4

    @test sv_itrf.v[1] ≈ -3.225636520  atol=8e-7
    @test sv_itrf.v[2] ≈ -2.872451450  atol=8e-7
    @test sv_itrf.v[3] ≈ +5.531924446  atol=8e-7
end

# J2000 to PEF
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
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#   v_j2000 =  -4.7432196000 i + 0.7905366000    j + 5.5337561900    k [km/s]
#
# one gets:
#
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#   v_pef  =    -3.2256327470 i -    2.8724425110 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Function sv_eci_to_ecef J2000 => PEF" begin
    JD_UT1 = date_to_jd(2004,4,6,7,51,28.386009) - 0.4399619/86400

    r_j2000  = [5102.50960000; 6123.01152000; 6378.13630000]
    v_j2000  = [-4.7432196000; 0.7905366000; 5.5337561900]
    sv_j2000 = orbsv(JD_UT1, r_j2000, v_j2000)
    sv_pef   = sv_eci_to_ecef(sv_j2000, J2000(), PEF(), JD_UT1)

    @test sv_pef.t === JD_UT1

    @test sv_pef.r[1] ≈ -1033.47503130 atol=1e-7
    @test sv_pef.r[2] ≈ +7901.30558560 atol=1e-7
    @test sv_pef.r[3] ≈ +6380.34453270 atol=1e-7

    @test sv_pef.v[1] ≈ -3.2256327470  atol=1e-7
    @test sv_pef.v[2] ≈ -2.8724425110  atol=1e-7
    @test sv_pef.v[3] ≈ +5.5319312880  atol=1e-7
end

################################################################################
#                                IAU-2006/2010
################################################################################

# GCRF to ITRF
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
#   r_gcrf = 5102.50895290   i + 6123.01139910   j + 6378.13693380   k [km]
#   v_gcrf =   -4.7432201610 i +    0.7905364950 j +    5.533755240  k [km/s]
#
# one gets:
#
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#   v_itrf =    -3.225636520  i -    2.872451450  j +    5.531924446  k [km/s]
#
################################################################################

@testset "Function sv_eci_to_ecef GCRF => ITRF" begin
    JD_UTC = date_to_jd(2004, 4, 6, 7, 51, 28.386009)

    r_gcrf = [5102.50895290; 6123.01139910; 6378.13693380]
    v_gcrf = [-4.7432201610; 0.7905364950; 5.5337557240]
    sv_gcrf = orbsv(JD_UTC, r_gcrf, v_gcrf)
    sv_itrf = sv_eci_to_ecef(sv_gcrf, GCRF(), ITRF(), JD_UTC, eop_iau2000a)

    @test sv_itrf.t === JD_UTC

    @test sv_itrf.r[1] ≈ -1033.4793830 atol=3e-4
    @test sv_itrf.r[2] ≈ +7901.2952754 atol=3e-4
    @test sv_itrf.r[3] ≈ +6380.3565958 atol=3e-4

    @test sv_itrf.v[1] ≈ -3.225636520  atol=8e-7
    @test sv_itrf.v[2] ≈ -2.872451450  atol=8e-7
    @test sv_itrf.v[3] ≈ +5.531924446  atol=8e-7
end

# GCRF to TIRS
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
#   r_gcrf = 5102.50895290   i + 6123.01139910   j + 6378.13693380   k [km]
#   v_gcrf =   -4.7432201610 i +    0.7905364950 j +    5.533755240  k [km/s]
#
# one gets:
#
#   r_tirs  = -1033.47503120   i + 7901.30558560   j + 6380.34453270   k [km]
#   v_tirs  =    -3.2256327470 i -    2.8724425110 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Function sv_eci_to_ecef TIRS => GCRF" begin
    JD_UT1 = date_to_jd(2004,4,6,7,51,28.386009) - 0.4399619/86400

    r_gcrf = [5102.50895290; 6123.01139910; 6378.13693380]
    v_gcrf = [-4.7432201610; 0.7905364950; 5.5337557240]
    sv_gcrf = orbsv(JD_UT1, r_gcrf, v_gcrf)
    sv_tirs = sv_eci_to_ecef(sv_gcrf, GCRF(), TIRS(), JD_UT1)

    @test sv_tirs.t === JD_UT1

    @test sv_tirs.r[1] ≈ -1033.47503120 atol=3e-4
    @test sv_tirs.r[2] ≈ +7901.30558560 atol=3e-4
    @test sv_tirs.r[3] ≈ +6380.34453270 atol=3e-4

    @test sv_tirs.v[1] ≈ -3.2256327470  atol=8e-7
    @test sv_tirs.v[2] ≈ -2.8724425110  atol=8e-7
    @test sv_tirs.v[3] ≈ +5.5319312880  atol=8e-7
end
