#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to ECI to ECEF transformations.
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

# File: ./src/transformations/eci_to_ecef.jl
# ==========================================

# Functions: rECItoECEF
# ---------------------

################################################################################
#                                  IAU-76/FK5
################################################################################

## GCRF to ITRF
## ============

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
#   r_gcrf = 5102.50895790  i + 6123.01140070   j + 6378.13692820   k [km]
#
# one gets:
#
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#
################################################################################

@testset "Function rECItoECEF GCRF => ITRF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]

    D_ITRF_GCRF = rECItoECEF(GCRF(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = D_ITRF_GCRF*r_gcrf

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    q_ITRF_GCRF = rECItoECEF(Quaternion, GCRF(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = vect(conj(q_ITRF_GCRF)*r_gcrf*q_ITRF_GCRF)

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    D_ITRF_GCRF = rECItoECEF(GCRF(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = D_ITRF_GCRF*r_gcrf

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    q_ITRF_GCRF = rECItoECEF(Quaternion, GCRF(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = vect(conj(q_ITRF_GCRF)*r_gcrf*q_ITRF_GCRF)

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4
end

## J2000 to ITRF
## =============

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
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#
################################################################################

@testset "Function rECItoECEF J2000 => ITRF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]

    D_ITRF_J2000 = rECItoECEF(J2000(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = D_ITRF_J2000*r_j2000

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    q_ITRF_J2000 = rECItoECEF(Quaternion, J2000(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = vect(conj(q_ITRF_J2000)*r_j2000*q_ITRF_J2000)

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    D_ITRF_J2000 = rECItoECEF(J2000(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = D_ITRF_J2000*r_j2000

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    q_ITRF_J2000 = rECItoECEF(Quaternion, J2000(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = vect(conj(q_ITRF_J2000)*r_j2000*q_ITRF_J2000)

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4
end

## TOD to ITRF
## ===========

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
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#
################################################################################

@testset "Function rECItoECEF TOD => ITRF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_tod = [5094.51620300; 6127.36527840; 6380.34453270]

    D_ITRF_TOD = rECItoECEF(TOD(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = D_ITRF_TOD*r_tod

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    q_ITRF_TOD = rECItoECEF(Quaternion, TOD(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = vect(conj(q_ITRF_TOD)*r_tod*q_ITRF_TOD)

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    D_ITRF_TOD = rECItoECEF(TOD(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = D_ITRF_TOD*r_tod

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    q_ITRF_TOD = rECItoECEF(Quaternion, TOD(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = vect(conj(q_ITRF_TOD)*r_tod*q_ITRF_TOD)

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4
end

## MOD to ITRF
## ===========

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
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#
################################################################################

@testset "Function rECItoECEF MOD => ITRF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_mod  = [5094.02837450; 6127.87081640; 6380.24851640]

    D_ITRF_MOD = rECItoECEF(MOD(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = D_ITRF_MOD*r_mod

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    q_ITRF_MOD = rECItoECEF(Quaternion, MOD(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = vect(conj(q_ITRF_MOD)*r_mod*q_ITRF_MOD)

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    D_ITRF_MOD = rECItoECEF(MOD(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = D_ITRF_MOD*r_mod

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    q_ITRF_MOD = rECItoECEF(Quaternion, MOD(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = vect(conj(q_ITRF_MOD)*r_mod*q_ITRF_MOD)

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4
end

## TEME to ITRF
## ============

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
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#
# one gets:
#
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#
################################################################################

@testset "Function rECItoECEF TEME => ITRF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]

    D_ITRF_TEME = rECItoECEF(TEME(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = D_ITRF_TEME*r_teme

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    q_ITRF_TEME = rECItoECEF(Quaternion, TEME(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = vect(conj(q_ITRF_TEME)*r_teme*q_ITRF_TEME)

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    D_ITRF_TEME = rECItoECEF(TEME(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = D_ITRF_TEME*r_teme

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    q_ITRF_TEME = rECItoECEF(Quaternion, TEME(), ITRF(), JD_UTC, eop_iau1980)
    r_itrf = vect(conj(q_ITRF_TEME)*r_teme*q_ITRF_TEME)

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4
end

## GCRF to PEF
## ===========

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
#   r_gcrf = 5102.50895790  i + 6123.01140070   j + 6378.13692820   k [km]
#
# one gets:
#
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECEF GCRF => PEF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]

    D_PEF_GCRF = rECItoECEF(GCRF(), PEF(), JD_UTC, eop_iau1980)
    r_pef = D_PEF_GCRF*r_gcrf

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    q_PEF_GCRF = rECItoECEF(Quaternion, GCRF(), PEF(), JD_UTC, eop_iau1980)
    r_pef = vect(conj(q_PEF_GCRF)*r_gcrf*q_PEF_GCRF)

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    D_PEF_GCRF = rECItoECEF(GCRF(), PEF(), JD_UTC, eop_iau1980)
    r_pef = D_PEF_GCRF*r_gcrf

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    q_PEF_GCRF = rECItoECEF(Quaternion, GCRF(), PEF(), JD_UTC, eop_iau1980)
    r_pef = vect(conj(q_PEF_GCRF)*r_gcrf*q_PEF_GCRF)

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4
end

## J2000 to PEF
## ============

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
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECEF J2000 => PEF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]

    D_PEF_J2000 = rECItoECEF(J2000(), PEF(), JD_UTC, eop_iau1980)
    r_pef = D_PEF_J2000*r_j2000

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    q_PEF_J2000 = rECItoECEF(Quaternion, J2000(), PEF(), JD_UTC, eop_iau1980)
    r_pef = vect(conj(q_PEF_J2000)*r_j2000*q_PEF_J2000)

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    D_PEF_J2000 = rECItoECEF(J2000(), PEF(), JD_UTC, eop_iau1980)
    r_pef = D_PEF_J2000*r_j2000

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    q_PEF_J2000 = rECItoECEF(Quaternion, J2000(), PEF(), JD_UTC, eop_iau1980)
    r_pef = vect(conj(q_PEF_J2000)*r_j2000*q_PEF_J2000)

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4
end

## TOD to PEF
## ==========

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
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECEF TOD => PEF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_tod  = [5094.51620300; 6127.36527840; 6380.34453270]

    D_PEF_TOD = rECItoECEF(TOD(), PEF(), JD_UTC, eop_iau1980)
    r_pef = D_PEF_TOD*r_tod

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    q_PEF_TOD = rECItoECEF(Quaternion, TOD(), PEF(), JD_UTC, eop_iau1980)
    r_pef = vect(conj(q_PEF_TOD)*r_tod*q_PEF_TOD)

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    D_PEF_TOD = rECItoECEF(TOD(), PEF(), JD_UTC, eop_iau1980)
    r_pef = D_PEF_TOD*r_tod

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    q_PEF_TOD = rECItoECEF(Quaternion, TOD(), PEF(), JD_UTC, eop_iau1980)
    r_pef = vect(conj(q_PEF_TOD)*r_tod*q_PEF_TOD)

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4
end

## PEF to MOD
## ==========

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
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#
# If not EOP correction is used, then we get the same result by using:
#
#   r_mod = 5094.02901670   i + 6127.87093630   j + 6380.24788850   k [km]
#
# Notice, however, that this result is obtained by using the correct UT1.
#
################################################################################

@testset "Function rECItoECEF MOD => PEF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_mod  = [5094.02837450; 6127.87081640; 6380.24851640]

    D_PEF_MOD = rECItoECEF(MOD(), PEF(), JD_UTC, eop_iau1980)
    r_pef = D_PEF_MOD*r_mod

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    q_PEF_MOD = rECItoECEF(Quaternion, MOD(), PEF(), JD_UTC, eop_iau1980)
    r_pef = vect(conj(q_PEF_MOD)*r_mod*q_PEF_MOD)

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    D_PEF_MOD = rECItoECEF(MOD(), PEF(), JD_UTC, eop_iau1980)
    r_pef = D_PEF_MOD*r_mod

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    q_PEF_MOD = rECItoECEF(Quaternion, MOD(), PEF(), JD_UTC, eop_iau1980)
    r_pef = vect(conj(q_PEF_MOD)*r_mod*q_PEF_MOD)

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    # No EOP corrections
    # ==================

    # Notice that if we use JD_UT1, then the correction to TT will be wrong.
    # However, this will lead to a much smaller error than assuming that UTC =
    # UT1.

    JD_UT1 = DatetoJD(2004,4,6,7,51,28.386009) - 0.4399619/86400
    r_mod  = [5094.02901670; 6127.87093630; 6380.24788850]

    D_PEF_MOD = rECItoECEF(MOD(), PEF(), JD_UT1)
    r_pef = D_PEF_MOD*r_mod

    @test r_pef[1] ≈ -1033.47503130 atol=8e-6
    @test r_pef[2] ≈ +7901.30558560 atol=8e-6
    @test r_pef[3] ≈ +6380.34453270 atol=8e-6

    q_PEF_MOD = rECItoECEF(Quaternion, MOD(), PEF(), JD_UT1)
    r_pef = vect(conj(q_PEF_MOD)*r_mod*q_PEF_MOD)

    @test r_pef[1] ≈ -1033.47503130 atol=8e-6
    @test r_pef[2] ≈ +7901.30558560 atol=8e-6
    @test r_pef[3] ≈ +6380.34453270 atol=8e-6
end

## TEME to PEF
## ===========

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
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#
# one gets:
#
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECEF TEME => PEF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]

    D_PEF_TEME = rECItoECEF(TEME(), PEF(), JD_UTC, eop_iau1980)
    r_pef = D_PEF_TEME*r_teme

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    q_PEF_TEME = rECItoECEF(Quaternion, TEME(), PEF(), JD_UTC, eop_iau1980)
    r_pef = vect(conj(q_PEF_TEME)*r_teme*q_PEF_TEME)

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    D_PEF_TEME = rECItoECEF(TEME(), PEF(), JD_UTC, eop_iau1980)
    r_pef = D_PEF_TEME*r_teme

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4

    q_PEF_TEME = rECItoECEF(Quaternion, TEME(), PEF(), JD_UTC, eop_iau1980)
    r_pef = vect(conj(q_PEF_TEME)*r_teme*q_PEF_TEME)

    @test r_pef[1] ≈ -1033.47503130 atol=3e-4
    @test r_pef[2] ≈ +7901.30558560 atol=3e-4
    @test r_pef[3] ≈ +6380.34453270 atol=3e-4
end

################################################################################
#                                IAU-2006/2010
################################################################################

## GCRF => TIRS
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
#   r_gcrf = 5102.50895290  i + 6123.01139910 j + 6378.13693380 k [km]
#
# one gets:
#
#   r_tirs = -1033.47503120   i + 7901.30558560   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECEF GCRF => TIRS" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_gcrf = [5102.50895290; 6123.01139910; 6378.13693380]

    D_TIRS_GCRF = rECItoECEF(GCRF(), TIRS(), JD_UTC, eop_iau2000a)
    r_tirs = D_TIRS_GCRF*r_gcrf

    @test r_tirs[1] ≈ -1033.47503120 atol=3e-4
    @test r_tirs[2] ≈ +7901.30558560 atol=3e-4
    @test r_tirs[3] ≈ +6380.34453270 atol=3e-4

    q_TIRS_GCRF = rECItoECEF(Quaternion, GCRF(), TIRS(), JD_UTC, eop_iau2000a)
    r_tirs = vect(q_TIRS_GCRF\r_gcrf*q_TIRS_GCRF)

    @test r_tirs[1] ≈ -1033.47503120 atol=3e-4
    @test r_tirs[2] ≈ +7901.30558560 atol=3e-4
    @test r_tirs[3] ≈ +6380.34453270 atol=3e-4
end

## GCRF => ITRF
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
#   r_gcrf = 5102.50895290  i + 6123.01139910 j + 6378.13693380 k [km]
#
# one gets:
#
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#
################################################################################

@testset "Function rECItoECEF GCRF => ITRF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_gcrf = [5102.50895290; 6123.01139910; 6378.13693380]

    D_ITRF_GCRF = rECItoECEF(GCRF(), ITRF(), JD_UTC, eop_iau2000a)
    r_itrf = D_ITRF_GCRF*r_gcrf

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    q_ITRF_GCRF = rECItoECEF(Quaternion, GCRF(), ITRF(), JD_UTC, eop_iau2000a)
    r_itrf = vect(q_ITRF_GCRF\r_gcrf*q_ITRF_GCRF)

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4
end

## CIRS => TIRS
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
#   r_cirs = -5100.01840470   i + 6122.78636480   j + 6380.34453270   k [km]
#
# one gets:
#
#   r_tirs = -1033.47503120   i + 7901.30558560   j + 6380.34453270   k [km]
#
################################################################################

@testset "Function rECItoECEF CIRS => TIRS" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_cirs = [+5100.01840470; +6122.78636480; +6380.34453270]

    D_TIRS_CIRS = rECItoECEF(CIRS(), TIRS(), JD_UTC, eop_iau2000a)
    r_tirs = D_TIRS_CIRS*r_cirs

    @test r_tirs[1] ≈ -1033.47503120 atol=3e-4
    @test r_tirs[2] ≈ +7901.30558560 atol=3e-4
    @test r_tirs[3] ≈ +6380.34453270 atol=3e-4

    q_TIRS_CIRS = rECItoECEF(Quaternion, CIRS(), TIRS(), JD_UTC, eop_iau2000a)
    r_tirs = vect(q_TIRS_CIRS\r_cirs*q_TIRS_CIRS)

    @test r_tirs[1] ≈ -1033.47503120 atol=3e-4
    @test r_tirs[2] ≈ +7901.30558560 atol=3e-4
    @test r_tirs[3] ≈ +6380.34453270 atol=3e-4
end

## CIRS => ITRF
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
#   r_cirs = -5100.01840470   i + 6122.78636480   j + 6380.34453270   k [km]
#
# one gets:
#
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#
################################################################################

@testset "Function rECItoECEF CIRS => ITRF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)
    r_cirs = [+5100.01840470; +6122.78636480; +6380.34453270]

    D_ITRF_CIRS = rECItoECEF(CIRS(), ITRF(), JD_UTC, eop_iau2000a)
    r_itrf = D_ITRF_CIRS*r_cirs

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4

    q_ITRF_CIRS = rECItoECEF(Quaternion, CIRS(), ITRF(), JD_UTC, eop_iau2000a)
    r_itrf = vect(q_ITRF_CIRS\r_cirs*q_ITRF_CIRS)

    @test r_itrf[1] ≈ -1033.4793830 atol=3e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=3e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=3e-4
end
