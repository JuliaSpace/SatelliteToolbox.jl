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

eop = read_iers_eop("./eop_IAU1980.txt", :IAU1980)

# File: ./src/transformations/ecef_to_ecef.jl
# ===========================================

# Functions: rECEFtoECEF
# ----------------------

# ITRF <=> PEF
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
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#   v_itrf =    -3.225636520  i -    2.872451450  j +    5.531924446  k [km/s]
#
# one gets:
#
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#   v_pef  =    -3.2256327470 i -    2.8724425110 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Function rECEFtoECEF ITRF <=> PEF" begin
    JD_UTC = DatetoJD(2004, 4, 6, 7, 51, 28.386009)

    ## ITRF => PEF
    ## ===========

    r_itrf = [-1033.4793830; 7901.2952754; 6380.3565958]
    v_itrf = [-3.225636520; -2.872451450; +5.531924446]

    ## DCM
    ## ---

    D_PEF_ITRF = rECEFtoECEF(ITRF(), PEF(), JD_UTC, eop)

    r_pef = D_PEF_ITRF*r_itrf
    v_pef = D_PEF_ITRF*v_itrf

    @test r_pef[1] ≈ -1033.47503130 atol=1e-4
    @test r_pef[2] ≈ +7901.30558560 atol=1e-4
    @test r_pef[3] ≈ +6380.34453270 atol=1e-4

    @test v_pef[1] ≈ -3.2256327470  atol=1e-7
    @test v_pef[2] ≈ -2.8724425110  atol=1e-7
    @test v_pef[3] ≈ +5.5319312880  atol=1e-7

    ## Quaternion
    ## ----------

    q_PEF_ITRF = rECEFtoECEF(Quaternion, ITRF(), PEF(), JD_UTC, eop)

    r_pef = vect(conj(q_PEF_ITRF)*r_itrf*q_PEF_ITRF)
    v_pef = vect(conj(q_PEF_ITRF)*v_itrf*q_PEF_ITRF)

    @test r_pef[1] ≈ -1033.47503130 atol=1e-4
    @test r_pef[2] ≈ +7901.30558560 atol=1e-4
    @test r_pef[3] ≈ +6380.34453270 atol=1e-4

    @test v_pef[1] ≈ -3.2256327470  atol=1e-7
    @test v_pef[2] ≈ -2.8724425110  atol=1e-7
    @test v_pef[3] ≈ +5.5319312880  atol=1e-7

    ## PEF => ITRF
    ## ===========

    r_pef  = [-1033.47503130; 7901.30558560; 6380.34453270]
    v_pef  = [-3.2256327470; -2.8724425110; +5.5319312880]

    ## DCM
    ## ---

    D_ITRF_PEF = rECEFtoECEF(PEF(), ITRF(), JD_UTC, eop)

    r_itrf = D_ITRF_PEF*r_pef
    v_itrf = D_ITRF_PEF*v_pef

    @test r_itrf[1] ≈ -1033.4793830 atol=1e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=1e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=1e-4

    @test v_itrf[1] ≈ -3.225636520  atol=1e-7
    @test v_itrf[2] ≈ -2.872451450  atol=1e-7
    @test v_itrf[3] ≈ +5.531924446  atol=1e-7

    ## Quaternion
    ## ----------

    q_ITRF_PEF = rECEFtoECEF(Quaternion, PEF(), ITRF(), JD_UTC, eop)

    r_itrf = vect(conj(q_ITRF_PEF)*r_pef*q_ITRF_PEF)
    v_itrf = vect(conj(q_ITRF_PEF)*v_pef*q_ITRF_PEF)

    @test r_itrf[1] ≈ -1033.4793830 atol=1e-4
    @test r_itrf[2] ≈ +7901.2952754 atol=1e-4
    @test r_itrf[3] ≈ +6380.3565958 atol=1e-4

    @test v_itrf[1] ≈ -3.225636520  atol=1e-7
    @test v_itrf[2] ≈ -2.872451450  atol=1e-7
    @test v_itrf[3] ≈ +5.531924446  atol=1e-7
end
