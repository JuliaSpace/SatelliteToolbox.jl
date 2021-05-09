# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to IAU-2006 transformations.
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

# File: ./src/transformations/iau2006/iau2006-cio.jl
# ==================================================

# Functions r_itrf_to_tirs_iau2006 and r_tirs_to_itrf_iau2006
# -----------------------------------------------------------

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
#   x_p    = -0.140682"
#   y_p    = +0.333309"
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#   v_itrf =    -3.225636520  i -    2.872451450  j +    5.531924446  k [km/s]
#
# one gets the following:
#
#   r_tirs  = -1033.47503120   i + 7901.30558560   j + 6380.34453270   k [km]
#   v_tirs  =    -3.2256327470 i -    2.8724425110 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Function r_itrf_to_tirs_iau2006 and r_tirs_to_itrf_iau2006" begin
    JD_TT = 2453101.828154745
    x_p   = -0.140682*pi/(180*3600)
    y_p   = +0.333309*pi/(180*3600)

    ## r_itrf_to_tirs_iau2006
    ## ===================

    r_itrf = [-1033.4793830; 7901.2952754; 6380.3565958]
    v_itrf = [-3.225636520; -2.872451450; +5.531924446]

    ## DCM
    ## ---

    D_TIRS_ITRF = r_itrf_to_tirs_iau2006(JD_TT, x_p, y_p)

    r_tirs = D_TIRS_ITRF*r_itrf
    v_tirs = D_TIRS_ITRF*v_itrf

    @test r_tirs[1] ≈ -1033.47503120 atol=1e-7
    @test r_tirs[2] ≈ +7901.30558560 atol=1e-7
    @test r_tirs[3] ≈ +6380.34453270 atol=1e-7

    @test v_tirs[1] ≈ -3.2256327470  atol=1e-9
    @test v_tirs[2] ≈ -2.8724425110  atol=1e-9
    @test v_tirs[3] ≈ +5.5319312880  atol=1e-9

    ## Quaternion
    ## ----------

    q_TIRS_ITRF = r_itrf_to_tirs_iau2006(Quaternion, JD_TT, x_p, y_p)

    r_tirs = vect(q_TIRS_ITRF \ r_itrf * q_TIRS_ITRF)
    v_tirs = vect(q_TIRS_ITRF \ v_itrf * q_TIRS_ITRF)

    @test r_tirs[1] ≈ -1033.47503120 atol=1e-7
    @test r_tirs[2] ≈ +7901.30558560 atol=1e-7
    @test r_tirs[3] ≈ +6380.34453270 atol=1e-7

    @test v_tirs[1] ≈ -3.2256327470  atol=1e-9
    @test v_tirs[2] ≈ -2.8724425110  atol=1e-9
    @test v_tirs[3] ≈ +5.5319312880  atol=1e-9

    ## r_tirs_to_itrf_iau2006
    ## ==============

    r_tirs  = [-1033.47503120; 7901.30558560; 6380.34453270]
    v_tirs  = [-3.2256327470; -2.8724425110; +5.5319312880]

    ## DCM
    ## ---

    D_ITRF_TIRS = r_tirs_to_itrf_iau2006(JD_TT, x_p, y_p)
    r_itrf = D_ITRF_TIRS*r_tirs
    v_itrf = D_ITRF_TIRS*v_tirs

    @test r_itrf[1] ≈ -1033.4793830 atol=1e-7
    @test r_itrf[2] ≈ +7901.2952754 atol=1e-7
    @test r_itrf[3] ≈ +6380.3565958 atol=1e-7

    @test v_itrf[1] ≈ -3.225636520  atol=1e-9
    @test v_itrf[2] ≈ -2.872451450  atol=1e-9
    @test v_itrf[3] ≈ +5.531924446  atol=1e-9

    ## Quaternion
    ## ----------

    q_ITRF_TIRS = r_tirs_to_itrf_iau2006(Quaternion, JD_TT, x_p, y_p)
    r_itrf = vect(q_ITRF_TIRS\r_tirs*q_ITRF_TIRS)
    v_itrf = vect(q_ITRF_TIRS\v_tirs*q_ITRF_TIRS)

    @test r_itrf[1] ≈ -1033.4793830 atol=1e-7
    @test r_itrf[2] ≈ +7901.2952754 atol=1e-7
    @test r_itrf[3] ≈ +6380.3565958 atol=1e-7

    @test v_itrf[1] ≈ -3.225636520  atol=1e-9
    @test v_itrf[2] ≈ -2.872451450  atol=1e-9
    @test v_itrf[3] ≈ +5.531924446  atol=1e-9
end

# Functions r_tirs_to_cirs_iau2006 and r_cirs_to_tirs_iau2006
# -----------------------------------------------------------

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
#   LOD    = 0.0015563 s
#   r_tirs = -1033.47503120   i + 7901.30558560   j + 6380.34453270   k [km]
#   v_tirs =    -3.2256327470 i -    2.8724425110 j +    5.5319312880 k [km/s]
#
# one gets the following:
#
#   r_cirs = -5100.01840470   i + 6122.78636480   j + 6380.34453270   k [km]
#   v_cirs =    -4.7453803300 i -    0.7903414530 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Functions r_tirs_to_cirs_iau2006 and r_cirs_to_tirs_iau2006" begin
    JD_UT1 = 2453101.827406783
    LOD    = 0.0015563
    w      = 7.292115146706979e-5*(1-LOD/86400)

    # r_tirs_to_cirs_iau2006
    # ===================

    r_tirs  = [-1033.47503120; 7901.30558560; 6380.34453270]
    v_tirs  = [-3.2256327470; -2.8724425110; +5.5319312880]

    # DCM
    # ---

    D_CIRS_TIRS = r_tirs_to_cirs_iau2006(JD_UT1)

    r_cirs = D_CIRS_TIRS * r_tirs
    v_cirs = D_CIRS_TIRS * (v_tirs + [0, 0, w] × r_tirs)

    @test r_cirs[1] ≈ +5100.01840470 atol=1e-7
    @test r_cirs[2] ≈ +6122.78636480 atol=1e-7
    @test r_cirs[3] ≈ +6380.34453270 atol=1e-7

    @test v_cirs[1] ≈ -4.7453803300  atol=1e-9
    @test v_cirs[2] ≈ +0.7903414530  atol=1e-9
    @test v_cirs[3] ≈ +5.5319312880  atol=1e-9

    # Quaternion
    # ----------

    q_CIRS_TIRS = r_tirs_to_cirs_iau2006(Quaternion, JD_UT1)

    r_cirs = vect(q_CIRS_TIRS\r_tirs * q_CIRS_TIRS)
    v_cirs = vect(q_CIRS_TIRS\(v_tirs + [0, 0, w] × r_tirs) * q_CIRS_TIRS)

    @test r_cirs[1] ≈ +5100.01840470 atol=1e-7
    @test r_cirs[2] ≈ +6122.78636480 atol=1e-7
    @test r_cirs[3] ≈ +6380.34453270 atol=1e-7

    @test v_cirs[1] ≈ -4.7453803300  atol=1e-9
    @test v_cirs[2] ≈ +0.7903414530  atol=1e-9
    @test v_cirs[3] ≈ +5.5319312880  atol=1e-9

    # r_cirs_to_tirs_iau2006
    # ===================

    r_cirs = [+5100.01840470; +6122.78636480; +6380.34453270]
    v_cirs = [-4.7453803300; +0.7903414530; +5.5319312880]

    # DCM
    # ---

    D_TIRS_CIRS = r_cirs_to_tirs_iau2006(JD_UT1)

    r_tirs = D_TIRS_CIRS * r_cirs
    v_tirs = D_TIRS_CIRS * v_cirs - [0, 0, w] × r_tirs

    @test r_tirs[1] ≈ -1033.47503120 atol=1e-7
    @test r_tirs[2] ≈ +7901.30558560 atol=1e-7
    @test r_tirs[3] ≈ +6380.34453270 atol=1e-7

    @test v_tirs[1] ≈ -3.2256327470  atol=1e-9
    @test v_tirs[2] ≈ -2.8724425110  atol=1e-9
    @test v_tirs[3] ≈ +5.5319312880  atol=1e-9

    # Quaternion
    # ----------

    q_TIRS_CIRS = r_cirs_to_tirs_iau2006(Quaternion, JD_UT1)

    r_tirs = vect(q_TIRS_CIRS \ r_cirs * q_TIRS_CIRS)
    v_tirs = vect(q_TIRS_CIRS \ v_cirs * q_TIRS_CIRS) - [0, 0, w] × r_tirs

    @test r_tirs[1] ≈ -1033.47503120 atol=1e-7
    @test r_tirs[2] ≈ +7901.30558560 atol=1e-7
    @test r_tirs[3] ≈ +6380.34453270 atol=1e-7

    @test v_tirs[1] ≈ -3.2256327470  atol=1e-9
    @test v_tirs[2] ≈ -2.8724425110  atol=1e-9
    @test v_tirs[3] ≈ +5.5319312880  atol=1e-9
end
#
# Functions r_cirs_to_gcrf_iau2006 and r_gcrf_to_cirs_iau2006
# -----------------------------------------------------------

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
#   JD_UT1 = 2453101.827406783
#   dX     = -0.000205"
#   dY     = -0.000136"
#   r_cirs = -5100.01840470   i + 6122.78636480   j + 6380.34453270   k [km]
#   v_cirs =    -4.7453803300 i -    0.7903414530 j +    5.5319312880 k [km/s]
#
# one gets the following (this is the result in Table 3-6):
#
#   r_gcrf = 5102.50895290  i + 6123.01139910 j + 6378.13693380 k [km]
#   v_gcrf =  -4.7432201610 i + 0.7905364950  j + 5.5337557240  k [km/s]
#
################################################################################

@testset "Functions r_cirs_to_gcrf_iau2006 and r_gcrf_to_cirs_iau2006" begin
    JD_TT = 2453101.827406783
    dX    = -0.000205 * pi / 180 / 3600
    dY    = -0.000136 * pi / 180 / 3600

    # r_cirs_to_gcrf_iau2006
    # ===================

    r_cirs = [+5100.01840470; +6122.78636480; +6380.34453270]
    v_cirs = [-4.7453803300; +0.7903414530; +5.5319312880]

    # DCM
    # ---

    D_GCRF_CIRS = r_cirs_to_gcrf_iau2006(JD_TT, dX, dY)
    r_gcrf = D_GCRF_CIRS * r_cirs
    v_gcrf = D_GCRF_CIRS * v_cirs

    @test r_gcrf[1] ≈ +5102.50895290 atol=9e-7
    @test r_gcrf[2] ≈ +6123.01139910 atol=9e-7
    @test r_gcrf[3] ≈ +6378.13693380 atol=9e-7

    @test v_gcrf[1] ≈ -4.7432201610  atol=1e-9
    @test v_gcrf[2] ≈ +0.7905364950  atol=1e-9
    @test v_gcrf[3] ≈ +5.5337557240  atol=1e-9

    # Quaternion
    # ----------

    q_GCRF_CIRS = r_cirs_to_gcrf_iau2006(Quaternion, JD_TT, dX, dY)

    r_gcrf = vect(q_GCRF_CIRS \ r_cirs * q_GCRF_CIRS)
    v_gcrf = vect(q_GCRF_CIRS \ v_cirs * q_GCRF_CIRS)

    @test r_gcrf[1] ≈ +5102.50895290 atol=9e-7
    @test r_gcrf[2] ≈ +6123.01139910 atol=9e-7
    @test r_gcrf[3] ≈ +6378.13693380 atol=9e-7

    @test v_gcrf[1] ≈ -4.7432201610  atol=1e-9
    @test v_gcrf[2] ≈ +0.7905364950  atol=1e-9
    @test v_gcrf[3] ≈ +5.5337557240  atol=1e-9

    # Function r_gcrf_to_cirs_iau2006
    # ============================

    r_gcrf = [5102.50895290; 6123.01139910; 6378.13693380]
    v_gcrf = [-4.7432201610; 0.7905364950; 5.5337557240]

    # DCM
    # ---

    D_CIRS_GCRF = r_gcrf_to_cirs_iau2006(JD_TT, dX, dY)
    r_cirs = D_CIRS_GCRF * r_gcrf
    v_cirs = D_CIRS_GCRF * v_gcrf

    @test r_cirs[1] ≈ +5100.01840470 atol=9e-7
    @test r_cirs[2] ≈ +6122.78636480 atol=9e-7
    @test r_cirs[3] ≈ +6380.34453270 atol=9e-7

    @test v_cirs[1] ≈ -4.7453803300  atol=1e-9
    @test v_cirs[2] ≈ +0.7903414530  atol=1e-9
    @test v_cirs[3] ≈ +5.5319312880  atol=1e-9

    # Quaternion
    # ----------

    q_CIRS_GCRF = r_gcrf_to_cirs_iau2006(Quaternion, JD_TT, dX, dY)

    r_cirs = vect(q_CIRS_GCRF \ r_gcrf * q_CIRS_GCRF)
    v_cirs = vect(q_CIRS_GCRF \ v_gcrf * q_CIRS_GCRF)

    @test r_cirs[1] ≈ +5100.01840470 atol=9e-7
    @test r_cirs[2] ≈ +6122.78636480 atol=9e-7
    @test r_cirs[3] ≈ +6380.34453270 atol=9e-7

    @test v_cirs[1] ≈ -4.7453803300  atol=1e-9
    @test v_cirs[2] ≈ +0.7903414530  atol=1e-9
    @test v_cirs[3] ≈ +5.5319312880  atol=1e-9
end
