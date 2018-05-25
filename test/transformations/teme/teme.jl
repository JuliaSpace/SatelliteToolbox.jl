#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to TEME transformations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-05-24: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/transformations/teme/teme.jl
# ========================================

# Functions rTEMEtoTOD and rTODtoTEME
# -----------------------------------

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
#   JD_TT = 2453101.828154745
#   δΔϵ   = -0.003875"
#   δΔψ   = -0.052195"
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#   v_teme =   -4.7461314870 i +    0.7858180410 j +    5.5319312880 k [km/s]
#
# one gets the following data:
#
#   r_tod = 5094.51620300   i + 6127.36527840   j + 6380.34453270   k [km]
#   v_tod =   -4.7460883850 i +    0.7860783240 j +    5.5319312880 k [km/s]
#
# Furthermore, using:
#
#   JD_TT = 2453101.828154745
#   δΔϵ   = 0"
#   δΔψ   = 0"
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#   v_teme =   -4.7461314870 i +    0.7858180410 j +    5.5319312880 k [km/s]
#
# one gets the following data:
#
#   r_tod = 5094.51478040   i + 6127.36646120   j + 6380.34453270   k [km]
#   v_tod =   -4.7460885670 i +    0.7860772220 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Functions rTEMEtoTOD and rTODtoTEME" begin
    JD_TT   = 2453101.828154745

    ## rTEMEtoTOD
    ## ==========

    ## -------------------------------------------------------------------------
    ##                             First Test
    ## -------------------------------------------------------------------------

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]
    v_teme = [-4.7461314870; 0.7858180410; 5.5319312880]

    ## DCM
    ## ---

    D_TOD_TEME = rTEMEtoTOD(JD_TT,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

    r_tod = D_TOD_TEME*r_teme
    v_tod = D_TOD_TEME*v_teme

    @test r_tod[1] ≈ +5094.51620300 atol=1e-7
    @test r_tod[2] ≈ +6127.36527840 atol=1e-7
    @test r_tod[3] ≈ +6380.34453270 atol=1e-7

    @test v_tod[1] ≈ -4.7460883850  atol=1e-9
    @test v_tod[2] ≈ +0.7860783240  atol=1e-9
    @test v_tod[3] ≈ +5.5319312880  atol=1e-9

    ## Quaternion
    ## ----------

    q_TOD_TEME = rTEMEtoTOD(Quaternion,
                            JD_TT,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

    r_tod = vect(conj(q_TOD_TEME)*r_teme*q_TOD_TEME)
    v_tod = vect(conj(q_TOD_TEME)*v_teme*q_TOD_TEME)

    @test r_tod[1] ≈ +5094.51620300 atol=1e-7
    @test r_tod[2] ≈ +6127.36527840 atol=1e-7
    @test r_tod[3] ≈ +6380.34453270 atol=1e-7

    @test v_tod[1] ≈ -4.7460883850  atol=1e-9
    @test v_tod[2] ≈ +0.7860783240  atol=1e-9
    @test v_tod[3] ≈ +5.5319312880  atol=1e-9

    ## -------------------------------------------------------------------------
    ##                            Second Test
    ## -------------------------------------------------------------------------

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]
    v_teme = [-4.7461314870; 0.7858180410; 5.5319312880]

    ## DCM
    ## ---

    D_TOD_TEME = rTEMEtoTOD(JD_TT)

    r_tod = D_TOD_TEME*r_teme
    v_tod = D_TOD_TEME*v_teme

    @test r_tod[1] ≈ +5094.51478040 atol=1e-7
    @test r_tod[2] ≈ +6127.36646120 atol=1e-7
    @test r_tod[3] ≈ +6380.34453270 atol=1e-7

    @test v_tod[1] ≈ -4.7460885670  atol=1e-9
    @test v_tod[2] ≈ +0.7860772220  atol=1e-9
    @test v_tod[3] ≈ +5.5319312880  atol=1e-9

    ## Quaternion
    ## ----------

    q_TOD_TEME = rTEMEtoTOD(Quaternion, JD_TT)

    r_tod = vect(conj(q_TOD_TEME)*r_teme*q_TOD_TEME)
    v_tod = vect(conj(q_TOD_TEME)*v_teme*q_TOD_TEME)

    @test r_tod[1] ≈ +5094.51478040 atol=1e-7
    @test r_tod[2] ≈ +6127.36646120 atol=1e-7
    @test r_tod[3] ≈ +6380.34453270 atol=1e-7

    @test v_tod[1] ≈ -4.7460885670  atol=1e-9
    @test v_tod[2] ≈ +0.7860772220  atol=1e-9
    @test v_tod[3] ≈ +5.5319312880  atol=1e-9

    ## rTODtoTEME
    ## ==========

    ## -------------------------------------------------------------------------
    ##                             First Test
    ## -------------------------------------------------------------------------

    r_tod = [5094.51620300; 6127.36527840; 6380.34453270]
    v_tod = [-4.7460883850; 0.7860783240; 5.5319312880]

    ## DCM
    ## ---

    D_TEME_TOD = rTODtoTEME(JD_TT,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

    r_teme = D_TEME_TOD*r_tod
    v_teme = D_TEME_TOD*v_tod

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9

    ## Quaternion
    ## ----------

    q_TEME_TOD = rTODtoTEME(Quaternion,
                            JD_TT,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

    r_teme = vect(conj(q_TEME_TOD)*r_tod*q_TEME_TOD)
    v_teme = vect(conj(q_TEME_TOD)*v_tod*q_TEME_TOD)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9

    ## -------------------------------------------------------------------------
    ##                            Second Test
    ## -------------------------------------------------------------------------

    r_tod = [5094.51478040; 6127.36646120; 6380.34453270]
    v_tod = [-4.7460885670; 0.7860772220; 5.5319312880]

    ## DCM
    ## ---

    D_TEME_TOD = rTODtoTEME(JD_TT)

    r_teme = D_TEME_TOD*r_tod
    v_teme = D_TEME_TOD*v_tod

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9

    ## Quaternion
    ## ----------

    q_TEME_TOD = rTODtoTEME(Quaternion, JD_TT)

    r_teme = vect(conj(q_TEME_TOD)*r_tod*q_TEME_TOD)
    v_teme = vect(conj(q_TEME_TOD)*v_tod*q_TEME_TOD)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9
end

# Functions rTEMEtoMOD and rMODtoTEME
# -----------------------------------

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
#   JD_TT = 2453101.828154745
#   δΔϵ   = -0.003875"
#   δΔψ   = -0.052195"
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#   v_teme =   -4.7461314870 i +    0.7858180410 j +    5.5319312880 k [km/s]
#
# one gets the following data:
#
#   r_mod = 5094.02837450   i + 6127.87081640   j + 6380.24851640   k [km]
#   v_mod =   -4.7462630520 i +    0.7860140450 j +    5.5317905620 k [km/s]
#
# Furthermore, using:
#
#   JD_TT = 2453101.828154745
#   δΔϵ   = 0"
#   δΔψ   = 0"
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#   v_teme =   -4.7461314870 i +    0.7858180410 j +    5.5319312880 k [km/s]
#
# one gets the following data:
#
#   r_mod = 5094.02901670   i + 6127.87093630   j + 6380.24788850   k [km]
#   v_mod =   -4.7462624950 i +    0.7860141490 j +    5.5317910250 k [km/s]
#
################################################################################

@testset "Functions rTEMEtoMOD and rMODtoTEME" begin
    JD_TT   = 2453101.828154745

    ## rTEMEtoMOD
    ## ==========

    ## -------------------------------------------------------------------------
    ##                             First Test
    ## -------------------------------------------------------------------------

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]
    v_teme = [-4.7461314870; 0.7858180410; 5.5319312880]

    ## DCM
    ## ---

    D_MOD_TEME = rTEMEtoMOD(JD_TT,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

    r_mod = D_MOD_TEME*r_teme
    v_mod = D_MOD_TEME*v_teme

    @test r_mod[1] ≈ +5094.02837450 atol=1e-7
    @test r_mod[2] ≈ +6127.87081640 atol=1e-7
    @test r_mod[3] ≈ +6380.24851640 atol=1e-7

    @test v_mod[1] ≈ -4.7462630520  atol=1e-9
    @test v_mod[2] ≈ +0.7860140450  atol=1e-9
    @test v_mod[3] ≈ +5.5317905620  atol=1e-9

    ## Quaternion
    ## ----------

    q_MOD_TEME = rTEMEtoMOD(Quaternion,
                            JD_TT,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

    r_mod = vect(conj(q_MOD_TEME)*r_teme*q_MOD_TEME)
    v_mod = vect(conj(q_MOD_TEME)*v_teme*q_MOD_TEME)

    @test r_mod[1] ≈ +5094.02837450 atol=1e-7
    @test r_mod[2] ≈ +6127.87081640 atol=1e-7
    @test r_mod[3] ≈ +6380.24851640 atol=1e-7

    @test v_mod[1] ≈ -4.7462630520  atol=1e-9
    @test v_mod[2] ≈ +0.7860140450  atol=1e-9
    @test v_mod[3] ≈ +5.5317905620  atol=1e-9

    ## -------------------------------------------------------------------------
    ##                            Second Test
    ## -------------------------------------------------------------------------

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]
    v_teme = [-4.7461314870; 0.7858180410; 5.5319312880]

    ## DCM
    ## ---

    D_MOD_TEME = rTEMEtoMOD(JD_TT)

    r_mod = D_MOD_TEME*r_teme
    v_mod = D_MOD_TEME*v_teme

    @test r_mod[1] ≈ +5094.02901670 atol=1e-7
    @test r_mod[2] ≈ +6127.87093630 atol=1e-7
    @test r_mod[3] ≈ +6380.24788850 atol=1e-7

    @test v_mod[1] ≈ -4.7462624950  atol=1e-9
    @test v_mod[2] ≈ +0.7860141490  atol=1e-9
    @test v_mod[3] ≈ +5.5317910250  atol=1e-9

    ## Quaternion
    ## ----------

    q_MOD_TEME = rTEMEtoMOD(Quaternion, JD_TT)

    r_mod = vect(conj(q_MOD_TEME)*r_teme*q_MOD_TEME)
    v_mod = vect(conj(q_MOD_TEME)*v_teme*q_MOD_TEME)

    @test r_mod[1] ≈ +5094.02901670 atol=1e-7
    @test r_mod[2] ≈ +6127.87093630 atol=1e-7
    @test r_mod[3] ≈ +6380.24788850 atol=1e-7

    @test v_mod[1] ≈ -4.7462624950  atol=1e-9
    @test v_mod[2] ≈ +0.7860141490  atol=1e-9
    @test v_mod[3] ≈ +5.5317910250  atol=1e-9

    ## rMODtoTEME
    ## ==========

    ## -------------------------------------------------------------------------
    ##                             First Test
    ## -------------------------------------------------------------------------

    r_mod = [5094.02837450; 6127.87081640; 6380.24851640]
    v_mod = [-4.7462630520; 0.7860140450; 5.5317905620]

    ## DCM
    ## ---

    D_TEME_MOD = rMODtoTEME(JD_TT,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

    r_teme = D_TEME_MOD*r_mod
    v_teme = D_TEME_MOD*v_mod

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9

    ## Quaternion
    ## ----------

    q_TEME_MOD = rMODtoTEME(Quaternion,
                            JD_TT,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

    r_teme = vect(conj(q_TEME_MOD)*r_mod*q_TEME_MOD)
    v_teme = vect(conj(q_TEME_MOD)*v_mod*q_TEME_MOD)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9

    ## -------------------------------------------------------------------------
    ##                            Second Test
    ## -------------------------------------------------------------------------

    r_mod = [5094.02901670; 6127.87093630; 6380.24788850]
    v_mod = [-4.7462624950; 0.7860141490; 5.5317910250]

    ## DCM
    ## ---

    D_TEME_MOD = rMODtoTEME(JD_TT)

    r_teme = D_TEME_MOD*r_mod
    v_teme = D_TEME_MOD*v_mod

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9

    ## Quaternion
    ## ----------

    q_TEME_MOD = rMODtoTEME(Quaternion, JD_TT)

    r_teme = vect(conj(q_TEME_MOD)*r_mod*q_TEME_MOD)
    v_teme = vect(conj(q_TEME_MOD)*v_mod*q_TEME_MOD)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9
end

# Functions rTEMEtoGCRF and rGCRFtoTEME
# -------------------------------------

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
#   JD_TT = 2453101.828154745
#   δΔϵ   = -0.003875"
#   δΔψ   = -0.052195"
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#   v_teme =   -4.7461314870 i +    0.7858180410 j +    5.5319312880 k [km/s]
#
# one gets the following data:
#
#   r_gcrf = 5102.50895790  i + 6123.01140070   j + 6378.13692820   k [km]
#   v_gcrf =  -4.7432201570 i + 0.7905364970    j + 5.5337557270    k [km/s]
#
# Furthermore, using:
#
#   JD_TT = 2453101.828154745
#   δΔϵ   = 0"
#   δΔψ   = 0"
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#   v_teme =   -4.7461314870 i +    0.7858180410 j +    5.5319312880 k [km/s]
#
# one gets the following data:
#
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#   v_j2000 =  -4.7432196000 i + 0.7905366000    j + 5.5337561900    k [km/s]
#
################################################################################

@testset "Functions rTEMEtoGCRF and rGCRFtoTEME" begin
    JD_TT   = 2453101.828154745

    ## rTEMEtoGCRF
    ## ===========

    ## -------------------------------------------------------------------------
    ##                             First Test
    ## -------------------------------------------------------------------------

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]
    v_teme = [-4.7461314870; 0.7858180410; 5.5319312880]

    ## DCM
    ## ---

    D_GCRF_TEME = rTEMEtoGCRF(JD_TT,
                              -0.003875*pi/(180*3600),
                              -0.052195*pi/(180*3600))

    r_gcrf = D_GCRF_TEME*r_teme
    v_gcrf = D_GCRF_TEME*v_teme

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-7
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-7
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-7

    @test v_gcrf[1] ≈ -4.7432201570  atol=1e-9
    @test v_gcrf[2] ≈ +0.7905364970  atol=1e-9
    @test v_gcrf[3] ≈ +5.5337557270  atol=1e-9

    ## Quaternion
    ## ----------

    q_GCRF_TEME = rTEMEtoGCRF(Quaternion,
                              JD_TT,
                              -0.003875*pi/(180*3600),
                              -0.052195*pi/(180*3600))

    r_gcrf = vect(conj(q_GCRF_TEME)*r_teme*q_GCRF_TEME)
    v_gcrf = vect(conj(q_GCRF_TEME)*v_teme*q_GCRF_TEME)

    @test r_gcrf[1] ≈ +5102.50895790 atol=1e-7
    @test r_gcrf[2] ≈ +6123.01140070 atol=1e-7
    @test r_gcrf[3] ≈ +6378.13692820 atol=1e-7

    @test v_gcrf[1] ≈ -4.7432201570  atol=1e-9
    @test v_gcrf[2] ≈ +0.7905364970  atol=1e-9
    @test v_gcrf[3] ≈ +5.5337557270  atol=1e-9

    ## -------------------------------------------------------------------------
    ##                            Second Test
    ## -------------------------------------------------------------------------

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]
    v_teme = [-4.7461314870; 0.7858180410; 5.5319312880]

    ## DCM
    ## ---

    D_J2000_TEME = rTEMEtoGCRF(JD_TT)

    r_j2000 = D_J2000_TEME*r_teme
    v_j2000 = D_J2000_TEME*v_teme

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-7
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-7
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-7

    @test v_j2000[1] ≈ -4.7432196000  atol=1e-9
    @test v_j2000[2] ≈ +0.7905366000  atol=1e-9
    @test v_j2000[3] ≈ +5.5337561900  atol=1e-9

    ## Quaternion
    ## ----------

    q_J2000_TEME = rTEMEtoGCRF(Quaternion, JD_TT)

    r_j2000 = vect(conj(q_J2000_TEME)*r_teme*q_J2000_TEME)
    v_j2000 = vect(conj(q_J2000_TEME)*v_teme*q_J2000_TEME)

    @test r_j2000[1] ≈ +5102.50960000 atol=1e-7
    @test r_j2000[2] ≈ +6123.01152000 atol=1e-7
    @test r_j2000[3] ≈ +6378.13630000 atol=1e-7

    @test v_j2000[1] ≈ -4.7432196000  atol=1e-9
    @test v_j2000[2] ≈ +0.7905366000  atol=1e-9
    @test v_j2000[3] ≈ +5.5337561900  atol=1e-9

    ## rGCRFtoTEME
    ## ===========

    ## -------------------------------------------------------------------------
    ##                             First Test
    ## -------------------------------------------------------------------------

    r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]
    v_gcrf = [-4.7432201570; 0.7905364970; 5.5337557270]

    ## DCM
    ## ---

    D_TEME_GCRF = rGCRFtoTEME(JD_TT,
                              -0.003875*pi/(180*3600),
                              -0.052195*pi/(180*3600))

    r_teme = D_TEME_GCRF*r_gcrf
    v_teme = D_TEME_GCRF*v_gcrf

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9

    ## Quaternion
    ## ----------

    q_TEME_GCRF = rGCRFtoTEME(Quaternion,
                              JD_TT,
                              -0.003875*pi/(180*3600),
                              -0.052195*pi/(180*3600))

    r_teme = vect(conj(q_TEME_GCRF)*r_gcrf*q_TEME_GCRF)
    v_teme = vect(conj(q_TEME_GCRF)*v_gcrf*q_TEME_GCRF)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9

    ## -------------------------------------------------------------------------
    ##                            Second Test
    ## -------------------------------------------------------------------------

    r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]
    v_j2000 = [-4.7432196000; 0.7905366000; 5.5337561900]

    ## DCM
    ## ---

    D_TEME_J2000 = rGCRFtoTEME(JD_TT)

    r_teme = D_TEME_J2000*r_j2000
    v_teme = D_TEME_J2000*v_j2000

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9

    ## Quaternion
    ## ----------

    q_TEME_J2000 = rGCRFtoTEME(Quaternion, JD_TT)

    r_teme = vect(conj(q_TEME_J2000)*r_j2000*q_TEME_J2000)
    v_teme = vect(conj(q_TEME_J2000)*v_j2000*q_TEME_J2000)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9
end

# Functions rTEMEtoPEF and rPEFtoTEME
# -----------------------------------

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
#   JD_UT1 = DatetoJD(2004,4,6,7,51,28.386009) - 0.4399619/86400
#   LOD    = 0.0015563 s
#   r_teme = 5094.18016210   i + 6127.64465950   j + 6380.34453270   k [km]
#   v_teme =   -4.7461314870 i +    0.7858180410 j +    5.5319312880 k [km/s]
#
# one gets the following data:
#
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#   v_pef  =    -3.2256327470 i -    2.8724425110 j +    5.5319312880 k [km/s]
#
################################################################################

@testset "Functions rTEMEtoPEF and rPEFtoTEME" begin
    JD_UT1 = DatetoJD(2004,4,6,7,51,28.386009) - 0.4399619/86400
    LOD    = 0.0015563
    w      = 7.292115146706979e-5*(1-LOD/86400)

    ## rTEMEtoPEF
    ## ==========

    r_teme = [5094.18016210; 6127.64465950; 6380.34453270]
    v_teme = [-4.7461314870; 0.7858180410; 5.5319312880]

    ## DCM
    ## ---

    D_TEME_PEF = rTEMEtoPEF(JD_UT1)

    r_pef = D_TEME_PEF*r_teme
    v_pef = D_TEME_PEF*v_teme - [0;0;w] × r_pef

    @test r_pef[1] ≈ -1033.47503130 atol=1e-7
    @test r_pef[2] ≈ +7901.30558560 atol=1e-7
    @test r_pef[3] ≈ +6380.34453270 atol=1e-7

    @test v_pef[1] ≈ -3.2256327470  atol=1e-9
    @test v_pef[2] ≈ -2.8724425110  atol=1e-9
    @test v_pef[3] ≈ +5.5319312880  atol=1e-9

    ## Quaternion
    ## ----------

    q_TEME_PEF = rTEMEtoPEF(Quaternion, JD_UT1)

    r_pef = vect(conj(q_TEME_PEF)*r_teme*q_TEME_PEF)
    v_pef = vect(conj(q_TEME_PEF)*v_teme*q_TEME_PEF) - [0;0;w] × r_pef

    @test r_pef[1] ≈ -1033.47503130 atol=1e-7
    @test r_pef[2] ≈ +7901.30558560 atol=1e-7
    @test r_pef[3] ≈ +6380.34453270 atol=1e-7

    @test v_pef[1] ≈ -3.2256327470  atol=1e-9
    @test v_pef[2] ≈ -2.8724425110  atol=1e-9
    @test v_pef[3] ≈ +5.5319312880  atol=1e-9

    ## rPEFtoTEME
    ## ==========

    r_pef  = [-1033.47503130; 7901.30558560; 6380.34453270]
    v_pef  = [-3.2256327470; -2.8724425110; +5.5319312880]

    ## DCM
    ## ---

    D_PEF_TEME = rPEFtoTEME(JD_UT1)

    r_teme = D_PEF_TEME*r_pef
    v_teme = D_PEF_TEME*(v_pef + [0;0;w] × r_pef)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9

    ## Quaternion
    ## ----------

    q_PEF_TEME = rPEFtoTEME(Quaternion, JD_UT1)

    r_teme = vect(conj(q_PEF_TEME)*r_pef*q_PEF_TEME)
    v_teme = vect(conj(q_PEF_TEME)*(v_pef + [0;0;w] × r_pef)*q_PEF_TEME)

    @test r_teme[1] ≈ +5094.18016210 atol=1e-7
    @test r_teme[2] ≈ +6127.64465950 atol=1e-7
    @test r_teme[3] ≈ +6380.34453270 atol=1e-7

    @test v_teme[1] ≈ -4.7461314870  atol=1e-9
    @test v_teme[2] ≈ +0.7858180410  atol=1e-9
    @test v_teme[3] ≈ +5.5319312880  atol=1e-9
end
