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
#   Tests related to IAU-76/FK5 transformations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S (2006). Revisiting
#       Spacetrack Report #3: Rev1. AIAA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-04-18: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/transformations/fk5/nutation.jl
# ===========================================

println("$(c)Testing functions in file: ./src/transformations/fk5/nutation.jl$d")
println("$(c)----------------------------------------------------------------$d")
println("")

# Function nutation_fk5
# ---------------------

println("    Testing function nutation_fk5...")

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# In this example, using JD_TT = 2453101.828154745, one gets the following data
# related to the nutation:
#
#   mϵ_1980 =  23.4387368˚
#   Δϵ_1980 =   0.0020316˚
#   Δψ_1980 =  -0.0034108˚
#
# SatToolbox provides the following results:
#
#   julia> (mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(2453101.828154745)
#   (0.40908312815877673, 3.5458763448549555e-5, -5.953027070867465e-5)
#
#   julia> mϵ_1980*180/pi
#   23.438736713507268
#
#   julia> Δϵ_1980*180/pi
#   0.0020316374923546377
#
#   julia> Δψ_1980*180/pi
#   -0.0034108332648783257
#
################################################################################

(mϵ_1980, Δϵ_1980, Δψ_1980) = nutation_fk5(2453101.828154745)

@test mϵ_1980*180/pi ≈ 23.4387368 atol=1e-7
@test Δϵ_1980*180/pi ≈  0.0020316 atol=1e-7
@test Δψ_1980*180/pi ≈ -0.0034108 atol=1e-7

println("        $(b)Test passed!$d")
println("")

# File: ./src/transformations/fk5/precession.jl
# =============================================

println("$(c)Testing functions in file: ./src/transformations/fk5/precession.jl$d")
println("$(c)------------------------------------------------------------------$d")
println("")

# Function precession_fk5
# -----------------------

println("    Testing function precession_fk5...")

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction.
#
# In this example, using JD_TT = 2453101.828154745, one gets the following data
# related to the precession:
#
#   ζ = 0.0273055˚
#   Θ = 0.0237306˚
#   z = 0.0273059˚
#
# SatToolbox provides the following results:
#
#   julia> (ζ,Θ,z) = precession_fk5(2453101.828154745)
#   julia> ζ*180/pi
#   0.027305539219877804
#
#   julia> Θ*180/pi
#   0.023730619896279084
#
#   julia> z*180/pi
#   0.02730593931829398
#
################################################################################

(ζ,Θ,z) = precession_fk5(2453101.828154745)

@test ζ*180/pi ≈ 0.0273055 atol=1e-7
@test Θ*180/pi ≈ 0.0237306 atol=1e-7
@test z*180/pi ≈ 0.0273059 atol=1e-7

println("        $(b)Test passed!$d")
println("")

# File: ./src/transformations/fk5/fk5.jl
# ======================================

println("$(c)Testing functions in file: ./src/transformations/fk5/fk5.jl$d")
println("$(c)-----------------------------------------------------------$d")
println("")

# Functions rITRFtoPEF_fk5 and rPEFtoITRF_fk5
# -------------------------------------------

println("    Testing functions rITRFtoPEF_fk5 and rPEFtoITRF_fk5...")

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
#   x_p = -0.140682"
#   y_p = +0.333309"
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#   v_itrf =    -3.225636520  i -    2.872451450  j +    5.531924446  k [km/s]
#
# one gets the following data:
#
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#   v_pef  =    -3.2256327470 i -    2.8724425110 j +    5.5319312880 k [km/s]
#
################################################################################

x_p = -0.140682*pi/(180*3600)
y_p = +0.333309*pi/(180*3600)

## rITRFtoPEF_fk5
## ==============

r_itrf = [-1033.4793830; 7901.2952754; 6380.3565958]
v_itrf = [-3.225636520; -2.872451450; +5.531924446]

## DCM
## ---

D_PEF_ITRF = rITRFtoPEF_fk5(x_p, y_p)

r_pef = D_PEF_ITRF*r_itrf
v_pef = D_PEF_ITRF*v_itrf

@test r_pef[1] ≈ -1033.47503130 atol=1e-7
@test r_pef[2] ≈ +7901.30558560 atol=1e-7
@test r_pef[3] ≈ +6380.34453270 atol=1e-7

@test v_pef[1] ≈ -3.2256327470  atol=1e-9
@test v_pef[2] ≈ -2.8724425110  atol=1e-9
@test v_pef[3] ≈ +5.5319312880  atol=1e-9

## Quaternion
## ----------

q_PEF_ITRF = rITRFtoPEF_fk5(Quaternion, x_p, y_p)

r_pef = vect(conj(q_PEF_ITRF)*r_itrf*q_PEF_ITRF)
v_pef = vect(conj(q_PEF_ITRF)*v_itrf*q_PEF_ITRF)

@test r_pef[1] ≈ -1033.47503130 atol=1e-7
@test r_pef[2] ≈ +7901.30558560 atol=1e-7
@test r_pef[3] ≈ +6380.34453270 atol=1e-7

@test v_pef[1] ≈ -3.2256327470  atol=1e-9
@test v_pef[2] ≈ -2.8724425110  atol=1e-9
@test v_pef[3] ≈ +5.5319312880  atol=1e-9

## rPEFtoITRF_fk5
## ==============

r_pef  = [-1033.47503130; 7901.30558560; 6380.34453270]
v_pef  = [-3.2256327470; -2.8724425110; +5.5319312880]

## DCM
## ---

D_ITRF_PEF = rPEFtoITRF_fk5(x_p, y_p)
r_itrf = D_ITRF_PEF*r_pef
v_itrf = D_ITRF_PEF*v_pef

@test r_itrf[1] ≈ -1033.4793830 atol=1e-7
@test r_itrf[2] ≈ +7901.2952754 atol=1e-7
@test r_itrf[3] ≈ +6380.3565958 atol=1e-7

@test v_itrf[1] ≈ -3.225636520  atol=1e-9
@test v_itrf[2] ≈ -2.872451450  atol=1e-9
@test v_itrf[3] ≈ +5.531924446  atol=1e-9

## Quaternion
## ----------

q_ITRF_PEF = rPEFtoITRF_fk5(Quaternion, x_p, y_p)
r_itrf = vect(conj(q_ITRF_PEF)*r_pef*q_ITRF_PEF)
v_itrf = vect(conj(q_ITRF_PEF)*v_pef*q_ITRF_PEF)

@test r_itrf[1] ≈ -1033.4793830 atol=1e-7
@test r_itrf[2] ≈ +7901.2952754 atol=1e-7
@test r_itrf[3] ≈ +6380.3565958 atol=1e-7

@test v_itrf[1] ≈ -3.225636520  atol=1e-9
@test v_itrf[2] ≈ -2.872451450  atol=1e-9
@test v_itrf[3] ≈ +5.531924446  atol=1e-9

println("        $(b)Test passed!$d")
println("")

# Functions rPEFtoTOD_fk5 and rTODtoPEF_fk5
# -----------------------------------------

println("    Testing functions rPEFtoTOD_fk5 and rTODtoPEF_fk5...")

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
#   JD_TT  = 2453101.828154745
#   LOD    = 0.0015563 s
#   δΔψ    = -0.052195"
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#   v_pef  =    -3.2256327470 i -    2.8724425110 j +    5.5319312880 k [km/s]
#
# one gets the following data:
#
#   r_tod = 5094.51620300   i + 6127.36527840   j + 6380.34453270   k [km]
#   v_tod =   -4.7460883850 i +    0.7860783240 j +    5.5319312880 k [km/s]
#
# Furthermore, using:
#
#   JD_TT = 2453101.828154745
#   δΔψ   = 0"
#   r_pef = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#   v_pef =    -3.2256327470 i -    2.8724425110 j +    5.5319312880 k [km/s]
#
# one gets the following data:
#
#   r_tod = 5094.51478040   i + 6127.36646120   j + 6380.34453280   k [km]
#   v_tod =   -4.7460885670 i +    0.7860772220 j +    5.5319312880 k [km/s]
#
################################################################################

JD_UT1 = DatetoJD(2004,4,6,7,51,28.386009) - 0.4399619/86400
LOD    = 0.0015563
w      = 7.292115146706979e-5*(1-LOD/86400)

## rPEFtoTOD_fk5
## =============

## -----------------------------------------------------------------------------
##                                 First Test
## -----------------------------------------------------------------------------

r_pef  = [-1033.47503130; 7901.30558560; 6380.34453270]
v_pef  = [-3.2256327470; -2.8724425110; +5.5319312880]

## DCM
## ---

D_TOD_PEF = rPEFtoTOD_fk5(JD_UT1, 2453101.828154745, -0.052195*pi/(180*3600))

r_tod = D_TOD_PEF*r_pef
v_tod = D_TOD_PEF*(v_pef + [0;0;w] × r_pef)

@test r_tod[1] ≈ +5094.51620300 atol=1e-7
@test r_tod[2] ≈ +6127.36527840 atol=1e-7
@test r_tod[3] ≈ +6380.34453270 atol=1e-7

@test v_tod[1] ≈ -4.7460883850  atol=1e-9
@test v_tod[2] ≈ +0.7860783240  atol=1e-9
@test v_tod[3] ≈ +5.5319312880  atol=1e-9

## Quaternion
## ----------

q_TOD_PEF = rPEFtoTOD_fk5(Quaternion,
                          JD_UT1,
                          2453101.828154745,
                          -0.052195*pi/(180*3600))

r_tod = vect(conj(q_TOD_PEF)*r_pef*q_TOD_PEF)
v_tod = vect(conj(q_TOD_PEF)*(v_pef + [0;0;w] × r_pef)*q_TOD_PEF)

@test r_tod[1] ≈ +5094.51620300 atol=1e-7
@test r_tod[2] ≈ +6127.36527840 atol=1e-7
@test r_tod[3] ≈ +6380.34453270 atol=1e-7

@test v_tod[1] ≈ -4.7460883850  atol=1e-9
@test v_tod[2] ≈ +0.7860783240  atol=1e-9
@test v_tod[3] ≈ +5.5319312880  atol=1e-9

## -----------------------------------------------------------------------------
##                                Second Test
## -----------------------------------------------------------------------------

## DCM
## ---

D_TOD_PEF = rPEFtoTOD_fk5(JD_UT1, 2453101.828154745)

r_tod = D_TOD_PEF*r_pef
v_tod = D_TOD_PEF*(v_pef + [0;0;w] × r_pef)

@test r_tod[1] ≈ +5094.51478040 atol=1e-7
@test r_tod[2] ≈ +6127.36646120 atol=1e-7
@test r_tod[3] ≈ +6380.34453280 atol=1e-7

@test v_tod[1] ≈ -4.7460885670  atol=1e-9
@test v_tod[2] ≈ +0.7860772220  atol=1e-9
@test v_tod[3] ≈ +5.5319312880  atol=1e-9

## Quaternion
## ----------

q_TOD_PEF = rPEFtoTOD_fk5(Quaternion,
                          JD_UT1,
                          2453101.828154745)

r_tod = vect(conj(q_TOD_PEF)*r_pef*q_TOD_PEF)
v_tod = vect(conj(q_TOD_PEF)*(v_pef + [0;0;w] × r_pef)*q_TOD_PEF)

@test r_tod[1] ≈ +5094.51478040 atol=1e-7
@test r_tod[2] ≈ +6127.36646120 atol=1e-7
@test r_tod[3] ≈ +6380.34453280 atol=1e-7

@test v_tod[1] ≈ -4.7460885670  atol=1e-9
@test v_tod[2] ≈ +0.7860772220  atol=1e-9
@test v_tod[3] ≈ +5.5319312880  atol=1e-9

## rTODtoPEF_fk5
## =============

## -----------------------------------------------------------------------------
##                                 First Test
## -----------------------------------------------------------------------------

r_tod = [5094.51620300; 6127.36527840; 6380.34453270]
v_tod = [-4.7460883850; 0.7860783240; 5.5319312880]

## DCM
## ---

D_PEF_TOD = rTODtoPEF_fk5(JD_UT1, 2453101.828154745, -0.052195*pi/(180*3600))

r_pef = D_PEF_TOD*r_tod
v_pef = D_PEF_TOD*v_tod - [0;0;w] × r_pef

@test r_pef[1] ≈ -1033.47503130 atol=1e-7
@test r_pef[2] ≈ +7901.30558560 atol=1e-7
@test r_pef[3] ≈ +6380.34453270 atol=1e-7

@test v_pef[1] ≈ -3.2256327470  atol=1e-9
@test v_pef[2] ≈ -2.8724425110  atol=1e-9
@test v_pef[3] ≈ +5.5319312880  atol=1e-9

## Quaternion
## ----------

q_PEF_TOD = rTODtoPEF_fk5(Quaternion,
                          JD_UT1,
                          2453101.828154745,
                          -0.052195*pi/(180*3600))

r_pef = vect(conj(q_PEF_TOD)*r_tod*q_PEF_TOD)
v_pef = vect(conj(q_PEF_TOD)*v_tod*q_PEF_TOD) - [0;0;w] × r_pef

@test r_pef[1] ≈ -1033.47503130 atol=1e-7
@test r_pef[2] ≈ +7901.30558560 atol=1e-7
@test r_pef[3] ≈ +6380.34453270 atol=1e-7

@test v_pef[1] ≈ -3.2256327470  atol=1e-9
@test v_pef[2] ≈ -2.8724425110  atol=1e-9
@test v_pef[3] ≈ +5.5319312880  atol=1e-9

## -----------------------------------------------------------------------------
##                                Second Test
## -----------------------------------------------------------------------------

r_tod = [5094.51478040; 6127.36646120; 6380.34453280]
v_tod = [-4.7460885670; 0.7860772220; 5.5319312880]

## DCM
## ---

D_PEF_TOD = rTODtoPEF_fk5(JD_UT1, 2453101.828154745)

r_pef = D_PEF_TOD*r_tod
v_pef = D_PEF_TOD*v_tod - [0;0;w] × r_pef

@test r_pef[1] ≈ -1033.47503130 atol=1e-7
@test r_pef[2] ≈ +7901.30558560 atol=1e-7
@test r_pef[3] ≈ +6380.34453270 atol=1e-7

@test v_pef[1] ≈ -3.2256327470  atol=1e-9
@test v_pef[2] ≈ -2.8724425110  atol=1e-9
@test v_pef[3] ≈ +5.5319312880  atol=1e-9

## Quaternion
## ----------

q_PEF_TOD = rTODtoPEF_fk5(Quaternion,
                          JD_UT1,
                          2453101.828154745)

r_pef = vect(conj(q_PEF_TOD)*r_tod*q_PEF_TOD)
v_pef = vect(conj(q_PEF_TOD)*v_tod*q_PEF_TOD) - [0;0;w] × r_pef

@test r_pef[1] ≈ -1033.47503130 atol=1e-7
@test r_pef[2] ≈ +7901.30558560 atol=1e-7
@test r_pef[3] ≈ +6380.34453270 atol=1e-7

@test v_pef[1] ≈ -3.2256327470  atol=1e-9
@test v_pef[2] ≈ -2.8724425110  atol=1e-9
@test v_pef[3] ≈ +5.5319312880  atol=1e-9

println("        $(b)Test passed!$d")
println("")

# Functions rTODtoMOD_fk5 and rMODtoTOD_fk5
# -----------------------------------------

println("    Testing functions rTODtoMOD_fk5 and rMODtoTOD_fk5...")

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
#   r_tod = 5094.51620300   i + 6127.36527840   j + 6380.34453270   k [km]
#   v_tod =   -4.7460883850 i +    0.7860783240 j +    5.5319312880 k [km/s]
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
#   r_tod = 5094.51478040   i + 6127.36646120   j + 6380.34453280   k [km]
#   v_tod =   -4.7460885670 i +    0.7860772220 j +    5.5319312880 k [km/s]
#
# one gets the following data:
#
#   r_mod = 5094.02901670   i + 6127.87093630   j + 6380.24788850   k [km]
#   v_mod =   -4.7462624950 i +    0.7860141490 j +    5.5317910250 k [km/s]
#
################################################################################

## rTODtoMOD_fk5
## =============

## -----------------------------------------------------------------------------
##                                 First Test
## -----------------------------------------------------------------------------

r_tod = [5094.51620300; 6127.36527840; 6380.34453270]
v_tod = [-4.7460883850; 0.7860783240; 5.5319312880]

## DCM
## ---

D_MOD_TOD = rTODtoMOD_fk5(2453101.828154745,
                          -0.003875*pi/(180*3600),
                          -0.052195*pi/(180*3600))
r_mod = D_MOD_TOD*r_tod
v_mod = D_MOD_TOD*v_tod

@test r_mod[1] ≈ +5094.02837450 atol=1e-7
@test r_mod[2] ≈ +6127.87081640 atol=1e-7
@test r_mod[3] ≈ +6380.24851640 atol=1e-7

@test v_mod[1] ≈ -4.7462630520  atol=1e-9
@test v_mod[2] ≈ +0.7860140450  atol=1e-9
@test v_mod[3] ≈ +5.5317905620  atol=1e-9

## Quaternion
## ----------

q_MOD_TOD = rTODtoMOD_fk5(Quaternion,
                          2453101.828154745,
                          -0.003875*pi/(180*3600),
                          -0.052195*pi/(180*3600))

r_mod = vect(conj(q_MOD_TOD)*r_tod*q_MOD_TOD)
v_mod = vect(conj(q_MOD_TOD)*v_tod*q_MOD_TOD)

@test r_mod[1] ≈ +5094.02837450 atol=1e-7
@test r_mod[2] ≈ +6127.87081640 atol=1e-7
@test r_mod[3] ≈ +6380.24851640 atol=1e-7

@test v_mod[1] ≈ -4.7462630520  atol=1e-9
@test v_mod[2] ≈ +0.7860140450  atol=1e-9
@test v_mod[3] ≈ +5.5317905620  atol=1e-9

## -----------------------------------------------------------------------------
##                                Second Test
## -----------------------------------------------------------------------------

r_tod = [5094.51478040; 6127.36646120; 6380.34453280]
v_tod = [-4.7460885670; 0.7860772220; 5.5319312880]

## DCM
## ---

D_MOD_TOD = rTODtoMOD_fk5(2453101.82815474)

r_mod = D_MOD_TOD*r_tod
v_mod = D_MOD_TOD*v_tod

@test r_mod[1] ≈ +5094.02901670 atol=1e-7
@test r_mod[2] ≈ +6127.87093630 atol=1e-7
@test r_mod[3] ≈ +6380.24788850 atol=1e-7

@test v_mod[1] ≈ -4.7462624950  atol=1e-9
@test v_mod[2] ≈ +0.7860141490  atol=1e-9
@test v_mod[3] ≈ +5.5317910250  atol=1e-9

## Quaternion
## ----------

q_MOD_TOD = rTODtoMOD_fk5(Quaternion, 2453101.828154745)

r_mod = vect(conj(q_MOD_TOD)*r_tod*q_MOD_TOD)
v_mod = vect(conj(q_MOD_TOD)*v_tod*q_MOD_TOD)

@test r_mod[1] ≈ +5094.02901670 atol=1e-7
@test r_mod[2] ≈ +6127.87093630 atol=1e-7
@test r_mod[3] ≈ +6380.24788850 atol=1e-7

@test v_mod[1] ≈ -4.7462624950  atol=1e-9
@test v_mod[2] ≈ +0.7860141490  atol=1e-9
@test v_mod[3] ≈ +5.5317910250  atol=1e-9

## rMODtoTOD
## =========

## -----------------------------------------------------------------------------
##                                 First Test
## -----------------------------------------------------------------------------

r_mod = [5094.02837450; 6127.87081640; 6380.24851640]
v_mod = [-4.7462630520; 0.7860140450; 5.5317905620]

## DCM
## ---

D_TOD_MOD = rMODtoTOD_fk5(2453101.828154745,
                          -0.003875*pi/(180*3600),
                          -0.052195*pi/(180*3600))
r_tod = D_TOD_MOD*r_mod
v_tod = D_TOD_MOD*v_mod

@test r_tod[1] ≈ +5094.51620300 atol=1e-7
@test r_tod[2] ≈ +6127.36527840 atol=1e-7
@test r_tod[3] ≈ +6380.34453270 atol=1e-7

@test v_tod[1] ≈ -4.7460883850  atol=1e-9
@test v_tod[2] ≈ +0.7860783240  atol=1e-9
@test v_tod[3] ≈ +5.5319312880  atol=1e-9

## Quaternion
## ----------

q_TOD_MOD = rMODtoTOD_fk5(Quaternion,
                          2453101.828154745,
                          -0.003875*pi/(180*3600),
                          -0.052195*pi/(180*3600))

r_tod = vect(conj(q_TOD_MOD)*r_mod*q_TOD_MOD)
v_tod = vect(conj(q_TOD_MOD)*v_mod*q_TOD_MOD)

@test r_tod[1] ≈ +5094.51620300 atol=1e-7
@test r_tod[2] ≈ +6127.36527840 atol=1e-7
@test r_tod[3] ≈ +6380.34453270 atol=1e-7

@test v_tod[1] ≈ -4.7460883850  atol=1e-9
@test v_tod[2] ≈ +0.7860783240  atol=1e-9
@test v_tod[3] ≈ +5.5319312880  atol=1e-9

## -----------------------------------------------------------------------------
##                                Second Test
## -----------------------------------------------------------------------------

r_mod = [5094.02901670; 6127.87093630; 6380.24788850]
v_mod = [-4.7462624950; 0.7860141490; 5.5317910250]

## DCM
## ---

D_TOD_MOD = rMODtoTOD_fk5(2453101.828154745)

r_tod = D_TOD_MOD*r_mod
v_tod = D_TOD_MOD*v_mod

@test r_tod[1] ≈ +5094.51478040 atol=1e-7
@test r_tod[2] ≈ +6127.36646120 atol=1e-7
@test r_tod[3] ≈ +6380.34453280 atol=1e-7

@test v_tod[1] ≈ -4.7460885670  atol=1e-9
@test v_tod[2] ≈ +0.7860772220  atol=1e-9
@test v_tod[3] ≈ +5.5319312880  atol=1e-9

## Quaternion
## ----------

q_TOD_MOD = rMODtoTOD_fk5(Quaternion, 2453101.828154745)

r_tod = vect(conj(q_TOD_MOD)*r_mod*q_TOD_MOD)
v_tod = vect(conj(q_TOD_MOD)*v_mod*q_TOD_MOD)

@test r_tod[1] ≈ +5094.51478040 atol=1e-7
@test r_tod[2] ≈ +6127.36646120 atol=1e-7
@test r_tod[3] ≈ +6380.34453280 atol=1e-7

@test v_tod[1] ≈ -4.7460885670  atol=1e-9
@test v_tod[2] ≈ +0.7860772220  atol=1e-9
@test v_tod[3] ≈ +5.5319312880  atol=1e-9

println("        $(b)Test passed!$d")
println("")

# Functions rMODtoGCRF_fk5 and rGCRFtoMOD_fk5
# -------------------------------------------

println("    Testing functions rMODtoGCRF_fk5 and rGCRFtoMOD_fk5...")

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
#   r_mod = 5094.02837450   i + 6127.87081640   j + 6380.24851640   k [km]
#   v_mod =   -4.7462630520 i +    0.7860140450 j +    5.5317905620 k [km/s]
#
# one gets:
#
#   r_gcrf = 5102.50895790  i + 6123.01140070   j + 6378.13692820   k [km]
#   v_gcrf =  -4.7432201570 i + 0.7905364970    j + 5.5337557270    k [km/s]
#
# Furthermore, using:
#
#   JD_TT = 2453101.828154745
#   r_mod = 5094.02901670   i + 6127.87093630   j + 6380.24788850   k [km]
#   v_mod =   -4.7462624950 i +    0.7860141490 j +    5.5317910250 k [km/s]
#
# one gets:
#
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#   v_j2000 =  -4.7432196000 i + 0.7905366000    j + 5.5337561900    k [km/s]
#
################################################################################

## rMODtoGCRF
## ==========

## -----------------------------------------------------------------------------
##                                 First Test
## -----------------------------------------------------------------------------

r_mod = [5094.02837450; 6127.87081640; 6380.24851640]
v_mod = [-4.7462630520; 0.7860140450; 5.5317905620]

## DCM
## ---

D_GCRF_MOD = rMODtoGCRF_fk5(2453101.828154745)

r_gcrf = D_GCRF_MOD*r_mod
v_gcrf = D_GCRF_MOD*v_mod

@test r_gcrf[1] ≈ +5102.50895790 atol=1e-7
@test r_gcrf[2] ≈ +6123.01140070 atol=1e-7
@test r_gcrf[3] ≈ +6378.13692820 atol=1e-7

@test v_gcrf[1] ≈ -4.7432201570  atol=1e-9
@test v_gcrf[2] ≈ +0.7905364970  atol=1e-9
@test v_gcrf[3] ≈ +5.5337557270  atol=1e-9

## Quaternion
## ----------

q_GCRF_MOD = rMODtoGCRF_fk5(Quaternion, 2453101.828154745)

r_gcrf = vect(conj(q_GCRF_MOD)*r_mod*q_GCRF_MOD)
v_gcrf = vect(conj(q_GCRF_MOD)*v_mod*q_GCRF_MOD)

@test r_gcrf[1] ≈ +5102.50895790 atol=1e-7
@test r_gcrf[2] ≈ +6123.01140070 atol=1e-7
@test r_gcrf[3] ≈ +6378.13692820 atol=1e-7

@test v_gcrf[1] ≈ -4.7432201570  atol=1e-9
@test v_gcrf[2] ≈ +0.7905364970  atol=1e-9
@test v_gcrf[3] ≈ +5.5337557270  atol=1e-9

## -----------------------------------------------------------------------------
##                                Second Test
## -----------------------------------------------------------------------------

r_mod = [5094.02901670; 6127.87093630; 6380.24788850]
v_mod = [-4.7462624950; 0.7860141490; 5.5317910250]

## DCM
## ---

D_J2000_MOD = rMODtoGCRF_fk5(2453101.828154745)

r_j2000 = D_J2000_MOD*r_mod
v_j2000 = D_J2000_MOD*v_mod

@test r_j2000[1] ≈ +5102.50960000 atol=1e-7
@test r_j2000[2] ≈ +6123.01152000 atol=1e-7
@test r_j2000[3] ≈ +6378.13630000 atol=1e-7

@test v_j2000[1] ≈ -4.7432196000  atol=1e-9
@test v_j2000[2] ≈ +0.7905366000  atol=1e-9
@test v_j2000[3] ≈ +5.5337561900  atol=1e-9

## Quaternion
## ----------

q_J2000_MOD = rMODtoGCRF_fk5(Quaternion, 2453101.828154745)

r_j2000 = vect(conj(q_J2000_MOD)*r_mod*q_J2000_MOD)
v_j2000 = vect(conj(q_J2000_MOD)*v_mod*q_J2000_MOD)

@test r_j2000[1] ≈ +5102.50960000 atol=1e-7
@test r_j2000[2] ≈ +6123.01152000 atol=1e-7
@test r_j2000[3] ≈ +6378.13630000 atol=1e-7

@test v_j2000[1] ≈ -4.7432196000  atol=1e-9
@test v_j2000[2] ≈ +0.7905366000  atol=1e-9
@test v_j2000[3] ≈ +5.5337561900  atol=1e-9

## rGCRFtoMOD
## ==========

## -----------------------------------------------------------------------------
##                                 First Test
## -----------------------------------------------------------------------------

r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]
v_gcrf = [-4.7432201570; 0.7905364970; 5.5337557270]

## DCM
## ---

D_MOD_GCRF = rGCRFtoMOD_fk5(2453101.828154745)

r_mod = D_MOD_GCRF*r_gcrf
v_mod = D_MOD_GCRF*v_gcrf

@test r_mod[1] ≈ +5094.02837450 atol=1e-7
@test r_mod[2] ≈ +6127.87081640 atol=1e-7
@test r_mod[3] ≈ +6380.24851640 atol=1e-7

@test v_mod[1] ≈ -4.7462630520  atol=1e-9
@test v_mod[2] ≈ +0.7860140450  atol=1e-9
@test v_mod[3] ≈ +5.5317905620  atol=1e-9

## Quaternion
## ----------

q_MOD_GCRF = rGCRFtoMOD_fk5(Quaternion, 2453101.828154745)

r_mod = vect(conj(q_MOD_GCRF)*r_gcrf*q_MOD_GCRF)
v_mod = vect(conj(q_MOD_GCRF)*v_gcrf*q_MOD_GCRF)

@test r_mod[1] ≈ +5094.02837450 atol=1e-7
@test r_mod[2] ≈ +6127.87081640 atol=1e-7
@test r_mod[3] ≈ +6380.24851640 atol=1e-7

@test v_mod[1] ≈ -4.7462630520  atol=1e-9
@test v_mod[2] ≈ +0.7860140450  atol=1e-9
@test v_mod[3] ≈ +5.5317905620  atol=1e-9

## -----------------------------------------------------------------------------
##                                Second Test
## -----------------------------------------------------------------------------

r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]
v_j2000 = [-4.7432196000; 0.7905366000; 5.5337561900]

## DCM
## ---

D_MOD_J2000 = rGCRFtoMOD_fk5(2453101.828154745)

r_mod = D_MOD_J2000*r_j2000
v_mod = D_MOD_J2000*v_j2000

@test r_mod[1] ≈ +5094.02901670 atol=1e-7
@test r_mod[2] ≈ +6127.87093630 atol=1e-7
@test r_mod[3] ≈ +6380.24788850 atol=1e-7

@test v_mod[1] ≈ -4.7462624950  atol=1e-9
@test v_mod[2] ≈ +0.7860141490  atol=1e-9
@test v_mod[3] ≈ +5.5317910250  atol=1e-9

## Quaternion
## ----------

q_J2000_MOD = rGCRFtoMOD_fk5(Quaternion, 2453101.828154745)

r_mod = vect(conj(q_J2000_MOD)*r_j2000*q_J2000_MOD)
v_mod = vect(conj(q_J2000_MOD)*v_j2000*q_J2000_MOD)

@test r_mod[1] ≈ +5094.02901670 atol=1e-7
@test r_mod[2] ≈ +6127.87093630 atol=1e-7
@test r_mod[3] ≈ +6380.24788850 atol=1e-7

@test v_mod[1] ≈ -4.7462624950  atol=1e-9
@test v_mod[2] ≈ +0.7860141490  atol=1e-9
@test v_mod[3] ≈ +5.5317910250  atol=1e-9

println("        $(b)Test passed!$d")
println("")

# Functions rITRFtoGCRF_fk5 and rGCRFtoITRF_fk5
# ---------------------------------------------

println("    Testing functions rITRFtoGCRF_fk5 and rGCRFtoITRF_fk5...")

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
#   JD_TT  = 2453101.828154745
#   δΔϵ    = -0.003875"
#   δΔψ    = -0.052195"
#   x_p    = -0.140682"
#   y_p    = +0.333309"
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#
# one gets the following data:
#
#   r_gcrf = 5102.50895790  i + 6123.01140070   j + 6378.13692820   k [km]
#
# Furthermore, using:
#
#   JD_UT1 = DatetoJD(2004,4,6,7,51,28.386009) - 0.4399619/86400
#   JD_TT  = 2453101.828154745
#   x_p    = -0.140682"
#   y_p    = +0.333309"
#   r_itrf = -1033.4793830    i + 7901.2952754    j + 6380.3565958    k [km]
#
# one gets the following data:
#
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#
################################################################################

JD_UT1   = DatetoJD(2004,4,6,7,51,28.386009) - 0.4399619/86400
JD_TT    = 2453101.828154745
δΔϵ_1980 = -0.003875*pi/(180*3600)
δΔψ_1980 = -0.052195*pi/(180*3600)
x_p      = -0.140682*pi/(180*3600)
y_p      = +0.333309*pi/(180*3600)

## rITRFtoGCRF_fk5
## ===============

r_itrf = [-1033.4793830; 7901.2952754; 6380.3565958]

## -----------------------------------------------------------------------------
##                                 First Test
## -----------------------------------------------------------------------------

## DCM
## ---

D_GCRF_ITRF = rITRFtoGCRF_fk5(JD_UT1, JD_TT, x_p, y_p, δΔϵ_1980, δΔψ_1980)

r_gcrf = D_GCRF_ITRF*r_itrf

@test r_gcrf[1] ≈ +5102.50895790 atol=1e-7
@test r_gcrf[2] ≈ +6123.01140070 atol=1e-7
@test r_gcrf[3] ≈ +6378.13692820 atol=1e-7

## Quaternion
## ----------

q_GCRF_ITRF = rITRFtoGCRF_fk5(Quaternion,
                              JD_UT1,
                              JD_TT,
                              x_p,
                              y_p,
                              δΔϵ_1980,
                              δΔψ_1980)

r_gcrf = vect(conj(q_GCRF_ITRF)*r_itrf*q_GCRF_ITRF)

@test r_gcrf[1] ≈ +5102.50895790 atol=1e-7
@test r_gcrf[2] ≈ +6123.01140070 atol=1e-7
@test r_gcrf[3] ≈ +6378.13692820 atol=1e-7

## -----------------------------------------------------------------------------
##                                Second Test
## -----------------------------------------------------------------------------

## DCM
## ---

D_GCRF_ITRF = rITRFtoGCRF_fk5(JD_UT1, JD_TT, x_p, y_p)

r_j2000 = D_GCRF_ITRF*r_itrf

@test r_j2000[1] ≈ +5102.50960000 atol=1e-7
@test r_j2000[2] ≈ +6123.01152000 atol=1e-7
@test r_j2000[3] ≈ +6378.13630000 atol=1e-7

## Quaternion
## ----------

q_J2000_ITRF = rITRFtoGCRF_fk5(Quaternion, JD_UT1, JD_TT, x_p, y_p)

r_j2000 = vect(conj(q_J2000_ITRF)*r_itrf*q_J2000_ITRF)

@test r_j2000[1] ≈ +5102.50960000 atol=1e-7
@test r_j2000[2] ≈ +6123.01152000 atol=1e-7
@test r_j2000[3] ≈ +6378.13630000 atol=1e-7

# rGCRFtoITRF_fk5
# ===============

## -----------------------------------------------------------------------------
##                                 First Test
## -----------------------------------------------------------------------------

r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]

## DCM
## ---

D_ITRF_GCRF = rGCRFtoITRF_fk5(JD_UT1, JD_TT, x_p, y_p, δΔϵ_1980, δΔψ_1980)

r_itrf = D_ITRF_GCRF*r_gcrf

@test r_itrf[1] ≈ -1033.4793830 atol=1e-7
@test r_itrf[2] ≈ +7901.2952754 atol=1e-7
@test r_itrf[3] ≈ +6380.3565958 atol=1e-7

## Quaternion
## ----------

q_ITRF_GCRF = rGCRFtoITRF_fk5(Quaternion,
                              JD_UT1,
                              JD_TT,
                              x_p,
                              y_p,
                              δΔϵ_1980,
                              δΔψ_1980)

r_itrf = vect(conj(q_ITRF_GCRF)*r_gcrf*q_ITRF_GCRF)

@test r_itrf[1] ≈ -1033.4793830 atol=1e-7
@test r_itrf[2] ≈ +7901.2952754 atol=1e-7
@test r_itrf[3] ≈ +6380.3565958 atol=1e-7

## -----------------------------------------------------------------------------
##                                Second Test
## -----------------------------------------------------------------------------

r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]

## DCM
## ---

D_ITRF_J2000 = rGCRFtoITRF_fk5(JD_UT1, JD_TT, x_p, y_p)

r_itrf = D_ITRF_J2000*r_j2000

@test r_itrf[1] ≈ -1033.4793830 atol=1e-7
@test r_itrf[2] ≈ +7901.2952754 atol=1e-7
@test r_itrf[3] ≈ +6380.3565958 atol=1e-7

## Quaternion
## ----------

q_ITRF_J2000 = rGCRFtoITRF_fk5(Quaternion, JD_UT1, JD_TT, x_p, y_p)

r_itrf = vect(conj(q_ITRF_J2000)*r_j2000*q_ITRF_J2000)

@test r_itrf[1] ≈ -1033.4793830 atol=1e-7
@test r_itrf[2] ≈ +7901.2952754 atol=1e-7
@test r_itrf[3] ≈ +6380.3565958 atol=1e-7

println("        $(b)Test passed!$d")
println("")

# Functions rPEFtoGCRF_fk5 and rGCRFtoPEF_fk5
# -------------------------------------------

println("    Testing functions rPEFtoGCRF_fk5 and rGCRFtoPEF_fk5...")

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
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#
# one gets the following data:
#
#   r_gcrf = 5102.50895790  i + 6123.01140070   j + 6378.13692820   k [km]
#
# Furthermore, using:
#
#   JD_TT = 2453101.828154745
#   δΔϵ   = 0"
#   δΔψ   = 0"
#   r_pef  = -1033.47503130   i + 7901.30558560   j + 6380.34453270   k [km]
#
# one gets the following data:
#
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#
# NOTE: The velocity cannot be easily tested because it is necessary to account
# for the Earth rotation when converting from/to PEF.
#
################################################################################

JD_UT1 = DatetoJD(2004,4,6,7,51,28.386009) - 0.4399619/86400

## rPEFtoGCRF_fk5
## ==============

r_pef  = [-1033.47503130; 7901.30558560; 6380.34453270]

## -----------------------------------------------------------------------------
##                                 First Test
## -----------------------------------------------------------------------------

## DCM
## ---

D_GCRF_PEF = rPEFtoGCRF_fk5(JD_UT1,
                            2453101.828154745,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

r_gcrf = D_GCRF_PEF*r_pef

@test r_gcrf[1] ≈ +5102.50895790 atol=1e-7
@test r_gcrf[2] ≈ +6123.01140070 atol=1e-7
@test r_gcrf[3] ≈ +6378.13692820 atol=1e-7

## Quaternion
## ----------

q_GCRF_PEF = rPEFtoGCRF_fk5(Quaternion,
                            JD_UT1,
                            2453101.828154745,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

r_gcrf = vect(conj(q_GCRF_PEF)*r_pef*q_GCRF_PEF)

@test r_gcrf[1] ≈ +5102.50895790 atol=1e-7
@test r_gcrf[2] ≈ +6123.01140070 atol=1e-7
@test r_gcrf[3] ≈ +6378.13692820 atol=1e-7

## -----------------------------------------------------------------------------
##                                Second Test
## -----------------------------------------------------------------------------

## DCM
## ---

D_J2000_PEF = rPEFtoGCRF_fk5(JD_UT1, 2453101.828154745)

r_j2000 = D_J2000_PEF*r_pef

@test r_j2000[1] ≈ +5102.50960000 atol=1e-7
@test r_j2000[2] ≈ +6123.01152000 atol=1e-7
@test r_j2000[3] ≈ +6378.13630000 atol=1e-7

## Quaternion
## ----------

q_J2000_PEF = rPEFtoGCRF_fk5(Quaternion, JD_UT1, 2453101.828154745)

r_j2000 = vect(conj(q_J2000_PEF)*r_pef*q_J2000_PEF)

@test r_j2000[1] ≈ +5102.50960000 atol=1e-7
@test r_j2000[2] ≈ +6123.01152000 atol=1e-7
@test r_j2000[3] ≈ +6378.13630000 atol=1e-7

## rGCRFtoPEF_fk5
## ==============

## -----------------------------------------------------------------------------
##                                 First Test
## -----------------------------------------------------------------------------

r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]

## DCM
## ---

D_PEF_GCRF = rGCRFtoPEF_fk5(JD_UT1,
                            2453101.828154745,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

r_pef = D_PEF_GCRF*r_gcrf

@test r_pef[1] ≈ -1033.47503130 atol=1e-7
@test r_pef[2] ≈ +7901.30558560 atol=1e-7
@test r_pef[3] ≈ +6380.34453270 atol=1e-7

## Quaternion
## ----------

q_PEF_GCRF = rGCRFtoPEF_fk5(Quaternion,
                            JD_UT1,
                            2453101.828154745,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

r_pef = vect(conj(q_PEF_GCRF)*r_gcrf*q_PEF_GCRF)

@test r_pef[1] ≈ -1033.47503130 atol=1e-7
@test r_pef[2] ≈ +7901.30558560 atol=1e-7
@test r_pef[3] ≈ +6380.34453270 atol=1e-7

## -----------------------------------------------------------------------------
##                                Second Test
## -----------------------------------------------------------------------------

r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]

## DCM
## ---

D_PEF_J2000 = rGCRFtoPEF_fk5(JD_UT1, 2453101.828154745)

r_pef = D_PEF_J2000*r_j2000

@test r_pef[1] ≈ -1033.47503130 atol=1e-7
@test r_pef[2] ≈ +7901.30558560 atol=1e-7
@test r_pef[3] ≈ +6380.34453270 atol=1e-7

## Quaternion
## ----------

q_PEF_J2000 = rGCRFtoPEF_fk5(Quaternion, JD_UT1, 2453101.828154745)

r_pef = vect(conj(q_PEF_J2000)*r_j2000*q_PEF_J2000)

@test r_pef[1] ≈ -1033.47503130 atol=1e-7
@test r_pef[2] ≈ +7901.30558560 atol=1e-7
@test r_pef[3] ≈ +6380.34453270 atol=1e-7

println("        $(b)Test passed!$d")
println("")

# Functions rTODtoGCRF_fk5 and rGCRFtoTOD_fk5
# -------------------------------------------

println("    Testing functions rTODtoGCRF_fk5 and rGCRFtoTOD_fk5...")

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
#   r_tod = 5094.51620300   i + 6127.36527840   j + 6380.34453270   k [km]
#   v_tod =   -4.7460883850 i +    0.7860783240 j +    5.5319312880 k [km/s]
#
# one gets:
#
#   r_gcrf = 5102.50895790  i + 6123.01140070   j + 6378.13692820   k [km]
#   v_gcrf =  -4.7432201570 i + 0.7905364970    j + 5.5337557270    k [km/s]
#
# Furthermore, using:
#
#   JD_TT = 2453101.828154745
#   δΔϵ   = 0
#   δΔψ   = 0
#   r_tod = 5094.51478040   i + 6127.36646120   j + 6380.34453280   k [km]
#   v_tod =   -4.7460885670 i +    0.7860772220 j +    5.5319312880 k [km/s]
#
# one gets:
#
#   r_j2000 = 5102.50960000  i + 6123.01152000   j + 6378.13630000   k [km]
#   v_j2000 =  -4.7432196000 i + 0.7905366000    j + 5.5337561900    k [km/s]
#
################################################################################

## rTODtoGCRF
## ==========

## -----------------------------------------------------------------------------
##                                 First Test
## -----------------------------------------------------------------------------

r_tod = [5094.51620300; 6127.36527840; 6380.34453270]
v_tod = [-4.7460883850; 0.7860783240; 5.5319312880]

## DCM
## ---

D_GCRF_TOD = rTODtoGCRF_fk5(2453101.828154745,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))
r_gcrf = D_GCRF_TOD*r_tod
v_gcrf = D_GCRF_TOD*v_tod

@test r_gcrf[1] ≈ +5102.50895790 atol=1e-7
@test r_gcrf[2] ≈ +6123.01140070 atol=1e-7
@test r_gcrf[3] ≈ +6378.13692820 atol=1e-7

@test v_gcrf[1] ≈ -4.7432201570  atol=1e-9
@test v_gcrf[2] ≈ +0.7905364970  atol=1e-9
@test v_gcrf[3] ≈ +5.5337557270  atol=1e-9

## Quaternion
## ----------

q_GCRF_TOD = rTODtoGCRF_fk5(Quaternion,
                            2453101.828154745,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

r_gcrf = vect(conj(q_GCRF_TOD)*r_tod*q_GCRF_TOD)
v_gcrf = vect(conj(q_GCRF_TOD)*v_tod*q_GCRF_TOD)

@test r_gcrf[1] ≈ +5102.50895790 atol=1e-7
@test r_gcrf[2] ≈ +6123.01140070 atol=1e-7
@test r_gcrf[3] ≈ +6378.13692820 atol=1e-7

@test v_gcrf[1] ≈ -4.7432201570  atol=1e-9
@test v_gcrf[2] ≈ +0.7905364970  atol=1e-9
@test v_gcrf[3] ≈ +5.5337557270  atol=1e-9

## -----------------------------------------------------------------------------
##                                Second Test
## -----------------------------------------------------------------------------

r_tod = [5094.51478040; 6127.36646120; 6380.34453280]
v_tod = [-4.7460885670; 0.7860772220; 5.5319312880]

## DCM
## ---

D_J2000_TOD = rTODtoGCRF_fk5(2453101.828154745)

r_j2000 = D_J2000_TOD*r_tod
v_j2000 = D_J2000_TOD*v_tod

@test r_j2000[1] ≈ +5102.50960000 atol=1e-7
@test r_j2000[2] ≈ +6123.01152000 atol=1e-7
@test r_j2000[3] ≈ +6378.13630000 atol=1e-7

@test v_j2000[1] ≈ -4.7432196000  atol=1e-9
@test v_j2000[2] ≈ +0.7905366000  atol=1e-9
@test v_j2000[3] ≈ +5.5337561900  atol=1e-9

## Quaternion
## ----------

q_J2000_TOD = rTODtoGCRF_fk5(Quaternion, 2453101.828154745)

r_j2000 = vect(conj(q_J2000_TOD)*r_tod*q_J2000_TOD)
v_j2000 = vect(conj(q_J2000_TOD)*v_tod*q_J2000_TOD)

@test r_j2000[1] ≈ +5102.50960000 atol=1e-7
@test r_j2000[2] ≈ +6123.01152000 atol=1e-7
@test r_j2000[3] ≈ +6378.13630000 atol=1e-7

@test v_j2000[1] ≈ -4.7432196000  atol=1e-9
@test v_j2000[2] ≈ +0.7905366000  atol=1e-9
@test v_j2000[3] ≈ +5.5337561900  atol=1e-9

## rGCRFtoTOD
## ==========

## -----------------------------------------------------------------------------
##                                 First Test
## -----------------------------------------------------------------------------

r_gcrf = [5102.50895790; 6123.01140070; 6378.13692820]
v_gcrf = [-4.7432201570; 0.7905364970; 5.5337557270]

## DCM
## ---

D_TOD_GCRF = rGCRFtoTOD_fk5(2453101.828154745,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

r_tod = D_TOD_GCRF*r_gcrf
v_tod = D_TOD_GCRF*v_gcrf

@test r_tod[1] ≈ +5094.51620300 atol=1e-7
@test r_tod[2] ≈ +6127.36527840 atol=1e-7
@test r_tod[3] ≈ +6380.34453270 atol=1e-7

@test v_tod[1] ≈ -4.7460883850  atol=1e-9
@test v_tod[2] ≈ +0.7860783240  atol=1e-9
@test v_tod[3] ≈ +5.5319312880  atol=1e-9

## Quaternion
## ----------

q_TOD_GCRF = rGCRFtoTOD_fk5(Quaternion,
                            2453101.828154745,
                            -0.003875*pi/(180*3600),
                            -0.052195*pi/(180*3600))

r_tod = vect(conj(q_TOD_GCRF)*r_gcrf*q_TOD_GCRF)
v_tod = vect(conj(q_TOD_GCRF)*v_gcrf*q_TOD_GCRF)

@test r_tod[1] ≈ +5094.51620300 atol=1e-7
@test r_tod[2] ≈ +6127.36527840 atol=1e-7
@test r_tod[3] ≈ +6380.34453270 atol=1e-7

@test v_tod[1] ≈ -4.7460883850  atol=1e-9
@test v_tod[2] ≈ +0.7860783240  atol=1e-9
@test v_tod[3] ≈ +5.5319312880  atol=1e-9

## -----------------------------------------------------------------------------
##                                Second Test
## -----------------------------------------------------------------------------

r_j2000 = [5102.50960000; 6123.01152000; 6378.13630000]
v_j2000 = [-4.7432196000; 0.7905366000; 5.5337561900]

## DCM
## ---

D_TOD_J2000 = rGCRFtoTOD_fk5(2453101.828154745)

r_tod = D_TOD_J2000*r_j2000
v_tod = D_TOD_J2000*v_j2000

@test r_tod[1] ≈ +5094.51478040 atol=1e-7
@test r_tod[2] ≈ +6127.36646120 atol=1e-7
@test r_tod[3] ≈ +6380.34453280 atol=1e-7

@test v_tod[1] ≈ -4.7460885670  atol=1e-9
@test v_tod[2] ≈ +0.7860772220  atol=1e-9
@test v_tod[3] ≈ +5.5319312880  atol=1e-9

## Quaternion
## ----------

q_TOD_J2000 = rGCRFtoTOD_fk5(Quaternion, 2453101.828154745)

r_tod = vect(conj(q_TOD_J2000)*r_j2000*q_TOD_J2000)
v_tod = vect(conj(q_TOD_J2000)*v_j2000*q_TOD_J2000)

@test r_tod[1] ≈ +5094.51478040 atol=1e-7
@test r_tod[2] ≈ +6127.36646120 atol=1e-7
@test r_tod[3] ≈ +6380.34453280 atol=1e-7

@test v_tod[1] ≈ -4.7460885670  atol=1e-9
@test v_tod[2] ≈ +0.7860772220  atol=1e-9
@test v_tod[3] ≈ +5.5319312880  atol=1e-9

println("        $(b)Test passed!$d")
println("")

