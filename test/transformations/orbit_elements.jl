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
#   Tests related orbit elements.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A., McClain, W. D (2013). Fundamentals of Astrodynamics
#       and Applications. Microcosm Press.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-05-26: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/transformations/orbit_elements.jl
# =============================================

# Functions: kepler_to_rv and rv_to_kepler
# ----------------------------------------

################################################################################
#                                 TEST RESULTS
################################################################################
#
# Scenario 01
# ===========
#
# Example 2-5: Finding position and velocity vectors (COE2RV Test Case) [1, p.
# 119-120].
#
# Cartesian representation:
#
#     r = 6525.344    I + 6861.535    J + 6449.125    K km
#     v =    4.902276 I +    5.533124 J -    1.975709 K km
#
# Orbit elements
#
#     ╔═════════════════╦══════════════╗
#     ║    Parameter    ║    Values    ║
#     ╠═════════════════╬══════════════╣
#     ║ p               ║ 11067.790 km ║
#     ║ Eccentricity    ║ 0.83285      ║
#     ║ Inclination     ║ 87.87°       ║
#     ║ RAAN            ║ 227.89°      ║
#     ║ Arg. of Perigee ║ 53.38°       ║
#     ║ True Anomaly    ║ 92.335°      ║
#     ╚═════════════════╩══════════════╝
#
################################################################################

@testset "Function kepler_to_rv and rv_to_kepler" begin

    ## kepler_to_rv
    ## ============

    p    = 11067.790*1000
    e    = 0.83285
    i    = 87.87*pi/180
    RAAN = 227.89*pi/180
    w    = 53.38*pi/180
    f    = 92.335*pi/180
    a    = p/(1-e^2)

    (r_i, v_i) = kepler_to_rv(a, e, i, RAAN, w, f)

    @test r_i[1]/1000 ≈ +6525.344 atol=5e-2
    @test r_i[2]/1000 ≈ +6861.535 atol=5e-2
    @test r_i[3]/1000 ≈ +6449.125 atol=5e-2

    @test v_i[1]/1000 ≈ +4.902276 atol=1e-4
    @test v_i[2]/1000 ≈ +5.533124 atol=1e-4
    @test v_i[3]/1000 ≈ -1.975709 atol=1e-4

    ## rv_to_kepler
    ## ============

    r_i = [6525.344; 6861.535; 6449.125]*1000
    v_i = [4.902276; 5.533124; -1.975709]*1000

    oe = rv_to_kepler(r_i..., v_i...)

    a, e, i, RAAN, w, f = oe.a, oe.e, oe.i, oe.Ω, oe.ω, oe.f
    p = a*(1-e^2)

    @test p/1000      ≈ 11067.790 atol=5e-2
    @test e           ≈ 0.83285   atol=1e-5
    @test i*180/pi    ≈ 87.87     atol=1e-2
    @test RAAN*180/pi ≈ 227.89    atol=1e-2
    @test w*180/pi    ≈ 53.38     atol=1e-2
    @test f*180/pi    ≈ 92.335    atol=1e-3
end
