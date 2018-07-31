#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# Get the current EOP Data.
eop = read_iers_eop("./eop_IAU1980_2.txt", :IAU1980)

# File: ./src/transformations/orbit_elements.jl
# =============================================

# Function: change_oe_frame
# -------------------------

################################################################################
#                                 TEST RESULTS
################################################################################
#
# Scenario 01
# ===========
#
# Using AGI STK, it was verified that:
#
#   ╔══════════════════════╦═════════════════════════════╦═══════════════╗
#   ║   Reference Frame    ║            TEME             ║     J2000     ║
#   ╠══════════════════════╬═════════════════════════════╬═══════════════╣
#   ║ Epoch                ║ 1-Jun-2016 11:00:00.000 UTC ║ -             ║
#   ║ Mean motion          ║ 14.4 revs/day               ║ 14.4 revs/day ║
#   ║ Eccentricity         ║ 0.001111                    ║ 0.001111      ║
#   ║ Inclination          ║ 98.405˚                     ║ 98.3365˚      ║
#   ║ Argument of Perigee  ║ 90˚                         ║ 90.0604˚      ║
#   ║ RAAN                 ║ 227.336˚                    ║ 227.134˚      ║
#   ║ Mean Anomaly         ║ 320˚                        ║ 320˚          ║
#   ╚══════════════════════╩═════════════════════════════╩═══════════════╝
#
################################################################################

@testset "Function change_oe_frame" begin
    # This is the Amazonia-1 orbit, which has a = 7130.982 km.
    a = 7130.982e3
    e = 0.001111
    i = 98.405*pi/180
    Ω = 227.336*pi/180
    ω = 90*pi/180
    f = M_to_f(e, 320*pi/180)
    epoch = DatetoJD(2016,6,1,11,0,0)

    oe_TEME  = Orbit(0, a, e, i, Ω, ω, f)
    oe_J2000 = change_oe_frame(oe_TEME, TEME(), J2000(), epoch, eop)
    m = f_to_M(oe_J2000)

    @test oe_J2000.a/1000   ≈ 7130.982 atol=1e-8
    @test oe_J2000.e        ≈ e        atol=1e-8
    @test oe_J2000.i*180/pi ≈ 98.3365  atol=1e-4
    @test oe_J2000.Ω*180/pi ≈ 227.134  atol=1e-3
    @test oe_J2000.ω*180/pi ≈ 90.0604  atol=1e-4
    @test m*180/pi          ≈ 320      atol=1e-8
end

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
