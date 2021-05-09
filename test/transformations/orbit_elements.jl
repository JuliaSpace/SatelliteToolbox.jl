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
    epoch = date_to_jd(2016,6,1,11,0,0)

    oe_TEME  = KeplerianElements(0, a, e, i, Ω, ω, f)
    oe_J2000 = change_oe_frame(oe_TEME, TEME(), J2000(), epoch, eop)
    m = get_M(oe_J2000)

    @test oe_J2000.a/1000   ≈ 7130.982 atol=1e-8
    @test oe_J2000.e        ≈ e        atol=1e-8
    @test oe_J2000.i*180/pi ≈ 98.3365  atol=1e-4
    @test oe_J2000.Ω*180/pi ≈ 227.134  atol=1e-3
    @test oe_J2000.ω*180/pi ≈ 90.0604  atol=1e-4
    @test m*180/pi          ≈ 320      atol=1e-8
end

