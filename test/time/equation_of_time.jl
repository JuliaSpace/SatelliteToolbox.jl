## Description #############################################################################
#
# Tests related to the equation of time.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorne, CA.
#
############################################################################################

# == File: ./src/time/equation_of_time.jl ==================================================

############################################################################################
#                                       Test Results                                       #
############################################################################################
#
# Figure 3-22: Equation of Time Variation [1, p. 179].
#
# From this figure, one can see that:
#
#    |     Day     | Equation of Time |
#    |-------------|------------------|
#    | February 11 | [-15, -14]       |
#    | May 11      | [ +3,  +4]       |
#    | July 26     | [ -7,  -6]       |
#    | November 2  | [+16,  17]       |
#
############################################################################################

@testset "Function equation_of_time" begin
    ang_to_min = 12 / π * 60

    # February 11.
    eot = SatelliteToolbox.equation_of_time(date_to_jd(2000, 2, 11, 0, 0, 0)) * ang_to_min
    @test -15 < eot < -14

    # May 11.
    eot = SatelliteToolbox.equation_of_time(date_to_jd(2000, 5, 11, 0, 0, 0)) * ang_to_min
    @test +3 < eot < +4

    # July 26.
    eot = SatelliteToolbox.equation_of_time(date_to_jd(2000, 7, 26, 0, 0, 0)) * ang_to_min
    @test -7 < eot < -6

    # November 2.
    eot = SatelliteToolbox.equation_of_time(date_to_jd(2000, 11, 2, 0, 0, 0)) * ang_to_min
    @test 16 < eot < 17

    # The DateTime method must match the Julian Day method.
    eot_dt = SatelliteToolbox.equation_of_time(DateTime(2000, 11, 2))
    eot_jd = SatelliteToolbox.equation_of_time(date_to_jd(2000, 11, 2, 0, 0, 0))
    @test eot_dt == eot_jd
end
