#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to the equation of time.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/sun/equation_of_time.jl
# ===================================

# Function equation_of_time
# -------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Figure 3-22: Equation of Time Variation [1, p. 179].
#
# From this figure, one can see that:
#
#    |     Day     | Equation of Time |
#    |-------------|------------------|
#    | February 11 | [-15,-14]        |
#    | May 11      | [+3,+4]          |
#    | July 26     | [-7,-6]          |
#    | November 2  | [16,17]          |
#
################################################################################

@testset "Function equation_of_time" begin
    # February 11.
    eot = equation_of_time(DatetoJD(2000, 2, 11, 0, 0, 0))*12*60/pi
    @test -15 < eot < -14

    # May 11.
    eot = equation_of_time(DatetoJD(2000, 5, 11, 0, 0, 0))*12*60/pi
    @test +3 < eot < +4

    # July 26.
    eot = equation_of_time(DatetoJD(2000, 7, 26, 0, 0, 0))*12*60/pi
    @test -7 < eot < -6

    # November 2.
    eot = equation_of_time(DatetoJD(2000, 11, 2, 0, 0, 0))*12*60/pi
    @test 16 < eot < 17
end
