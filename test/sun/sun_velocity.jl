#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to the velocity of the Sun.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/sun/velocity.jl
# ===========================

# Function sun_velocity_i
# -----------------------

################################################################################
#                                 Test Results
################################################################################
#
# The time-derivative of the Sun position vector computed by the function
# `sun_velocity_i` is compared to the numerical differentiation of the Sun
# position vector using the function `sun_position_i`, which was already
# validated.
#
################################################################################

@testset "Function sun_velocity_i" begin
    JD_start = DatetoJD(1950, 1, 1, 0, 0, 0)
    JD_stop  = DatetoJD(2019, 1, 1, 0, 0, 0)

    for trial = 1:100
        JD_UT1 = rand(JD_start:JD_stop)
        Δt     = 0.1
        s_t₁   = sun_position_i(JD_UT1)
        s_t₂   = sun_position_i(JD_UT1 + Δt/86400)
        v_n    = (s_t₂ - s_t₁)/Δt
        v      = sun_velocity_i(JD_UT1)

        @test norm(v-v_n)/norm(v)*100 < 0.055
    end
end
