#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to the position of the Moon.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#   [2] Meeus, J (1998). Astronomical algorithms. Willmann-Bell, Inc, Richmond,
#       VA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/sun/moon_position.jl
# ================================

# Function moon_position_i
# ------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Example 47.a: Calculate the geocentric longitude, latitude, and equatorial
#               horizontal parallax of the Moon for 1992 April 12, at 0h TD.
#
#   In the ecliptic plane w.r.t. the mean equator of date:
#
#       λ = 133°.162655 (longitude)
#       β = -3°.229126 (latitude)
#       Δ = 368409.7 km
#
################################################################################

@testset "Function moon_position_i (Meeus)" begin
    JD_TDB = DatetoJD(1992, 04, 12, 0, 0, 0)

    r_moon_mod = moon_position_i(JD_TDB)

    # We need to compute the mean ecliptic to obtain the variables in the
    # example.
    T_TDB = (JD_TDB - JD_J2000)/36525.0
    ϵ = @evalpoly(T_TDB, 23.439_291, -0.013_004_2, -1.64e-7, +5.04e-7)
    ϵ = mod2pi(deg2rad(ϵ))

    r_moon_ecl = create_rotation_matrix(ϵ, :X)*r_moon_mod

    λ = atand(r_moon_ecl[2], r_moon_ecl[1])
    β = atand(r_moon_ecl[3], √(r_moon_ecl[1]^2 + r_moon_ecl[2]^2))
    Δ = norm(r_moon_ecl)/1000

    @test λ ≈ 133.162655 atol = 1e-6
    @test β ≈ -3.229126  atol = 1e-6
    @test Δ ≈ 368409.7   atol = 1e-1
end

################################################################################
#                                 Test Results
################################################################################
#
# Example 5-3: Finding the Moon position vector [1, p. 289]
#
#   In April 28, 1994, 00:00 UTC, one gets:
#
#             | -134_240.626 |
#   r_J2000 = | -311_571.590 | km
#             | -126_693.785 |
#
#   Note: it assumes that UTC ≈ UT1 ≈ TDB.
#
################################################################################

@testset "Function moon_position_i (fast)" begin
    JD_UTC = DatetoJD(1994,4,28,0,0,0)
    JD_TBD = JD_UTC

    r_J2000 = moon_position_i(JD_TBD, Val(:fast))

    @test r_J2000[1]/1000 ≈ -134_240.626 atol = 1e-3
    @test r_J2000[2]/1000 ≈ -311_571.590 atol = 1e-3
    @test r_J2000[3]/1000 ≈ -126_693.785 atol = 1e-3
end
