#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to the position of the Sun in equatorial and local frames.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/sun/sunpos.jl
# ===============================

# Function sun_position_el
# -----------------------

################################################################################
#                                 Test Results
################################################################################
#
# Finding the Sun position vector in equatorial and local frames
#
#   Using: JD = 2460065.945463, Latitude = 51.0 °N, Longitude = 10.0 °E
#           Pressure = 1.0 atm, Temperature = 20.0 °C
#
#   One must get:
#   RA = 02h 34m 09s; ~ 38.5375 °
#   DEC = +15° 06’ 42”; ~ 15.111667 °
#
################################################################################

@testset "Function sun_position_el" begin
    JD = 2460065.945463
    Latitude = 51.0
    Longitude = 10.0
    Pressure = 1.0
    Temperature = 20.0

    # Tolerances for various algorithms
    alg = ('1', '2', '3', '4', '5')
    tol = (1e-1; 1e-1; 1e-1; 1e-1; 1e-1)

    for algorithm in 1:5
        s_eq = sun_position_el(JD, Latitude, Longitude, Pressure, Temperature, 'e', alg[algorithm])

        @test s_eq[1] ≈ 38.5375 atol = tol[algorithm]
        @test s_eq[2] ≈ 15.111667 atol = tol[algorithm]
    end
end
