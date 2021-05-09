#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to the Jacchia-Robert 1971 model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/earth/atmospheric_models/jr1971
# ===========================================

# Functions: jr1971
# -----------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
#   Values obtained from GMAT R2018a using the following inputs:
#
#       Date:      2017-01-01 00:00:00 UTC
#       Latitude:  45 deg
#       Longitude: 0 deg
#       F10.7:     100
#       F10.7ₐ:    100
#       Kp:        4
#
#   Result:
#
#       | Altitude [km]  | Density [g/cm³] |
#       |----------------|------------------|
#       |             92 | 3.87506e-09      |
#       |            100 | 3.96585e-09      |
#       |          100.1 | 6.8354e-10       |
#       |          110.5 | 1.2124e-10       |
#       |            125 | 1.60849e-11      |
#       |          125.1 | 1.58997e-11      |
#       |            300 | 1.30609e-14      |
#       |            700 | 1.34785e-17      |
#       |           1500 | 4.00464e-19      |
#
# Scenario 02
# ===========
#
#   Values obtained from GMAT R2018a using the following inputs:
#
#       Date:      2017-01-01 00:00:00 UTC
#       Latitude:  45 deg
#       Longitude: 0 deg
#       F10.7:     100
#       F10.7ₐ:    100
#       Kp:        1
#
#   Result:
#
#       | Altitude [km]  | Density [g/cm³] |
#       |----------------|------------------|
#       |             92 | 3.56187e-09      |
#       |            100 | 3.65299e-09      |
#       |          100.1 | 6.28941e-10      |
#       |          110.5 | 1.11456e-10      |
#       |            125 | 1.46126e-11      |
#       |          125.1 | 1.44428e-11      |
#       |            300 | 9.08858e-15      |
#       |            700 | 8.51674e-18      |
#       |           1500 | 2.86915e-19      |
#
# Scenario 03
# ===========
#
#   Values obtained from GMAT R2018a using the following inputs:
#
#       Date:      2017-01-01 00:00:00 UTC
#       Latitude:  45 deg
#       Longitude: 0 deg
#       F10.7:     100
#       F10.7ₐ:    100
#       Kp:        9
#
#   Result:
#
#       | Altitude [km]  | Density [g/cm³] |
#       |----------------|------------------|
#       |             92 | 5.55597e-09      |
#       |            100 | 5.63386e-09      |
#       |          100.1 | 9.75634e-10      |
#       |          110.5 | 1.73699e-10      |
#       |            125 | 2.41828e-11      |
#       |          125.1 | 2.39097e-11      |
#       |            300 | 3.52129e-14      |
#       |            700 | 1.28622e-16      |
#       |           1500 | 1.97775e-18      |
#
################################################################################

@testset "Function jr1971" begin
    # Common inputs to all scenarios.
    JD   = date_to_jd(2017,1,1,0,0,0)
    glat = 45*pi/180
    glon = 0.0
    F10  = 100
    F10ₐ = 100
    h    = [92; 100; 100.1; 110.5; 125; 125.1; 300; 700; 1500]*1000

    # Scenario 01
    # ===========

    Kp = 4

    # Results in [kg/m³]
    results = [3.87506e-09;
               3.96585e-09;
               6.83540e-10;
               1.21240e-10;
               1.60849e-11;
               1.58997e-11;
               1.30609e-14;
               1.34785e-17;
               4.00464e-19]*1000

    for i = 1:length(h)
        ret = jr1971(JD, glat, glon, h[i], F10, F10ₐ, Kp)
        @test ret.rho ≈ results[i] rtol = 5e-4
    end

    # Scenario 02
    # ===========

    Kp = 1

    # Results in [kg/m³]
    results = [3.56187e-09;
               3.65299e-09;
               6.28941e-10;
               1.11456e-10;
               1.46126e-11;
               1.44428e-11;
               9.08858e-15;
               8.51674e-18;
               2.86915e-19]*1000

    for i = 1:length(h)
        ret = jr1971(JD, glat, glon, h[i], F10, F10ₐ, Kp)
        @test ret.rho ≈ results[i] rtol = 5e-4
    end

    # Scenario 03
    # ===========

    Kp = 9

    # Results in [kg/m³]
    results = [5.55597e-09
               5.63386e-09
               9.75634e-10
               1.73699e-10
               2.41828e-11
               2.39097e-11
               3.52129e-14
               1.28622e-16
               1.97775e-18]*1000

    for i = 1:length(h)
        ret = jr1971(JD, glat, glon, h[i], F10, F10ₐ, Kp)
        @test ret.rho ≈ results[i] rtol = 5e-4
    end
end
