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
#   Tests related to IGRF model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-05-15: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/earth/geomagnetic_field_models
# ==========================================

# Load the test files.
#
# The test values were created using the original `igrf12syn` function. The
# source code of the program that generated those files can be seen on:
#
#   https://github.com/SatelliteToolbox/IGRF_Test
#

igrf12_geocentric_test = readdlm("./IGRF12_test_geocentric.txt")
igrf12_geodetic_test   = readdlm("./IGRF12_test_geodetic.txt")

# Function igrf12syn
# ------------------

@testset "Function igrf12syn" begin
    # Testing the geocentric part of the algorithm.
    for i = 1:size(igrf12_geocentric_test, 1)
        date  = igrf12_geocentric_test[i,1]
        r     = igrf12_geocentric_test[i,2]
        colat = igrf12_geocentric_test[i,3]
        elong = igrf12_geocentric_test[i,4]
        xt    = igrf12_geocentric_test[i,5]
        yt    = igrf12_geocentric_test[i,6]
        zt    = igrf12_geocentric_test[i,7]
        ft    = igrf12_geocentric_test[i,8]

        # Call IGRF with the same inputs as those in the test.
        (x, y, z, f) = igrf12syn(0, date, 2, r, colat, elong)

        # Test the values.
        @test x ≈ xt atol=3e-1
        @test y ≈ yt atol=3e-1
        @test z ≈ zt atol=3e-1
        @test f ≈ ft atol=3e-1
    end

    # Testing the geodetic part of the algorithm.
    for i = 1:size(igrf12_geodetic_test, 1)
        date  = igrf12_geodetic_test[i,1]
        h     = igrf12_geodetic_test[i,2]
        colat = igrf12_geodetic_test[i,3]
        elong = igrf12_geodetic_test[i,4]
        xt    = igrf12_geodetic_test[i,5]
        yt    = igrf12_geodetic_test[i,6]
        zt    = igrf12_geodetic_test[i,7]
        ft    = igrf12_geodetic_test[i,8]

        # Call IGRF with the same inputs as those in the test.
        (x, y, z, f) = igrf12syn(0, date, 1, h, colat, elong)

        # Test the values.
        @test x ≈ xt atol=3e-1
        @test y ≈ yt atol=3e-1
        @test z ≈ zt atol=3e-1
        @test f ≈ ft atol=3e-1
    end
end

# Function igrf12
# ---------------

@testset "Function igrf12" begin
    # Testing the geocentric part of the algorithm.
    for i = 1:size(igrf12_geocentric_test, 1)
        date  = igrf12_geocentric_test[i,1]
        r     = igrf12_geocentric_test[i,2]
        colat = igrf12_geocentric_test[i,3]
        elong = igrf12_geocentric_test[i,4]
        xt    = igrf12_geocentric_test[i,5]
        yt    = igrf12_geocentric_test[i,6]
        zt    = igrf12_geocentric_test[i,7]
        ft    = igrf12_geocentric_test[i,8]

        # Call IGRF with the same inputs as those in the test.
        (elong > 180) && (elong = elong-360)
        B = igrf12(date, r*1000, (90-colat)*pi/180, elong*pi/180)

        x, y, z = B[:]
        f       = norm(B)

        # Test the values.
        @test x ≈ xt atol=3e-1
        @test y ≈ yt atol=3e-1
        @test z ≈ zt atol=3e-1
        @test f ≈ ft atol=3e-1
    end

    # Testing the geodetic part of the algorithm.
    for i = 1:size(igrf12_geodetic_test, 1)
        date  = igrf12_geodetic_test[i,1]
        h     = igrf12_geodetic_test[i,2]
        colat = igrf12_geodetic_test[i,3]
        elong = igrf12_geodetic_test[i,4]
        xt    = igrf12_geodetic_test[i,5]
        yt    = igrf12_geodetic_test[i,6]
        zt    = igrf12_geodetic_test[i,7]
        ft    = igrf12_geodetic_test[i,8]

        # Call IGRF with the same inputs as those in the test.
        (elong > 180) && (elong = elong-360)
        B = igrf12(date, h*1000, (90-colat)*pi/180, elong*pi/180, Val{:geodetic})

        x, y, z = B[:]
        f       = norm(B)

        # Test the values.
        @test x ≈ xt atol=3e-1
        @test y ≈ yt atol=3e-1
        @test z ≈ zt atol=3e-1
        @test f ≈ ft atol=3e-1
    end
end

