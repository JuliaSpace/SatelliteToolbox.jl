#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to IGRF model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/earth/geomagnetic_field_models/igrf
# ===============================================

# Load the test files.
#
# The test values were created using the original `igrfXXsyn` functions, where
# `XX` is the IGRF version. The source code of the program that generated those
# files can be seen on:
#
#   https://github.com/SatelliteToolbox/IGRF_Test
#

const igrf12_geocentric_test = readdlm("./IGRF12_test_geocentric.txt")
const igrf12_geodetic_test   = readdlm("./IGRF12_test_geodetic.txt")
const igrf13_geocentric_test = readdlm("./IGRF13_test_geocentric.txt")
const igrf13_geodetic_test   = readdlm("./IGRF13_test_geodetic.txt")

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
        if date <= 2020
            (x, y, z, f) = igrf12syn(0, date, 2, r, colat, elong)
        else
            (x, y, z, f) = @test_logs (:warn, "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2020.") #=
                =# igrf12syn(0, date, 2, r, colat, elong)
        end

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
        if date <= 2020
            (x, y, z, f) = igrf12syn(0, date, 1, h, colat, elong)
        else
            (x, y, z, f) = @test_logs (:warn, "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2020.") #=
                =# igrf12syn(0, date, 1, h, colat, elong)
        end

        # Test the values.
        @test x ≈ xt atol=3e-1
        @test y ≈ yt atol=3e-1
        @test z ≈ zt atol=3e-1
        @test f ≈ ft atol=3e-1
    end
end

# Function igrf13syn
# ------------------

@testset "Function igrf13syn" begin
    # Testing the geocentric part of the algorithm.
    for i = 1:size(igrf13_geocentric_test, 1)
        date  = igrf13_geocentric_test[i,1]
        r     = igrf13_geocentric_test[i,2]
        colat = igrf13_geocentric_test[i,3]
        elong = igrf13_geocentric_test[i,4]
        xt    = igrf13_geocentric_test[i,5]
        yt    = igrf13_geocentric_test[i,6]
        zt    = igrf13_geocentric_test[i,7]
        ft    = igrf13_geocentric_test[i,8]

        # Call IGRF with the same inputs as those in the test.
        if date <= 2025
            (x, y, z, f) = igrf13syn(0, date, 2, r, colat, elong)
        else
            (x, y, z, f) = @test_logs (:warn, "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.") #=
                =# igrf13syn(0, date, 2, r, colat, elong)
        end

        # Test the values.
        @test x ≈ xt atol=3e-1
        @test y ≈ yt atol=3e-1
        @test z ≈ zt atol=3e-1
        @test f ≈ ft atol=3e-1
    end

    # Testing the geodetic part of the algorithm.
    for i = 1:size(igrf13_geodetic_test, 1)
        date  = igrf13_geodetic_test[i,1]
        h     = igrf13_geodetic_test[i,2]
        colat = igrf13_geodetic_test[i,3]
        elong = igrf13_geodetic_test[i,4]
        xt    = igrf13_geodetic_test[i,5]
        yt    = igrf13_geodetic_test[i,6]
        zt    = igrf13_geodetic_test[i,7]
        ft    = igrf13_geodetic_test[i,8]

        # Call IGRF with the same inputs as those in the test.
        if date <= 2025
            (x, y, z, f) = igrf13syn(0, date, 1, h, colat, elong)
        else
            (x, y, z, f) = @test_logs (:warn, "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.") #=
                =# igrf13syn(0, date, 1, h, colat, elong)
        end

        # Test the values.
        @test x ≈ xt atol=3e-1
        @test y ≈ yt atol=3e-1
        @test z ≈ zt atol=3e-1
        @test f ≈ ft atol=3e-1
    end
end

# Function igrf
# -------------

@testset "Function igrf" begin
    # Auxiliary variables to use the version without allocations.
    P  = Matrix{Float64}(undef, 14, 14)
    dP = similar(P)

    # Testing the geocentric part of the algorithm.
    for i = 1:size(igrf13_geocentric_test, 1)
        date  = igrf13_geocentric_test[i,1]
        r     = igrf13_geocentric_test[i,2]
        colat = igrf13_geocentric_test[i,3]
        elong = igrf13_geocentric_test[i,4]
        xt    = igrf13_geocentric_test[i,5]
        yt    = igrf13_geocentric_test[i,6]
        zt    = igrf13_geocentric_test[i,7]
        ft    = igrf13_geocentric_test[i,8]

        # Call IGRF with the same inputs as those in the test.
        (elong > 180) && (elong = elong-360)

        if date <= 2025
            Ba = igrf(date, r*1000, (90-colat)*pi/180, elong*pi/180)
            Bn = igrf(date, r*1000, (90-colat)*pi/180, elong*pi/180, P, dP)
        else
            Ba = @test_logs (:warn, "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.") #=
                 =# igrf(date, r*1000, (90-colat)*pi/180, elong*pi/180)
            Bn = @test_logs (:warn, "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.") #=
                 =# igrf(date, r*1000, (90-colat)*pi/180, elong*pi/180, P, dP)
        end

        x, y, z = Ba[:]
        f       = norm(Ba)

        # Test the values.
        @test x   ≈ xt atol=3e-1
        @test y   ≈ yt atol=3e-1
        @test z   ≈ zt atol=3e-1
        @test f   ≈ ft atol=3e-1
        @test Ba == Bn
    end

    # Testing the geodetic part of the algorithm.
    for i = 1:size(igrf13_geodetic_test, 1)
        date  = igrf13_geodetic_test[i,1]
        h     = igrf13_geodetic_test[i,2]
        colat = igrf13_geodetic_test[i,3]
        elong = igrf13_geodetic_test[i,4]
        xt    = igrf13_geodetic_test[i,5]
        yt    = igrf13_geodetic_test[i,6]
        zt    = igrf13_geodetic_test[i,7]
        ft    = igrf13_geodetic_test[i,8]

        # Call IGRF with the same inputs as those in the test.
        (elong > 180) && (elong = elong-360)

        if date <= 2025
            Ba = igrf(date, h*1000, (90-colat)*pi/180, elong*pi/180, Val(:geodetic))
            Bn = igrf(date, h*1000, (90-colat)*pi/180, elong*pi/180, Val(:geodetic))
        else
            Ba = @test_logs (:warn, "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.") #=
                 =# igrf(date, h*1000, (90-colat)*pi/180, elong*pi/180, Val(:geodetic))
            Bn = @test_logs (:warn, "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.") #=
                 =# igrf(date, h*1000, (90-colat)*pi/180, elong*pi/180, Val(:geodetic), P, dP)
        end

        x, y, z = Ba[:]
        f       = norm(Ba)

        # Test the values.
        @test x   ≈ xt atol=3e-1
        @test y   ≈ yt atol=3e-1
        @test z   ≈ zt atol=3e-1
        @test f   ≈ ft atol=3e-1
        @test Ba == Bn
    end

    # Reduced degree
    # --------------------------------------------------------------------------

    # Geocentric
    # ----------

    # If `max_degree` is equal or lower than 0, then the results must be the
    # same as if `max_degree` is 1.
    Bref = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree =  1)
    B1   = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree =  0)
    B2   = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree = -1)

    @test Bref == B1
    @test Bref == B2

    # If `max_degree` is higher than the number of available coefficients, then
    # it must be clamped.
    Bref = igrf(2020.4452, 6515e3, 0.45, -1.34)
    B1   = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree = 13)
    B2   = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree = 22)

    @test Bref == B1
    @test Bref == B2

    # The difference between the field computed using 13 or 12 coefficients must
    # be small.
    B1 = igrf(2020.4452, 6515e3, 0.45, -1.34)
    B2 = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree = 12)

    @test acosd(dot(B1 / norm(B1), B2 / norm(B2))) < 0.01

    # Testing the version in which we provide the matrices to the associated
    # Legendre functions.
    Pred  = zeros(5, 5)
    dPred = zeros(5, 5)

    Bref = igrf(2020.4452, 6515e3, 0.45, -1.34; max_degree = 4)
    B1   = igrf(2020.4452, 6515e3, 0.45, -1.34, Pred, dPred; max_degree = 4)

    @test B1 == Bref

    # Geodetic
    # --------

    # If `max_degree` is equal or lower than 0, then the results must be the
    # same as if `max_degree` is 1.
    Bref = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree =  1)
    B1   = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree =  0)
    B2   = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree = -1)

    @test Bref == B1
    @test Bref == B2

    # If `max_degree` is higher than the number of available coefficients, then
    # it must be clamped.
    Bref = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic))
    B1   = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree = 13)
    B2   = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree = 22)

    @test Bref == B1
    @test Bref == B2

    # The difference between the field computed using 13 or 12 coefficients must
    # be small.
    B1 = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic))
    B2 = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree = 12)

    @test acosd(dot(B1 / norm(B1), B2 / norm(B2))) < 0.01

    # Testing the version in which we provide the matrices to the associated
    # Legendre functions.
    Pred  = zeros(5, 5)
    dPred = zeros(5, 5)

    Bref = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic); max_degree = 4)
    B1   = igrf(2020.4452, 752e3, 0.45, -1.34, Val(:geodetic), Pred, dPred; max_degree = 4)

    @test B1 == Bref

    # Errors
    # --------------------------------------------------------------------------

    P₀  = zeros(10,10)
    dP₀ = zeros(10,10)
    P₁  = zeros(15,15)
    dP₁ = zeros(15,15)

    @test_throws ErrorException igrf(2020, 140e3, pi/4, pi/2, Val(:geodetic), P₀, dP₀)
    @test_throws ErrorException igrf(2020, 140e3, pi/4, pi/2, Val(:geodetic), P₀, dP₁)
    @test_throws ErrorException igrf(2020, 140e3, pi/4, pi/2, Val(:geodetic), P₁, dP₀)
    @test_nowarn                igrf(2020, 140e3, pi/4, pi/2, Val(:geodetic), P₁, dP₁)

    @test_throws ErrorException igrf(2020, R0+140e3, pi/4, pi/2, Val(:geocentric), P₀, dP₀)
    @test_throws ErrorException igrf(2020, R0+140e3, pi/4, pi/2, Val(:geocentric), P₀, dP₁)
    @test_throws ErrorException igrf(2020, R0+140e3, pi/4, pi/2, Val(:geocentric), P₁, dP₀)
    @test_nowarn                igrf(2020, R0+140e3, pi/4, pi/2, Val(:geocentric), P₁, dP₁)

    # Issues
    # --------------------------------------------------------------------------

    # Calculation close to the geographic pole.
    B = igrf(2019, 7150e3, π / 2 - 1e-15, 0.55)
    @test B[1] ≈ 907.752507827486
    @test B[2] ≈ 173.19657970935657
    @test B[3] ≈ 41139.95114358637

    B = igrf(2019, 7150e3, π / 2, 0.55)
    @test B[1] ≈ 907.7525078274671
    @test B[2] ≈ 173.1965797093516
    @test B[3] ≈ 41139.95114358636
end

# Function igrfd
# --------------

@testset "Function igrfd" begin
        # Auxiliary variables to use the version without allocations.
    P  = Matrix{Float64}(undef, 14, 14)
    dP = similar(P)

    # Testing the geocentric part of the algorithm.
    for i = 1:size(igrf13_geocentric_test, 1)
        date  = igrf13_geocentric_test[i,1]
        r     = igrf13_geocentric_test[i,2]
        colat = igrf13_geocentric_test[i,3]
        elong = igrf13_geocentric_test[i,4]
        xt    = igrf13_geocentric_test[i,5]
        yt    = igrf13_geocentric_test[i,6]
        zt    = igrf13_geocentric_test[i,7]
        ft    = igrf13_geocentric_test[i,8]

        # Call IGRF with the same inputs as those in the test.
        (elong > 180) && (elong = elong-360)

        if date <= 2025
            Ba = igrfd(date, r*1000, 90-colat, elong)
            Bn = igrfd(date, r*1000, 90-colat, elong, P, dP)
        else
            Ba = @test_logs (:warn, "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.") #=
                 =# igrfd(date, r*1000, 90-colat, elong)
            Bn = @test_logs (:warn, "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.") #=
                 =# igrfd(date, r*1000, 90-colat, elong, P, dP)
        end

        x, y, z = Ba[:]
        f       = norm(Ba)

        # Test the values.
        @test x   ≈ xt atol=3e-1
        @test y   ≈ yt atol=3e-1
        @test z   ≈ zt atol=3e-1
        @test f   ≈ ft atol=3e-1
        @test Ba == Bn
    end

    # Testing the geodetic part of the algorithm.
    for i = 1:size(igrf13_geodetic_test, 1)
        date  = igrf13_geodetic_test[i,1]
        h     = igrf13_geodetic_test[i,2]
        colat = igrf13_geodetic_test[i,3]
        elong = igrf13_geodetic_test[i,4]
        xt    = igrf13_geodetic_test[i,5]
        yt    = igrf13_geodetic_test[i,6]
        zt    = igrf13_geodetic_test[i,7]
        ft    = igrf13_geodetic_test[i,8]

        # Call IGRF with the same inputs as those in the test.
        (elong > 180) && (elong = elong-360)

        if date <= 2025
            Ba = igrfd(date, h*1000, 90-colat, elong, Val(:geodetic))
            Bn = igrfd(date, h*1000, 90-colat, elong, Val(:geodetic), P, dP)
        else
            Ba = @test_logs (:warn, "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.") #=
                 =# igrfd(date, h*1000, 90-colat, elong, Val(:geodetic))
            Bn = @test_logs (:warn, "The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.") #=
                 =# igrfd(date, h*1000, 90-colat, elong, Val(:geodetic), P, dP)
        end

        x, y, z = Ba[:]
        f       = norm(Ba)

        # Test the values.
        @test x   ≈ xt atol=3e-1
        @test y   ≈ yt atol=3e-1
        @test z   ≈ zt atol=3e-1
        @test f   ≈ ft atol=3e-1
        @test Ba == Bn
    end

    # Reduced degree
    # --------------------------------------------------------------------------

    # Geocentric
    # ----------

    # If `max_degree` is equal or lower than 0, then the results must be the
    # same as if `max_degree` is 1.
    Bref = igrfd(2020.4452, 6515e3, -19, 106; max_degree =  1)
    B1   = igrfd(2020.4452, 6515e3, -19, 106; max_degree =  0)
    B2   = igrfd(2020.4452, 6515e3, -19, 106; max_degree = -1)

    @test Bref == B1
    @test Bref == B2

    # If `max_degree` is higher than the number of available coefficients, then
    # it must be clamped.
    Bref = igrfd(2020.4452, 6515e3, -19, 106)
    B1   = igrfd(2020.4452, 6515e3, -19, 106; max_degree = 13)
    B2   = igrfd(2020.4452, 6515e3, -19, 106; max_degree = 22)

    @test Bref == B1
    @test Bref == B2

    # The difference between the field computed using 13 or 12 coefficients must
    # be small.
    B1 = igrfd(2020.4452, 6515e3, -19, 106)
    B2 = igrfd(2020.4452, 6515e3, -19, 106; max_degree = 12)

    @test acosd(dot(B1 / norm(B1), B2 / norm(B2))) < 0.01

    # Testing the version in which we provide the matrices to the associated
    # Legendre functions.
    Pred  = zeros(5, 5)
    dPred = zeros(5, 5)

    Bref = igrfd(2020.4452, 6515e3, -19, 106; max_degree = 4)
    B1   = igrfd(2020.4452, 6515e3, -19, 106, Pred, dPred; max_degree = 4)

    @test B1 == Bref

    # Geodetic
    # --------

    # If `max_degree` is equal or lower than 0, then the results must be the
    # same as if `max_degree` is 1.
    Bref = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree =  1)
    B1   = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree =  0)
    B2   = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree = -1)

    @test Bref == B1
    @test Bref == B2

    # If `max_degree` is higher than the number of available coefficients, then
    # it must be clamped.
    Bref = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic))
    B1   = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree = 13)
    B2   = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree = 22)

    @test Bref == B1
    @test Bref == B2

    # The difference between the field computed using 13 or 12 coefficients must
    # be small.
    B1 = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic))
    B2 = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree = 12)

    @test acosd(dot(B1 / norm(B1), B2 / norm(B2))) < 0.01

    # Testing the version in which we provide the matrices to the associated
    # Legendre functions.
    Pred  = zeros(5, 5)
    dPred = zeros(5, 5)

    Bref = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic); max_degree = 4)
    B1   = igrfd(2020.4452, 752e3, -19, 106, Val(:geodetic), Pred, dPred; max_degree = 4)

    @test B1 == Bref

    # Errors
    # --------------------------------------------------------------------------

    @test_throws ErrorException igrfd(2010, 7000, -91,    0)
    @test_throws ErrorException igrfd(2010, 7000, +91,    0)
    @test_throws ErrorException igrfd(2010, 7000, 0, -181)
    @test_throws ErrorException igrfd(2010, 7000, 0, +181)

    P₀  = zeros(10,10)
    dP₀ = zeros(10,10)
    P₁  = zeros(15,15)
    dP₁ = zeros(15,15)

    @test_throws ErrorException igrfd(2020, 140e3, 45, 90, Val(:geodetic), P₀, dP₀)
    @test_throws ErrorException igrfd(2020, 140e3, 45, 90, Val(:geodetic), P₀, dP₁)
    @test_throws ErrorException igrfd(2020, 140e3, 45, 90, Val(:geodetic), P₁, dP₀)
    @test_nowarn                igrfd(2020, 140e3, 45, 90, Val(:geodetic), P₁, dP₁)

    @test_throws ErrorException igrfd(2020, R0+140e3, 45, 90, Val(:geocentric), P₀, dP₀)
    @test_throws ErrorException igrfd(2020, R0+140e3, 45, 90, Val(:geocentric), P₀, dP₁)
    @test_throws ErrorException igrfd(2020, R0+140e3, 45, 90, Val(:geocentric), P₁, dP₀)
    @test_nowarn                igrfd(2020, R0+140e3, 45, 90, Val(:geocentric), P₁, dP₁)
end
