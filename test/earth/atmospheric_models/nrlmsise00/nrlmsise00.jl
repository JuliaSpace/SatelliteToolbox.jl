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
#   Tests related to the NRLMSISE-00 Atmospheric Model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/nrlmsise00_output.txt
#   [2] https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-06-15: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/earth/atmospheric_models/nrlmsise00
# ===============================================

# Functions: gtd7 and gtd7d
# -------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Tests available in the text file `output.txt` at [1].
#
#  DAY           172          81         172         172         172
#  UT         29000.      29000.      75000.      29000.      29000.
#  ALT          400.        400.       1000.        100.        400.
#  LAT           60.         60.         60.         60.          0.
#  LONG         -70.        -70.        -70.        -70.        -70.
#  LST           16.         16.         16.         16.         16.
#  F107A        150.        150.        150.        150.        150.
#  F107         150.        150.        150.        150.        150.
#  AP             4.          4.          4.          4.          4.
#
#  TINF      1250.54     1166.75     1239.89     1027.32     1212.40
#  TG        1241.42     1161.71     1239.89      206.89     1208.14
#  HE      6.665E+05   3.407E+06   1.124E+05   5.412E+07   1.851E+06
#  O       1.139E+08   1.586E+08   6.934E+04   1.919E+11   1.477E+08
#  N2      1.998E+07   1.391E+07   4.247E+01   6.116E+12   1.579E+07
#  O2      4.023E+05   3.263E+05   1.323E-01   1.225E+12   2.634E+05
#  AR      3.557E+03   1.560E+03   2.619E-05   6.023E+10   1.589E+03
#  H       3.475E+04   4.854E+04   2.017E+04   1.060E+07   5.816E+04
#  N       4.096E+06   4.381E+06   5.741E+03   2.616E+05   5.479E+06
#  ANM O   2.667E+04   6.957E+03   2.374E+04   0.000E+00   1.264E+03
#  RHO     4.075E-15   5.002E-15   2.757E-18   3.584E-10   4.810E-15
#
#
#  DAY           172         172         172         172         172
#  UT         29000.      29000.      29000.      29000.      29000.
#  ALT          400.        400.        400.        400.        400.
#  LAT           60.         60.         60.         60.         60.
#  LONG           0.        -70.        -70.        -70.        -70.
#  LST           16.          4.         16.         16.         16.
#  F107A        150.        150.         70.        150.        150.
#  F107         150.        150.        150.        180.        150.
#  AP             4.          4.          4.          4.         40.
#
#  TINF      1220.15     1116.39     1031.25     1306.05     1361.87
#  TG        1212.71     1113.00     1024.85     1293.37     1347.39
#  HE      8.673E+05   5.776E+05   3.740E+05   6.748E+05   5.529E+05
#  O       1.279E+08   6.979E+07   4.783E+07   1.245E+08   1.198E+08
#  N2      1.823E+07   1.237E+07   5.240E+06   2.369E+07   3.496E+07
#  O2      2.922E+05   2.493E+05   1.760E+05   4.912E+05   9.340E+05
#  AR      2.403E+03   1.406E+03   5.502E+02   4.579E+03   1.096E+04
#  H       3.686E+04   5.292E+04   8.897E+04   3.245E+04   2.686E+04
#  N       3.897E+06   1.070E+06   1.980E+06   5.371E+06   4.890E+06
#  ANM O   2.667E+04   2.667E+04   9.122E+03   2.667E+04   2.805E+04
#  RHO     4.356E-15   2.471E-15   1.572E-15   4.564E-15   4.975E-15
#
#
#  DAY           172         172         172         172         172
#  UT         29000.      29000.      29000.      29000.      29000.
#  ALT            0.         10.         30.         50.         70.
#  LAT           60.         60.         60.         60.         60.
#  LONG         -70.        -70.        -70.        -70.        -70.
#  LST           16.         16.         16.         16.         16.
#  F107A        150.        150.        150.        150.        150.
#  F107         150.        150.        150.        150.        150.
#  AP             4.          4.          4.          4.          4.
#
#  TINF      1027.32     1027.32     1027.32     1027.32     1027.32
#  TG         281.46      227.42      237.44      279.56      219.07
#  HE      1.375E+14   4.427E+13   2.128E+12   1.412E+11   1.255E+10
#  O       0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00
#  N2      2.050E+19   6.598E+18   3.171E+17   2.104E+16   1.875E+15
#  O2      5.499E+18   1.770E+18   8.506E+16   5.645E+15   4.923E+14
#  AR      2.452E+17   7.892E+16   3.793E+15   2.517E+14   2.240E+13
#  H       0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00
#  N       0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00
#  ANM O   0.000E+00   0.000E+00   0.000E+00   0.000E+00   0.000E+00
#  RHO     1.261E-03   4.059E-04   1.951E-05   1.295E-06   1.148E-07
#
################################################################################

@testset "Function gtd7d" begin
    # Scenario 01
    # ===========

    test_list = ["nrlmsise00_test01.txt",
                 "nrlmsise00_test02.txt",
                 "nrlmsise00_test03.txt"]

    for file in test_list
        raw = readdlm(file)

        # Loop for each test case.
        for i = 2:size(raw,2)
            # Parse the inputs.
            doy       = raw[1,i]
            sec       = raw[2,i]
            alt       = raw[3,i]
            g_lat     = raw[4,i]
            g_long    = raw[5,i]
            lst       = raw[6,i]
            f107A     = raw[7,i]
            f107      = raw[8,i]
            ap        = raw[9,i]

            # Parse the outputs.
            tinf      = raw[10,i]
            tg        = raw[11,i]
            den_He    = raw[12,i]
            den_O     = raw[13,i]
            den_N2    = raw[14,i]
            den_O2    = raw[15,i]
            den_Ar    = raw[16,i]
            den_H     = raw[17,i]
            den_N     = raw[18,i]
            den_aO    = raw[19,i]
            den_Total = raw[20,i]

            # Configure and run NRLMSISE-00.
            nrlmsise00d = conf_nrlmsise00(0,
                                          doy,
                                          sec,
                                          alt,
                                          g_lat,
                                          g_long,
                                          lst,
                                          f107A,
                                          f107,
                                          ap)

            out = gtd7(nrlmsise00d)

            # Test the results.
            @test out.T_exo     ≈ tinf      rtol = 1e-3 atol = 1e-9
            @test out.T_alt     ≈ tg        rtol = 1e-3 atol = 1e-9
            @test out.den_He    ≈ den_He    rtol = 1e-3 atol = 1e-9
            @test out.den_O     ≈ den_O     rtol = 1e-3 atol = 1e-9
            @test out.den_N2    ≈ den_N2    rtol = 1e-3 atol = 1e-9
            @test out.den_O2    ≈ den_O2    rtol = 1e-3 atol = 1e-9
            @test out.den_Ar    ≈ den_Ar    rtol = 1e-3 atol = 1e-9
            @test out.den_H     ≈ den_H     rtol = 1e-3 atol = 1e-9
            @test out.den_N     ≈ den_N     rtol = 1e-3 atol = 1e-9
            @test out.den_aO    ≈ den_aO    rtol = 1e-3 atol = 1e-9

            @test out.den_Total ≈ den_Total rtol = 1e-3

            # Compare the results with `gtd7d`, in which the only modification
            # is that the density of anomalous oxygen is summed in `den_Total`.

            outd = gtd7d(nrlmsise00d)

            @test outd.T_exo     ≈ tinf      rtol = 1e-3 atol = 1e-9
            @test outd.T_alt     ≈ tg        rtol = 1e-3 atol = 1e-9
            @test outd.den_He    ≈ den_He    rtol = 1e-3 atol = 1e-9
            @test outd.den_O     ≈ den_O     rtol = 1e-3 atol = 1e-9
            @test outd.den_N2    ≈ den_N2    rtol = 1e-3 atol = 1e-9
            @test outd.den_O2    ≈ den_O2    rtol = 1e-3 atol = 1e-9
            @test outd.den_Ar    ≈ den_Ar    rtol = 1e-3 atol = 1e-9
            @test outd.den_H     ≈ den_H     rtol = 1e-3 atol = 1e-9
            @test outd.den_N     ≈ den_N     rtol = 1e-3 atol = 1e-9
            @test outd.den_aO    ≈ den_aO    rtol = 1e-3 atol = 1e-9

            @test outd.den_Total ≈ den_Total + 1.66e-24*16den_aO rtol = 1e-3
        end
    end
end

# Function: nrlmsise00
# ---------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Simulation using only the daily AP.
#
# Data obtained from the online version of NRLMSISE-00 [2]:
#
# Input parameters
#
# year= 1986, month= 6, day= 19, hour=21.50,
# Time_type = Universal
# Coordinate_type = Geographic
# latitude= -16.00, longitude= 312.00, height= 100.00
# Prof. parameters: start= 0.00 stop= 1000.00 step= 100.00
#
# Optional parametes: F10.7(daily) =121.; F10.7(3-month avg) =80.; ap(daily) = 7.
#
#    Selected parameters are:
# 1 Height, km
# 2 O, cm-3
# 3 N2, cm-3
# 4 O2, cm-3
# 5 Mass_density, g/cm-3
# 6 Temperature_neutral, K
# 7 Temperature_exospheric, K
# 8 He, cm-3
# 9 Ar, cm-3
# 10 H, cm-3
# 11 N, cm-3
# 12 Anomalous_Oxygen, cm-3
#
#       1          2          3          4         5      6     7          8          9         10         11         12
#     0.0  0.000E+00  1.918E+19  5.145E+18 1.180E-03  297.7  1027  1.287E+14  2.294E+17  0.000E+00  0.000E+00  0.000E+00
#   100.0  4.244E+11  9.498E+12  2.240E+12 5.783E-10  165.9  1027  1.061E+08  9.883E+10  2.209E+07  3.670E+05  0.000E+00
#   200.0  2.636E+09  2.248E+09  1.590E+08 1.838E-13  829.3   909  5.491E+06  1.829E+06  2.365E+05  3.125E+07  1.802E-09
#   300.0  3.321E+08  6.393E+07  2.881E+06 1.212E-14  900.6   909  3.173E+06  1.156E+04  1.717E+05  6.708E+06  4.938E+00
#   400.0  5.088E+07  2.414E+06  6.838E+04 1.510E-15  907.7   909  1.979E+06  1.074E+02  1.518E+05  1.274E+06  1.237E+03
#   500.0  8.340E+06  1.020E+05  1.839E+03 2.410E-16  908.5   909  1.259E+06  1.170E+00  1.355E+05  2.608E+05  4.047E+03
#   600.0  1.442E+06  4.729E+03  5.500E+01 4.543E-17  908.6   909  8.119E+05  1.455E-02  1.214E+05  5.617E+04  4.167E+03
#   700.0  2.622E+05  2.394E+02  1.817E+00 1.097E-17  908.6   909  5.301E+05  2.050E-04  1.092E+05  1.264E+04  3.173E+03
#   800.0  4.999E+04  1.317E+01  6.608E-02 3.887E-18  908.6   909  3.503E+05  3.254E-06  9.843E+04  2.964E+03  2.246E+03
#   900.0  9.980E+03  7.852E-01  2.633E-03 1.984E-18  908.6   909  2.342E+05  5.794E-08  8.900E+04  7.237E+02  1.571E+03
#  1000.0  2.082E+03  5.055E-02  1.146E-04 1.244E-18  908.6   909  1.582E+05  1.151E-09  8.069E+04  1.836E+02  1.103E+03
#
# Scenario 02
# ===========
#
# Simulation using the AP vector.
#
# Data obtained from the online version of NRLMSISE-00 [2]:
#
# Input parameters
# year= 2016, month= 6, day= 1, hour=11.00,
# Time_type = Universal
# Coordinate_type = Geographic
# latitude= -23.00, longitude= 315.00, height= 100.00
# Prof. parameters: start= 0.00 stop= 1000.00 step= 100.00
#
# Optional parametes: F10.7(daily) =not specified; F10.7(3-month avg) =not specified; ap(daily) = not specified
#
#    Selected parameters are:
# 1 Height, km
# 2 O, cm-3
# 3 N2, cm-3
# 4 O2, cm-3
# 5 Mass_density, g/cm-3
# 6 Temperature_neutral, K
# 7 Temperature_exospheric, K
# 8 He, cm-3
# 9 Ar, cm-3
# 10 H, cm-3
# 11 N, cm-3
# 12 Anomalous_Oxygen, cm-3
# 13 F10_7_daily
# 14 F10_7_3_month_avg
# 15 ap_daily
# 16 ap_00_03_hr_prior
# 17 ap_03_06_hr_prior
# 18 ap_06_09_hr_prior
# 19 ap_09_12_hr_prior
# 20 ap_12_33_hr_prior
# 21 ap_33_59_hr_prior
#
#       1          2          3          4         5      6     7          8          9         10         11         12      13     14     15      16     17     18      19     20     21
#     0.0  0.000E+00  1.919E+19  5.149E+18 1.181E-03  296.0  1027  1.288E+14  2.296E+17  0.000E+00  0.000E+00  0.000E+00   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   100.0  6.009E+11  9.392E+12  2.213E+12 5.764E-10  165.2  1027  1.339E+08  9.580E+10  2.654E+07  2.737E+05  0.000E+00   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   200.0  3.882E+09  1.978E+09  1.341E+08 2.028E-13  687.5   706  2.048E+07  9.719E+05  5.172E+05  1.632E+07  1.399E-09   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   300.0  3.140E+08  2.481E+07  9.340E+05 9.666E-15  705.5   706  1.082E+07  1.879E+03  3.860E+05  2.205E+06  3.876E+00   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   400.0  2.855E+07  3.737E+05  7.743E+03 8.221E-16  706.0   706  5.938E+06  4.688E+00  3.317E+05  2.633E+05  9.738E+02   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   500.0  2.787E+06  6.372E+03  7.382E+01 9.764E-17  706.0   706  3.319E+06  1.397E-02  2.868E+05  3.424E+04  3.186E+03   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   600.0  2.910E+05  1.222E+02  8.047E-01 2.079E-17  706.0   706  1.887E+06  4.919E-05  2.490E+05  4.742E+03  3.281E+03   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   700.0  3.240E+04  2.622E+00  9.973E-03 8.474E-18  706.0   706  1.090E+06  2.034E-07  2.171E+05  6.946E+02  2.498E+03   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   800.0  3.835E+03  6.264E-02  1.398E-04 4.665E-18  706.0   706  6.393E+05  9.809E-10  1.900E+05  1.074E+02  1.768E+03   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#   900.0  4.816E+02  1.659E-03  2.204E-06 2.817E-18  706.0   706  3.806E+05  5.481E-12  1.669E+05  1.747E+01  1.236E+03   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#  1000.0  6.400E+01  4.853E-05  3.891E-08 1.772E-18  706.0   706  2.298E+05  3.528E-14  1.471E+05  2.988E+00  8.675E+02   89.0   89.4    4.9    3.0    4.0    4.0    6.0    9.9   10.6
#
#
################################################################################

@testset "Function nrlmsise00" begin
    ## Scenario 01
    ## ===========

    # Test outputs.
    test_out = [
       0.0  0.000E+00  1.918E+19  5.145E+18 1.180E-03  297.7  1027  1.287E+14 2.294E+17  0.000E+00  0.000E+00  0.000E+00;
     100.0  4.244E+11  9.498E+12  2.240E+12 5.783E-10  165.9  1027  1.061E+08 9.883E+10  2.209E+07  3.670E+05  0.000E+00;
     200.0  2.636E+09  2.248E+09  1.590E+08 1.838E-13  829.3   909  5.491E+06 1.829E+06  2.365E+05  3.125E+07  1.802E-09;
     300.0  3.321E+08  6.393E+07  2.881E+06 1.212E-14  900.6   909  3.173E+06 1.156E+04  1.717E+05  6.708E+06  4.938E+00;
     400.0  5.088E+07  2.414E+06  6.838E+04 1.510E-15  907.7   909  1.979E+06 1.074E+02  1.518E+05  1.274E+06  1.237E+03;
     500.0  8.340E+06  1.020E+05  1.839E+03 2.410E-16  908.5   909  1.259E+06 1.170E+00  1.355E+05  2.608E+05  4.047E+03;
     600.0  1.442E+06  4.729E+03  5.500E+01 4.543E-17  908.6   909  8.119E+05 1.455E-02  1.214E+05  5.617E+04  4.167E+03;
     700.0  2.622E+05  2.394E+02  1.817E+00 1.097E-17  908.6   909  5.301E+05 2.050E-04  1.092E+05  1.264E+04  3.173E+03;
     800.0  4.999E+04  1.317E+01  6.608E-02 3.887E-18  908.6   909  3.503E+05 3.254E-06  9.843E+04  2.964E+03  2.246E+03;
     900.0  9.980E+03  7.852E-01  2.633E-03 1.984E-18  908.6   909  2.342E+05 5.794E-08  8.900E+04  7.237E+02  1.571E+03;
    1000.0  2.082E+03  5.055E-02  1.146E-04 1.244E-18  908.6   909  1.582E+05 1.151E-09  8.069E+04  1.836E+02  1.103E+03;
    ]

    # Constant input parameters.
    year   = 1986
    month  = 6
    day    = 19
    hour   = 21
    minute = 30
    second = 00
    g_lat  = -16*pi/180
    g_long = 312*pi/180
    f107   = 121
    f107A  = 80
    ap     = 7
    JD     = DatetoJD(year, month, day, hour, minute, second)

    for i = 1:size(test_out,1)
        # Run the NRLMSISE-00 model wih the input parameters.
        out = nrlmsise00(JD,
                         test_out[i,1]*1000,
                         g_lat,
                         g_long,
                         f107A,
                         f107,
                         ap;
                         output_si = false,
                         dversion = false)

        @test out.den_O     ≈ test_out[i, 2]  rtol = 1e-3 atol = 1e-9
        @test out.den_N2    ≈ test_out[i, 3]  rtol = 1e-3 atol = 1e-9
        @test out.den_O2    ≈ test_out[i, 4]  rtol = 1e-3 atol = 1e-9
        @test out.T_alt     ≈ test_out[i, 6]  rtol = 1e-1 atol = 1e-9
        @test out.T_exo     ≈ test_out[i, 7]  rtol = 1e-1 atol = 1e-9
        @test out.den_He    ≈ test_out[i, 8]  rtol = 1e-3 atol = 1e-9
        @test out.den_Ar    ≈ test_out[i, 9]  rtol = 1e-3 atol = 1e-9
        @test out.den_H     ≈ test_out[i,10]  rtol = 1e-3 atol = 1e-9
        @test out.den_N     ≈ test_out[i,11]  rtol = 1e-3 atol = 1e-9
        @test out.den_aO    ≈ test_out[i,12]  rtol = 1e-3 atol = 1e-9
        @test out.den_Total ≈ test_out[i, 5]  rtol = 1e-3

        # Test the version that calls `gtd7d` instead of `gtd7`.
        outd = nrlmsise00(JD,
                          test_out[i,1]*1000,
                          g_lat,
                          g_long,
                          f107A,
                          f107,
                          ap;
                          output_si = false,
                          dversion = true)

        @test outd.den_O     ≈ test_out[i, 2]  rtol = 1e-3 atol = 1e-9
        @test outd.den_N2    ≈ test_out[i, 3]  rtol = 1e-3 atol = 1e-9
        @test outd.den_O2    ≈ test_out[i, 4]  rtol = 1e-3 atol = 1e-9
        @test outd.T_alt     ≈ test_out[i, 6]  rtol = 1e-1 atol = 1e-9
        @test outd.T_exo     ≈ test_out[i, 7]  rtol = 1e-1 atol = 1e-9
        @test outd.den_He    ≈ test_out[i, 8]  rtol = 1e-3 atol = 1e-9
        @test outd.den_Ar    ≈ test_out[i, 9]  rtol = 1e-3 atol = 1e-9
        @test outd.den_H     ≈ test_out[i,10]  rtol = 1e-3 atol = 1e-9
        @test outd.den_N     ≈ test_out[i,11]  rtol = 1e-3 atol = 1e-9
        @test outd.den_aO    ≈ test_out[i,12]  rtol = 1e-3 atol = 1e-9

        @test outd.den_Total ≈ test_out[i, 5] + 1.66e-24*16outd.den_aO rtol = 1e-3
    end

    ## Scenario 02
    ## ===========

    # Test outputs.
    test_out = [
       0.0  0.000E+00  1.919E+19  5.149E+18 1.181E-03  296.0  1027  1.288E+14 2.296E+17  0.000E+00  0.000E+00  0.000E+00   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6;
     100.0  6.009E+11  9.392E+12  2.213E+12 5.764E-10  165.2  1027  1.339E+08 9.580E+10  2.654E+07  2.737E+05  0.000E+00   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6;
     200.0  3.882E+09  1.978E+09  1.341E+08 2.028E-13  687.5   706  2.048E+07 9.719E+05  5.172E+05  1.632E+07  1.399E-09   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6;
     300.0  3.140E+08  2.481E+07  9.340E+05 9.666E-15  705.5   706  1.082E+07 1.879E+03  3.860E+05  2.205E+06  3.876E+00   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6;
     400.0  2.855E+07  3.737E+05  7.743E+03 8.221E-16  706.0   706  5.938E+06 4.688E+00  3.317E+05  2.633E+05  9.738E+02   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6;
     500.0  2.787E+06  6.372E+03  7.382E+01 9.764E-17  706.0   706  3.319E+06 1.397E-02  2.868E+05  3.424E+04  3.186E+03   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6;
     600.0  2.910E+05  1.222E+02  8.047E-01 2.079E-17  706.0   706  1.887E+06 4.919E-05  2.490E+05  4.742E+03  3.281E+03   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6;
     700.0  3.240E+04  2.622E+00  9.973E-03 8.474E-18  706.0   706  1.090E+06 2.034E-07  2.171E+05  6.946E+02  2.498E+03   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6;
     800.0  3.835E+03  6.264E-02  1.398E-04 4.665E-18  706.0   706  6.393E+05 9.809E-10  1.900E+05  1.074E+02  1.768E+03   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6;
     900.0  4.816E+02  1.659E-03  2.204E-06 2.817E-18  706.0   706  3.806E+05 5.481E-12  1.669E+05  1.747E+01  1.236E+03   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6;
    1000.0  6.400E+01  4.853E-05  3.891E-08 1.772E-18  706.0   706  2.298E+05 3.528E-14  1.471E+05  2.988E+00  8.675E+02   89.0   89.4    4.9    3.0 4.0    4.0    6.0    9.9   10.6;
    ]

    # Constant input parameters.
    year   = 2016
    month  = 6
    day    = 1
    hour   = 11
    minute = 00
    second = 00
    g_lat  = -23*pi/180
    g_long = -45*pi/180
    JD     = DatetoJD(year, month, day, hour, minute, second)

    # TODO: The tolerances here are much higher than in the previous test with
    # the daily AP only. This must be analyzed. However, it was observed that in
    # the NRLMSISE-00 version, the last component of the AP array is the
    # "Average of eight 3 hour AP indices from 36 to 57 hours prior to current
    # time." Whereas the label from the online version is "ap_33_59_hr_prior",
    # which seems to be different. On the other hand, the MSISE90 describe the
    # last vector of AP array as "Average of eight 3 hr AP indicies [sic] from
    # 36 to 59 hrs prior to current time." Hence, it seems that the online
    # version is using the old algorithm related to the AP array. We must change
    # the algorithm in our implementation to the old one and see if the result
    # is more similar to the online version.
    #
    # The information can be obtained at:
    #
    #   https://ccmc.gsfc.nasa.gov/modelweb/
    #   https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/msise90/
    #

    for i = 1:size(test_out,1)
        f107  = test_out[i, 13]
        f107A = test_out[i, 14]
        ap_a  = test_out[i, 15:21]

        # Run the NRLMSISE-00 model wih the input parameters.
        out = nrlmsise00(JD,
                         test_out[i,1]*1000,
                         g_lat,
                         g_long,
                         f107A,
                         f107,
                         ap_a;
                         output_si = false,
                         dversion = false)

        @test out.den_O     ≈ test_out[i, 2]  rtol = 5e-2 atol = 1e-9
        @test out.den_N2    ≈ test_out[i, 3]  rtol = 5e-2 atol = 1e-9
        @test out.den_O2    ≈ test_out[i, 4]  rtol = 5e-2 atol = 1e-9
        @test out.T_alt     ≈ test_out[i, 6]  rtol = 1e-1 atol = 1e-9
        @test out.T_exo     ≈ test_out[i, 7]  rtol = 1e-1 atol = 1e-9
        @test out.den_He    ≈ test_out[i, 8]  rtol = 5e-2 atol = 1e-9
        @test out.den_Ar    ≈ test_out[i, 9]  rtol = 5e-2 atol = 1e-9
        @test out.den_H     ≈ test_out[i,10]  rtol = 5e-2 atol = 1e-9
        @test out.den_N     ≈ test_out[i,11]  rtol = 5e-2 atol = 1e-9
        @test out.den_aO    ≈ test_out[i,12]  rtol = 9e-2 atol = 1e-9
        @test out.den_Total ≈ test_out[i, 5]  rtol = 5e-2

        # Test the version that calls `gtd7d` instead of `gtd7`.
        outd = nrlmsise00(JD,
                          test_out[i,1]*1000,
                          g_lat,
                          g_long,
                          f107A,
                          f107,
                          ap_a;
                          output_si = false,
                          dversion = true)

        @test outd.den_O     ≈ test_out[i, 2]  rtol = 5e-2 atol = 1e-9
        @test outd.den_N2    ≈ test_out[i, 3]  rtol = 5e-2 atol = 1e-9
        @test outd.den_O2    ≈ test_out[i, 4]  rtol = 5e-2 atol = 1e-9
        @test outd.T_alt     ≈ test_out[i, 6]  rtol = 1e-1 atol = 1e-9
        @test outd.T_exo     ≈ test_out[i, 7]  rtol = 1e-1 atol = 1e-9
        @test outd.den_He    ≈ test_out[i, 8]  rtol = 5e-2 atol = 1e-9
        @test outd.den_Ar    ≈ test_out[i, 9]  rtol = 5e-2 atol = 1e-9
        @test outd.den_H     ≈ test_out[i,10]  rtol = 5e-2 atol = 1e-9
        @test outd.den_N     ≈ test_out[i,11]  rtol = 5e-2 atol = 1e-9
        @test outd.den_aO    ≈ test_out[i,12]  rtol = 9e-2 atol = 1e-9

        @test outd.den_Total ≈ test_out[i, 5] + 1.66e-24*16outd.den_aO rtol = 5e-2
    end
end

