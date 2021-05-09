# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to Date and Time conversion.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] http://aa.usno.navy.mil/data/docs/JulianDate.php
#   [2] https://www.ietf.org/timezones/data/leap-seconds.list
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/time/time.jl
# ========================

# Function get_Δat
# ----------------

@testset "Function get_Δat"  begin
    # Leap seconds values obtained from [2].
    leap_secs = [2272060800	10	1 "Jan" 1972;
                 2303683200	12	1 "Jan" 1973;
                 2287785600	11	1 "Jul" 1972;
                 2335219200	13	1 "Jan" 1974;
                 2366755200	14	1 "Jan" 1975;
                 2398291200	15	1 "Jan" 1976;
                 2429913600	16	1 "Jan" 1977;
                 2461449600	17	1 "Jan" 1978;
                 2492985600	18	1 "Jan" 1979;
                 2524521600	19	1 "Jan" 1980;
                 2571782400	20	1 "Jul" 1981;
                 2603318400	21	1 "Jul" 1982;
                 2634854400	22	1 "Jul" 1983;
                 2698012800	23	1 "Jul" 1985;
                 2776982400	24	1 "Jan" 1988;
                 2840140800	25	1 "Jan" 1990;
                 2871676800	26	1 "Jan" 1991;
                 2918937600	27	1 "Jul" 1992;
                 2950473600	28	1 "Jul" 1993;
                 2982009600	29	1 "Jul" 1994;
                 3029443200	30	1 "Jan" 1996;
                 3076704000	31	1 "Jul" 1997;
                 3124137600	32	1 "Jan" 1999;
                 3345062400	33	1 "Jan" 2006;
                 3439756800	34	1 "Jan" 2009;
                 3550089600	35	1 "Jul" 2012;
                 3644697600	36	1 "Jul" 2015;
                 3692217600	37	1 "Jan" 2017;
                ]

    for i = 1:size(leap_secs,1)
        ΔAT   = leap_secs[i,2]
        day   = leap_secs[i,3]
        month = begin
            if leap_secs[i,4] == "Jan"
                1
            elseif leap_secs[i,4] == "Jul"
                7
            else
                error("Invalid array `leap_secs`.")
            end
        end
        year  = leap_secs[i,5]

        @test get_Δat(DatetoJD(year,month,day,0,0,0))         == ΔAT

        # One second before the entry in `leap_secs`, we must have `ΔAT-1` leap
        # seconds. However, this is not true for the very first line.
        (i == 1) && continue
        @test get_Δat(DatetoJD(year,month,day,0,0,0)-1/86400) == ΔAT-1
    end
end

# Functions jd_ut1_to_utc and jd_utc_to_ut1
# -----------------------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-7: Calculating Dynamical Time [1].
#
# The following date:
#
#   Mountain Standard Time Zone: May 14, 2004, 10:43
#
# is converted into the following one using ΔUT1 = -0.463326 s:
#
#   UT1: May 14, 2004, 16:42:59.5367
#
# Scenario 02
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction [1].
#
# The following UTC date:
#
#   UTC: April 6, 2004, 07:51:28.386009
#
# is converted into the following one using ΔUT1 = -0.4399619 s:
#
#   UT1: April 6, 2004, 07:51:27.946047
#
################################################################################

@testset "Functions jd_ut1_to_utc and jd_utc_to_ut1" begin

    ## Scenario 01
    ## ===========

    ΔUT1 = -0.463326

    ## jd_utc_to_ut1
    ## -----------

    # At the mentioned date, Mountain Standard Time is 6h behind UTC.
    JD_UTC = DatetoJD(2004, 5, 14, 10+6, 43, 0)
    JD_UT1 = jd_utc_to_ut1(JD_UTC, ΔUT1)

    (year, month, day, hour, minute, second) = JDtoDate(JD_UT1)

    @test year   == 2004
    @test month  == 5
    @test day    == 14
    @test hour   == 16
    @test minute == 42
    @test second  ≈ 59.5367 atol=1e-4

    ## jd_ut1_to_utc
    ## -------------

    JD_UT1 = DatetoJD(2004, 5, 14, 16, 42, 59.5367)
    JD_UTC = jd_ut1_to_utc(JD_UT1, ΔUT1)

    (year, month, day, hour, minute, second) = JDtoDate(JD_UTC)

    @test year   == 2004
    @test month  == 5
    @test day    == 14
    @test hour   == 10+6
    @test minute == 43
    @test second  ≈ 0.0000 atol=1e-4

    ## Scenario 02
    ## ===========

    ΔUT1 = -0.4399619

    ## jd_utc_to_ut1
    ## -------------

    JD_UTC = DatetoJD(2004, 4, 6, 07, 51, 28.386009)
    JD_UT1 = jd_utc_to_ut1(JD_UTC, ΔUT1)

    (year, month, day, hour, minute, second) = JDtoDate(JD_UT1)

    @test year   == 2004
    @test month  == 4
    @test day    == 6
    @test hour   == 7
    @test minute == 51
    @test second  ≈ 27.946047 atol=1e-4

    ## jd_ut1_to_utc
    ## -------------

    JD_UT1 = DatetoJD(2004, 4, 6, 07, 51, 27.946047)
    JD_UTC = jd_ut1_to_utc(JD_UT1, ΔUT1)

    (year, month, day, hour, minute, second) = JDtoDate(JD_UTC)

    @test year   == 2004
    @test month  == 4
    @test day    == 6
    @test hour   == 7
    @test minute == 51
    @test second  ≈ 28.386009 atol=1e-4
end

# Functions jd_tt_to_utc and jd_utc_to_tt
# ---------------------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-7: Calculating Dynamical Time [1].
#
# The following date:
#
#   Mountain Standard Time Zone: May 14, 2004, 10:43
#
# is converted into:
#
#   TT: May 14, 2004, 16:44:04.1840
#
# Scenario 02
# ===========
#
# Example 3-15: Performing IAU-76/FK5 reduction [1].
#
# The following UTC date:
#
#   UTC: April 6, 2004, 07:51:28.386009
#
# is converted into:
#
#   JD [TT]: 2453101.828154745
#
################################################################################

@testset "Functions jd_tt_to_utc and jd_utc_to_tt" begin

    ## Scenario 01
    ## ===========

    ## jd_utc_to_tt
    ## ------------

    # At the mentioned date, Mountain Standard Time is 6h behind UTC.
    JD_UTC = DatetoJD(2004, 5, 14, 10+6, 43, 0)
    JD_TT  = jd_utc_to_tt(JD_UTC)

    (year, month, day, hour, minute, second) = JDtoDate(JD_TT)

    @test year   == 2004
    @test month  == 5
    @test day    == 14
    @test hour   == 16
    @test minute == 44
    @test second  ≈ 04.1840 atol=1e-4

    ## jd_tt_to_utc
    ## ------------

    JD_TT  = DatetoJD(2004, 5, 14, 16, 44, 04.1840)
    JD_UTC = jd_tt_to_utc(JD_TT)

    (year, month, day, hour, minute, second) = JDtoDate(JD_UTC)

    @test year   == 2004
    @test month  == 5
    @test day    == 14
    @test hour   == 10+6

    # TODO: Fix this rounding problem.
    if minute == 42
        @test minute == 42
        @test second  ≈ 60.0000 atol=1e-4
    else
        @test minute == 43
        @test second  ≈ 0.0000 atol=1e-4
    end

    ## Scenario 02
    ## ===========

    ## jd_utc_to_tt
    ## ------------

    # At the mentioned date, Mountain Standard Time is 6h behind UTC.
    JD_UTC = DatetoJD(2004, 4, 6, 07, 51, 28.386009)
    JD_TT  = jd_utc_to_tt(JD_UTC)

    @test JD_TT ≈ 2453101.828154745 atol=1e-9

    ## jd_tt_to_utc
    ## ------------

    JD_TT  = 2453101.828154745
    JD_UTC = jd_tt_to_utc(JD_TT)

    (year, month, day, hour, minute, second) = JDtoDate(JD_UTC)

    @test year   == 2004
    @test month  == 4
    @test day    == 6
    @test hour   == 07
    @test minute == 51
    @test second  ≈ 28.386009 atol=1e-4

end

# File: ./src/time/julian_day.jl
# ==============================

# Function JDtoGMST
# -----------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Data obtained from [1]:
#
#   ╔═════════════════════╦════════════════╗
#   ║    Gregorian Day    ║   Julian Day   ║
#   ╠═════════════════════╬════════════════╣
#   ║ 1986-06-19 21:35:22 ║ 2446601.399560 ║
#   ║ 1987-05-19 04:00:00 ║ 2446934.666667 ║
#   ║ 2018-04-16 20:19:37 ║ 2458225.346956 ║
#   ║ 1822-09-07 12:00:00 ║ 2386781.000000 ║
#   ║ 1900-01-01 00:00:00 ║ 2415020.500000 ║
#   ║ 2000-01-01 00:00:00 ║ 2451544.500000 ║
#   ║ 2013-10-19 22:00:00 ║ 2456585.416667 ║
#   ╚═════════════════════╩════════════════╝
#
# Data obtained from SatelliteToolbox:
#
#   ╔═════════════════════╦════════════════╗
#   ║    Gregorian Day    ║   Julian Day   ║
#   ╠═════════════════════╬════════════════╣
#   ║ 1986-06-19 21:35:22 ║ 2446601.399560 ║
#   ║ 1987-05-19 04:00:00 ║ 2446934.666667 ║
#   ║ 2018-04-16 20:19:37 ║ 2458225.346956 ║
#   ║ 1822-09-07 12:00:00 ║ 2386781.000000 ║
#   ║ 1900-01-01 00:00:00 ║ 2415020.500000 ║
#   ║ 2000-01-01 00:00:00 ║ 2451544.500000 ║
#   ║ 2013-10-19 22:00:00 ║ 2456585.416667 ║
#   ╚═════════════════════╩════════════════╝
#
################################################################################

@testset "Function JDtoGMST" begin
    @test DatetoJD(1986,06,19,21,35,22) ≈ 2446601.399560 atol=1e-6
    @test DatetoJD(1987,05,19,04,00,00) ≈ 2446934.666667 atol=1e-6
    @test DatetoJD(2018,04,16,20,19,37) ≈ 2458225.346956 atol=1e-6
    @test DatetoJD(1822,09,07,12,00,00) ≈ 2386781.000000 atol=1e-6
    @test DatetoJD(1900,01,01,00,00,00) ≈ 2415020.500000 atol=1e-6
    @test DatetoJD(2000,01,01,00,00,00) ≈ 2451544.500000 atol=1e-6
    @test DatetoJD(2013,10,19,22,00,00) ≈ 2456585.416667 atol=1e-6
end

# Function JDtoDate
# -----------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Data obtained from [1]:
#
#   ╔═════════════════════╦════════════════╗
#   ║    Gregorian Day    ║   Julian Day   ║
#   ╠═════════════════════╬════════════════╣
#   ║ 1986-06-19 21:35:22 ║ 2446601.399560 ║
#   ║ 1987-05-19 04:00:00 ║ 2446934.666667 ║
#   ║ 2018-04-16 20:19:37 ║ 2458225.346956 ║
#   ║ 1822-09-07 12:00:00 ║ 2386781.000000 ║
#   ║ 1900-01-01 00:00:00 ║ 2415020.500000 ║
#   ║ 2000-01-01 00:00:00 ║ 2451544.500000 ║
#   ║ 2013-10-19 22:00:00 ║ 2456585.416667 ║
#   ╚═════════════════════╩════════════════╝
#
# Data obtained from SatelliteToolbox:
#
#   ╔═════════════════════╦════════════════╗
#   ║    Gregorian Day    ║   Julian Day   ║
#   ╠═════════════════════╬════════════════╣
#   ║ 1986-06-19 21:35:22 ║ 2446601.399560 ║
#   ║ 1987-05-19 04:00:00 ║ 2446934.666667 ║
#   ║ 2018-04-16 20:19:37 ║ 2458225.346956 ║
#   ║ 1822-09-07 12:00:00 ║ 2386781.000000 ║
#   ║ 1900-01-01 00:00:00 ║ 2415020.500000 ║
#   ║ 2000-01-01 00:00:00 ║ 2451544.500000 ║
#   ║ 2013-10-19 22:00:00 ║ 2456585.416667 ║
#   ╚═════════════════════╩════════════════╝
#
################################################################################

@testset "Function JDtoDate" begin
    @test JDtoDate(Int,2446601.399560) == (1986,06,19,21,35,22)
    @test JDtoDate(Int,2446934.666667) == (1987,05,19,04,00,00)
    @test JDtoDate(Int,2458225.346956) == (2018,04,16,20,19,37)
    @test JDtoDate(Int,2386781.000000) == (1822,09,07,12,00,00)
    @test JDtoDate(Int,2415020.500000) == (1900,01,01,00,00,00)
    @test JDtoDate(Int,2451544.500000) == (2000,01,01,00,00,00)
    @test JDtoDate(Int,2456585.416667) == (2013,10,19,22,00,00)
end

# Function DatetoJD
# -----------------

# Issue #19
# ==============================================================================

@testset "Issue #19" begin
    year_i = rand(1900:2100)

    year, month, day, hour, min, sec = JDtoDate(DatetoJD(year_i,12, 31, 23, 59, 59))

    @test year  == year_i
    @test month == 12
    @test day   == 31
    @test hour  == 23
    @test min   == 59
    @test sec    ≈ 59 rtol = 1e-5

    year, month, day, hour, min, sec = JDtoDate(DatetoJD(year_i, 1, 1, 0, 0, 0))

    @test year  == year_i
    @test month == 1
    @test day   == 1
    @test hour  == 0
    @test min   == 0
    @test sec    ≈ 0 rtol = 1e-5

    date_time_1 = JDtoDate(DateTime, DatetoJD(year_i, 12, 31, 23, 59, 59))
    date_time_2 = DateTime(year_i, 12, 31, 23, 59, 59)

    @test date_time_1 === date_time_2
end

# Issue #30
# ==============================================================================

@testset "Issue #30" begin
    year, month, day, h, m, s = JDtoDate(Int, 2.45754095833333e6)

    @test year  == 2016
    @test month == 6
    @test day   == 1
    @test h     == 11
    @test m     == 0
    @test s     == 0
end

# Function DatetoJD(::DateTime)
# ==============================================================================

@testset "DatetoJD with DateTime" begin
    JDdatetime = DatetoJD(DateTime(2020, 8, 14, 12, 4, 1))
    JDnums     = DatetoJD(2020, 8, 14, 12, 4, 1)
    @test JDdatetime == JDnums
end
