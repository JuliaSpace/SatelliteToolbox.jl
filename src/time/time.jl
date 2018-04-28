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
#   Functions related Data and Time conversion.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] https://support.microsoft.com/en-us/help/214019/method-to-determine-whether-a-year-is-a-leap-year
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-04-28: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Add conversions between Julian Day UTC, UT1, and TT.
#
# 2018-04-09: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export is_leap_year
export JD_UT1toUTC, JD_UTCtoUT1
export JD_UTCtoTT,  JD_TTtoUTC
export get_ΔAT

################################################################################
#                                  Constants
################################################################################

# Table containing the leap seconds between UTC and the International Atomic
# Time (TAI). Notice that only the dates in which an increment in the leap
# seconds occurred should be added to this table.
#
#           Julian Day         ΔAT [s]
const ΔAT_Data = [2441499.500000     11.0;
                  2441683.500000     12.0;
                  2442048.500000     13.0;
                  2442413.500000     14.0;
                  2442778.500000     15.0;
                  2443144.500000     16.0;
                  2443509.500000     17.0;
                  2443874.500000     18.0;
                  2444239.500000     19.0;
                  2444786.500000     20.0;
                  2445151.500000     21.0;
                  2445516.500000     22.0;
                  2446247.500000     23.0;
                  2447161.500000     24.0;
                  2447892.500000     25.0;
                  2448257.500000     26.0;
                  2448804.500000     27.0;
                  2449169.500000     28.0;
                  2449534.500000     29.0;
                  2450083.500000     30.0;
                  2450630.500000     31.0;
                  2451179.500000     32.0;
                  2453736.500000     33.0;
                  2454832.500000     34.0;
                  2456109.500000     35.0;
                  2457204.500000     36.0;
                  2457754.500000     37.0;
                 ]

"""
### function get_ΔAT(JD::Number)

Get the accumulated leap seconds (ΔAT) between UTC and International Atomic Time
(TAI) in the given `JD`. This function search for ΔAT in the array `ΔAT_Data`.

##### Args

* JD: Julian Day in which ΔAT will be computed.

##### Returns

The accumulated leap seconds (ΔAT) between UTC and TAI.

##### Remarks

If `JD` is before `ΔAT_Data[1,1]`, then 10 will be returned. **Notice that this
can lead to errors.**

If `JD` is after `ΔAT_Data[end,1]`, then `ΔAT_Data[end,2]` will be returned,
because it is not possible yet to predict when leap seconds will be added.

"""

function get_ΔAT(JD::Number)
    # If `JD` is before `ΔAT_Data[1,1]`, then return 10.0.
    if JD < ΔAT_Data[1,1]
        return 10.0
    else
        for i = 2:size(ΔAT_Data,1)
            @inbounds (JD < ΔAT_Data[i,1]) && return ΔAT_Data[i-1,2]
        end
    end

    # In this case, `JD` is after `ΔAT_Data[end,1]`.
    return ΔAT_Data[end,2]
end

################################################################################
#                                  Functions
################################################################################

"""
### function JD_UTCtoUT1(JD_UTC::Number, ΔUT1::Number)

Convert the Julian Day in UTC to the Julian Day in UT1 using the accumulated
difference `ΔUT1`, which is provided by IERS EOP Data.

##### Args

* JD_UTC: Julian Day [UTC].
* ΔUT1: Accumulated difference between UTC and UT1 [s].

##### Returns

The Julian Day in UT1.

"""

JD_UTCtoUT1(JD_UTC::Number, ΔUT1::Number) = JD_UTC + ΔUT1

"""
### function JD_UT1toUTC(JD_UT1::Number, ΔUT1::Number)

Convert the Julian Day in UT1 to the Julian Day in UTC using the accumulated
difference `ΔUT1`, which is provided by IERS EOP Data.

##### Args

* JD_UT1: Julian Day [UT1].
* ΔUT1: Accumulated difference between UTC and UT1 [s].

##### Returns

The Julian Day in UT1.

"""

JD_UT1toUTC(JD_UT1::Number, ΔUT1::Number) = JD_UT1 - ΔUT1

"""
### function JD_UTCtoUT1(JD_UTC::Number, eop::Union{EOPData_IAU1980,EOPData_IAU2000A})

Convert the Julian Day in UTC to the Julian Day in UT1 using the accumulated
difference given by the EOP Data `eop`. Notice that the accumulated difference
will be interpolated.

##### Args

* JD_UTC: Julian Day [UTC].
* eop: IERS EOP Data (see `get_iers_eop`).

##### Returns

The Julian Day in UT1.

"""

function JD_UTCtoUT1(JD_UTC::Number, eop::Union{EOPData_IAU1980,EOPData_IAU2000A})
	JD_UTCtoUT1(JD_UTC, eop.UT1_UTC[JD_UTC])
end

"""
### function JD_UTCtoUT1(JD_UTC::Number, eop::Union{EOPData_IAU1980,EOPData_IAU2000A})

Convert the Julian Day in UT1 to the Julian Day in UTC using the accumulated
difference given by the EOP Data `eop`. Notice that the accumulated difference
will be interpolated.

##### Args

* JD_UTC: Julian Day [UTC].
* eop: IERS EOP Data (see `get_iers_eop`).

##### Returns

The Julian Day in UTC.

"""

function JD_UT1toUTC(JD_UT1::Number, eop::Union{EOPData_IAU1980,EOPData_IAU2000A})
	JD_UT1toUTC(JD_UT1, eop.UT1_UTC[JD_UT1])
end

"""
### function JD_UTCtoTT(JD_UTC::Number [, ΔAT::Number])

Convert the Julian Day in UTC to the Julian Day in TT (Terrestrial Time) using
the accumulated difference `ΔAT` between UTC and the International Atomic Time
(TAI). If no value is provided, then the leap seconds will be obtained from the
table `ΔAT_Data`. **Notice that, in this case, if a date previous to 1973 is
provided, then a fixed value of 10 will be used, leading to wrong
computations.**

##### Args

* JD_UTC: Julian Day [UTC].
* ΔAT: (OPTIONAL) Accumulated difference between UTC and TAI [s].

##### Returns

The Julian Day in TT.

"""

JD_UTCtoTT(JD_UTC::Number, ΔAT::Number) = JD_UTC + (ΔAT + 32.184)/86400

function JD_UTCtoTT(JD_UTC::Number)
    ΔAT = get_ΔAT(JD_UTC)
    JD_UTCtoTT(JD_UTC,ΔAT)
end

"""
### function JD_TTtoUTC(JD_TT::Number, ΔAT::Number = 37)

Convert the Julian Day in TT (Terrestrial Time) to the Julian Day in UTC
(Terrestrial Time) using the accumulated difference `ΔAT` between UTC and the
International Atomic Time (TAI). If no value is provided, then the leap seconds
will be obtained from the table `ΔAT_Data`. **Notice that, in this case, if a
date previous to 1973 is provided, then a fixed value of 10 will be used,
leading to wrong computations.**

##### Args

* JD_TT: Julian Day [UTC].
* ΔAT: (OPTIONAL) Accumulated difference between UTC and TAI [s].

##### Returns

The Julian Day in UTC.

"""

JD_TTtoUTC(JD_TT::Number, ΔAT::Number) = JD_TT - (ΔAT + 32.184)/86400

function JD_TTtoUTC(JD_UTC::Number)
    ΔAT = get_ΔAT(JD_UTC)
    JD_TTtoUTC(JD_UTC,ΔAT)
end

"""
### function is_leap_year(year::Int)

Check if the year `year` is a leap year.

##### Args

* year: Year, must be > 0.

##### Returns

* **TRUE**: `year` is a leap year.
* **FALSE**: `year` is not a leap year.

##### Remarks

This algorithm was based on [3].

"""

function is_leap_year(year::Int)
    # Check if `year` is positive. This algorithm does not handle negative
    # years.
    (year < 0) && throw(ArgumentError("The year must be positive."))

    if (year % 4) == 0
        if (year % 100) == 0
            if (year % 400) == 0
                return true
            else
                return false
            end
        else
            return true
        end
    else
        return false
    end
end

