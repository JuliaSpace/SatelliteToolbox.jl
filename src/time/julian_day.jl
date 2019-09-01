#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
#   [2] https://quasar.as.utexas.edu/BillInfo/JulianDatesG.html
#
#   [3] https://support.microsoft.com/en-us/help/214019/method-to-determine-whether-a-year-is-a-leap-year
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Remarks
#
# Information about the Julian Day obtained from [2] (Accessed on 2018-04-11).
# ============================================================================
#
# The Julian Day Count is a uniform count of days from a remote epoch in the
# past (-4712 January 1, 12 hours Greenwich Mean Time (Julian proleptic
# Calendar) = 4713 BCE January 1, 12 hours GMT (Julian proleptic Calendar) =
# 4714 BCE November 24, 12 hours GMT (Gregorian proleptic Calendar)). At this
# instant, the Julian Day Number is 0. It is convenient for astronomers to use
# since it is not necessary to worry about odd numbers of days in a month,
# leap years, etc. Once you have the Julian Day Number of a particular date
# in history, it is easy to calculate time elapsed between it and any other
# Julian Day Number.
#
# The Julian Day Count has nothing to do with the Julian Calendar introduced by
# Julius Caesar. It is named for Julius Scaliger, the father of Josephus Justus
# Scaliger, who invented the concept. It can also be thought of as a logical
# follow-on to the old Egyptian civil calendar, which also used years of
# constant lengths.
#
# Scaliger chose the particular date in the remote past because it was before
# recorded history and because in that year, three important cycles coincided
# with their first year of the cycle: The 19-year Metonic Cycle, the 15-year
# Indiction Cycle (a Roman Taxation Cycle) and the 28-year Solar Cycle (the
# length of time for the old Julian Calendar to repeat exactly).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export DatetoJD, JDtoDate

################################################################################
#                                  Functions
################################################################################

"""
    function DatetoJD(Y::Int, M::Int, D::Int, h::Int, m::Int, s::Number)

Convert a date represented using the Gregorian Calendar (Year = `y`, Month = `M`
(1-12), Day = `D`, Hour = `h` (0-24), minute = `m`, and second = `s`) to Julian
Day.

# Remarks

The algorithm was obtained from \\[2] (Accessed on 2018-04-11).

"""
function DatetoJD(Y::Int, M::Int, D::Int, h::Int, m::Int, s::Number)
    # Check the input.
    ( (M < 1) || (M > 12) ) && throw(ArgumentError("Invalid month. It must be an integer between 1 and 12."))
    ( (D < 1) || (D > 31) ) && throw(ArgumentError("Invalid day. It must be an integer between 1 and 31."))
    ( (h < 0) || (h > 23) ) && throw(ArgumentError("Invalid hour. It must be an integer between 0 and 23."))
    ( (m < 0) || (m > 59) ) && throw(ArgumentError("Invalid minute. It must be an integer between 0 and 59."))
    ( (s < 0) || (s > 60) ) && throw(ArgumentError("Invalid second. It must be an integer between 0 and 60."))

    # Check if the date is valid in terms of number of days in a month.
    if M == 2
        if is_leap_year(Y)
            (D > 29) && throw(ArgumentError("Wrong day number given the year and the month."))
        else
            (D > 28) && throw(ArgumentError("Wrong day number given the year and the month."))
        end
    elseif M in [4, 6, 9, 11]
        (D > 30) && throw(ArgumentError("Wrong day number given the year and the month."))
    end

    # If the month is January / February, then consider it as the 13rd / 14th
    # month of the last year.
    if (M == 1) || (M == 2)
        Y -= 1
        M += 12
    end

    a = div(Y,100)
    b = div(a,4)
    c = 2-a+b
    e = floor(Int,365.25*(Y+4716))
    f = floor(Int,30.6001*(M+1))

    # Compute the Julian Day considering the time of day.
    #
    # Notice that the algorithm in [2] always return the Julian day at 00:00
    # GMT.
    c+D+e+f-1524.5 + ((h*60 + m)*60 + s)/86400
end

"""
    function DatetoJD(date::Date)

Convert the date `date` to Julian Day.

"""
function DatetoJD(date::Date)
    return DatetoJD(Dates.year(date), Dates.month(date), Dates.day(date),
                    0, 0, 0)
end

"""
    function DatetoJD(dateTime::DateTime)

Convert the date and time `dateTime` to Julian Day.

"""
function DatetoJD(dateTime::DateTime)
    return DatetoJD(Dates.year(dateTime),
                    Dates.month(dateTime),
                    Dates.day(dateTime),
                    Dates.hour(dateTime),
                    Dates.minute(dateTime),
                    Dates.second(dateTime))
end

"""
    function JDtoDate([T,] JD::Number)

Convert a date represented in Julian Day `JD` to Gregorian Calendar. The
optional parameter `T` defines the return type. If `T` is omitted, then it
defaults to `Int`.

# Returns

If `T` is omitted or `Int`, then a tuple with the following data will be
returned:

* Year.
* Month (`1` => **January**, `2` => **February**, ...).
* Day.
* Hour (0 - 24).
* Minute (0 - 59).
* Second (0 - 59).

Notice that if `T` is `Int`, then the seconds field will be Integer. Otherwise,
it will be floating point.

If `T` is `Date`, then it will return the Julia structure `Date`. Notice that
the hours, minutes, and seconds will be neglected because the structure `Date`
does not handle them.

If `T` is `DateTime`, then it will return the Julia structure `DateTime`.

# Remarks

The algorithm was obtained from \\[2] (Accessed on 2018-04-11). In [2], there is
the following warning:

> Note: This method will not give dates accurately on the Gregorian Proleptic
> Calendar, i.e., the calendar you get by extending the Gregorian calendar
> backwards to years earlier than 1582. using the Gregorian leap year rules.
> In particular, the method fails if Y<400.

"""
function JDtoDate(JD::Number)
    Q = JD + 0.5
    Z = floor(Int, Q)
    W = div(Z - 1867216.25,36524.25)
    X = div(W,4)
    A = Z+1+W-X
    B = A+1524
    C = div(B-122.1,365.25)
    D = floor(Int,365.25*C)
    E = div(B-D,30.6001)
    F = floor(Int,30.6001*E)

    # In this case, `dayf` will have the fractional part of the day.
    dayf   = B-D-F+(Q-Z)
    monthf = (E-1 <= 12) ? E-1 : E-13
    yearf  = ( (monthf == 1) || (monthf == 2) ) ? C-4715 : C-4716

    # Get the hour, minute, and second from the day.
    hf = (dayf % 1)*24
    mf = (hf   % 1)*60
    sf = (mf   % 1)*60

    # Transform everything in integers.
    year  = floor(Int, yearf)
    month = floor(Int, monthf)
    day   = floor(Int, dayf)
    h     = floor(Int, hf)
    m     = floor(Int, mf)
    s     = sf

    # Return.
    (year, month, day, h, m, s)
end

function JDtoDate(::Type{Int}, JD::Number)
    (year, month, day, h, m, s) = JDtoDate(JD)

    # If the seconds are large than 59.5, then we must take care of the
    # rounding. In this case, we should advance a minute that can trigger an
    # advance in all the other parameters. The safest way here is to just call
    # the function `JDtoDate` again with an offset.

    if s >= 59.5
        Δs = 60.01 - s
        (year, month, day, h, m, s) = JDtoDate(JD + Δs/86400)
    end

    (year, month, day, h, m, round(Int,s))
end

function JDtoDate(::Type{Date}, JD::Number)
    (Y, M, D, ~, ~, ~) = JDtoDate(JD)
    Date(Y,M,D)
end

JDtoDate(::Type{DateTime}, JD::Number) = julian2datetime(JD)
