## Constants
const PI = 3.14159265358979
const PI2 = 6.28318530717959
const PIM = 1.57079632679490

export sun_position_el

## type definition of struct
struct sunpos
    # Inputs
    UT::Float64
    Day::Int64
    Month::Int64
    Year::Int64
    Dt::Float64
    Longitude::Float64
    Latitude::Float64
    Pressure::Float64
    Temperature::Float64
    
    # Outputs
    RightAscension::Union{Missing, Float64}
    Declination::Union{Missing, Float64}
    HourAngle::Union{Missing, Float64}
    Zenith::Union{Missing, Float64}
    Azimuth::Union{Missing, Float64}
end

## Outer Constructor - Multiple Dispatch
function sunpos(UT::Float64=0.0, Day::Int64=1, Month::Int64=1, Year::Int64=2010,
    Dt::Union{Bool, Float64}=65.0, Longitude::Float64=0.0, Latitude::Float64=0.0,
    Pressure::Float64=1.0, Temperature::Float64=20.0)

    # Outputs
    RightAscension = missing;
    Declination = missing;
    HourAngle = missing;
    Zenith = missing;
    Azimuth = missing;
    
    if Dt == true
        Dt = diff_tt_utc(Year);
    end

    return sunpos(UT, Day, Month, Year, Dt, Longitude, Latitude, Pressure, Temperature, 
    RightAscension, Declination, HourAngle, Zenith, Azimuth);
end

function sunpos(JD::Number, Dt::Union{Bool, Float64}=65.0, Longitude::Float64=0.0, 
    Latitude::Float64=0.0, Pressure::Float64=1.0, Temperature::Float64=20.0)

    Year, Month, Day, h, m, s = jd_to_date(JD);
    UT = hms_to_h(h, m, s);

    return sunpos(UT, Day, Month, Year, Dt, Longitude, Latitude, Pressure, Temperature);
end

## Algorithms
#= 
The Julia code was transcribed from C++ obtained at
http://www.solaritaly.enea.it/StrSunPosition/SunPositionEn.php

by Stephan Buchert, Swedish Institute of Space Physics, scb@irfu.se

source: https://ui.adsabs.harvard.edu/abs/2012SoEn...86.1323G/abstract

The Julia code was transcribed from C++ obtained at
http://www.solaritaly.enea.it/StrSunPosition/SunPositionEn.php

by Stephan Buchert, Swedish Institute of Space Physics, scb@irfu.se

The maximum errors of the algorithms (in the long versions), in absolute angular
distance from the true position, obtained on 20 million trials, are:
Algorithm 1: 0.19 deg
Algorithm 2: 0.034 deg
Algorithm 3: 0.0094 deg
Algorithm 4: 0.0094 deg
Algorithm 5: 0.0027 deg.

Inputs:
UT: Universal time UT (Greenwich time), in hours;
Day, Month, and Year;
Dt: Difference TT-UT (terrestrial time - universal time), in seconds;
Longitude and Latitude of the observer, in radians;
Pressure in atm;
Temperature in Celsius degrees.

=#

function sunpos_Algorithm1(s::sunpos)
    # implementation of Algorithm1

    # Inputs
    UT = s.UT;
    Day = s.Day;
    Month = s.Month;
    Year = s.Year;
    Dt = s.Dt;
    Longitude = s.Longitude;
    Latitude = s.Latitude;
    Pressure = s.Pressure;
    Temperature = s.Temperature;

    # Algorithm1
    if s.Month <= 2
        mt = Month + 12;
        yt = Year - 1;
    else
        mt = Month;
        yt = Year;
    end
    
    t = convert(Float64 , floor(365.25*(yt-2000)) + floor(30.6001*(mt+1)) - floor(0.01*(yt)) + Day ) + 0.0416667*UT - 21958.0;
    te = t + 1.1574e-5*Dt;

    wte = 0.017202786*te;

    s1 = sin(wte);
    c1 = cos(wte);
    s2 = 2.0*s1*c1;
    c2 = (c1+s1)*(c1-s1);

    RightAscension = -1.38880 + 1.72027920e-2*te + 3.199e-2*s1 - 2.65e-3*c1 + 4.050e-2*s2 + 1.525e-2*c2;
    RightAscension = rem(RightAscension, PI2);

    Declination = 6.57e-3 + 7.347e-2*s1 - 3.9919e-1*c1 + 7.3e-4*s2 - 6.60e-3*c2;

    HourAngle = 1.75283 + 6.3003881*t + Longitude - RightAscension;
    HourAngle = rem(HourAngle + PI, PI2) - PI;

    # Outputs
    sp = sin(Latitude);
    cp = sqrt((1-sp*sp));
    sd = sin(Declination);
    cd = sqrt(1-sd*sd);
    sH = sin(HourAngle);
    cH = cos(HourAngle);
    se0 = sp*sd + cp*cd*cH;
    ep = asin(se0) - 4.26e-5*sqrt(1.0-se0*se0);

    if (ep > 0.0)
        De = (0.08422*Pressure) / ((273.0+Temperature)*tan(ep + 0.003138/(ep + 0.08919)));
    else
        De = 0.0;
    end

    Zenith = PIM - ep - De;
    Azimuth = atan(sH, cH*sp - sd*cp/cd);

    snew = sunpos(UT, Day, Month, Year, Dt, Longitude, Latitude, Pressure, Temperature,
    RightAscension, Declination, HourAngle, Zenith, Azimuth);

    return snew;
end

function sunpos_Algorithm2(s::sunpos)
    # implementation of Algorithm2

    # Inputs
    UT = s.UT;
    Day = s.Day;
    Month = s.Month;
    Year = s.Year;
    Dt = s.Dt;
    Longitude = s.Longitude;
    Latitude = s.Latitude;
    Pressure = s.Pressure;
    Temperature = s.Temperature;

    # Algorithm2
    if s.Month <= 2
        mt = Month + 12;
        yt = Year - 1;
    else
        mt = Month;
        yt = Year;
    end
    
    t = convert(Float64 , floor(365.25*(yt-2000)) + floor(30.6001*(mt+1)) - floor(0.01*(yt)) + Day ) + 0.0416667*UT - 21958.0;
    te = t + 1.1574e-5*Dt;

    wte = 0.017202786*te;

    s1 = sin(wte);
    c1 = cos(wte);
    s2 = 2.0*s1*c1;
    c2 = (c1+s1)*(c1-s1);
    s3 = s2*c1 + c2*s1;
    c3 = c2*c1 - s2*s1;
    s4 = 2.0*s2*c2;
    c4 = (c2+s2)*(c2-s2);

    RightAscension = -1.38880 + 1.72027920e-2*te + 3.199e-2*s1 - 2.65e-3*c1 + 4.050e-2*s2 + 1.525e-2*c2 + 1.33e-3*s3 + 3.8e-4*c3 + 7.3e-4*s4 + 6.2e-4*c4;
    RightAscension = rem(RightAscension, PI2);
  
    Declination = 6.57e-3 + 7.347e-2*s1 - 3.9919e-1*c1 + 7.3e-4*s2 - 6.60e-3*c2 + 1.50e-3*s3 - 2.58e-3*c3 + 6e-5*s4 - 1.3e-4*c4;
  
    HourAngle = 1.75283 + 6.3003881*t + Longitude - RightAscension;
    HourAngle = rem(HourAngle + PI, PI2) - PI;

    # Outputs
    sp = sin(Latitude);
    cp = sqrt((1-sp*sp));
    sd = sin(Declination);
    cd = sqrt(1-sd*sd);
    sH = sin(HourAngle);
    cH = cos(HourAngle);
    se0 = sp*sd + cp*cd*cH;
    ep = asin(se0) - 4.26e-5*sqrt(1.0-se0*se0);

    if (ep > 0.0)
        De = (0.08422*Pressure) / ((273.0+Temperature)*tan(ep + 0.003138/(ep + 0.08919)));
    else
        De = 0.0;
    end

    Zenith = PIM - ep - De;
    Azimuth = atan(sH, cH*sp - sd*cp/cd);

    snew = sunpos(UT, Day, Month, Year, Dt, Longitude, Latitude, Pressure, Temperature,
    RightAscension, Declination, HourAngle, Zenith, Azimuth);

    return snew;
end

function sunpos_Algorithm3(s::sunpos)
    # implementation of Algorithm3
    
    # Inputs
    UT = s.UT;
    Day = s.Day;
    Month = s.Month;
    Year = s.Year;
    Dt = s.Dt;
    Longitude = s.Longitude;
    Latitude = s.Latitude;
    Pressure = s.Pressure;
    Temperature = s.Temperature;

    # Algorithm3
    if s.Month <= 2
        mt = Month + 12;
        yt = Year - 1;
    else
        mt = Month;
        yt = Year;
    end
    
    t = convert(Float64 , floor(365.25*(yt-2000)) + floor(30.6001*(mt+1)) - floor(0.01*(yt)) + Day ) + 0.0416667*UT - 21958.0;
    te = t + 1.1574e-5*Dt;

    wte = 0.0172019715*te;

    lambda = -1.388803 + 1.720279216e-2*te + 3.3366e-2*sin(wte - 0.06172) + 3.53e-4*sin(2.0*wte - 0.1163);

    epsi = 4.089567e-1 - 6.19e-9*te;

    sl = sin(lambda);
    cl = cos(lambda);
    se = sin(epsi);
    ce = sqrt(1-se*se);

    RightAscension = atan(sl*ce, cl);
    if (RightAscension < 0.0) 
      RightAscension += PI2;
    end
  
    Declination = asin(sl*se);
  
    HourAngle = 1.7528311 + 6.300388099*t + Longitude - RightAscension;
    HourAngle = rem(HourAngle + PI, PI2) - PI;

    # Outputs
    sp = sin(Latitude);
    cp = sqrt((1-sp*sp));
    sd = sin(Declination);
    cd = sqrt(1-sd*sd);
    sH = sin(HourAngle);
    cH = cos(HourAngle);
    se0 = sp*sd + cp*cd*cH;
    ep = asin(se0) - 4.26e-5*sqrt(1.0-se0*se0);

    if (ep > 0.0)
        De = (0.08422*Pressure) / ((273.0+Temperature)*tan(ep + 0.003138/(ep + 0.08919)));
    else
        De = 0.0;
    end

    Zenith = PIM - ep - De;
    Azimuth = atan(sH, cH*sp - sd*cp/cd);

    snew = sunpos(UT, Day, Month, Year, Dt, Longitude, Latitude, Pressure, Temperature,
    RightAscension, Declination, HourAngle, Zenith, Azimuth);

    return snew;
end

function sunpos_Algorithm4(s::sunpos)
    # implementation of Algorithm4
    
    # Inputs
    UT = s.UT;
    Day = s.Day;
    Month = s.Month;
    Year = s.Year;
    Dt = s.Dt;
    Longitude = s.Longitude;
    Latitude = s.Latitude;
    Pressure = s.Pressure;
    Temperature = s.Temperature;

    # Algorithm4
    if s.Month <= 2
        mt = Month + 12;
        yt = Year - 1;
    else
        mt = Month;
        yt = Year;
    end
    
    t = convert(Float64 , floor(365.25*(yt-2000)) + floor(30.6001*(mt+1)) - floor(0.01*(yt)) + Day ) + 0.0416667*UT - 21958.0;
    te = t + 1.1574e-5*Dt;

    wte = 0.0172019715*te;

    L = 1.752790 + 1.720279216e-2*te + 3.3366e-2*sin(wte - 0.06172) + 3.53e-4*sin(2.0*wte - 0.1163);

    nu = 9.282e-4*te - 0.8;
    Dlam = 8.34e-5*sin(nu);
    lambda = L + PI + Dlam;

    epsi = 4.089567e-1 - 6.19e-9*te + 4.46e-5*cos(nu);

    sl = sin(lambda);
    cl = cos(lambda);
    se = sin(epsi);
    ce = sqrt(1-se*se);

    RightAscension = atan(sl*ce, cl);
    if (RightAscension < 0.0) 
        RightAscension += PI2;
    end

    Declination = asin(sl*se);

    HourAngle = 1.7528311 + 6.300388099*t + Longitude - RightAscension + 0.92*Dlam;
    HourAngle = rem(HourAngle + PI, PI2) - PI;

    # Outputs
    sp = sin(Latitude);
    cp = sqrt((1-sp*sp));
    sd = sin(Declination);
    cd = sqrt(1-sd*sd);
    sH = sin(HourAngle);
    cH = cos(HourAngle);
    se0 = sp*sd + cp*cd*cH;
    ep = asin(se0) - 4.26e-5*sqrt(1.0-se0*se0);

    if (ep > 0.0)
        De = (0.08422*Pressure) / ((273.0+Temperature)*tan(ep + 0.003138/(ep + 0.08919)));
    else
        De = 0.0;
    end

    Zenith = PIM - ep - De;
    Azimuth = atan(sH, cH*sp - sd*cp/cd);

    snew = sunpos(UT, Day, Month, Year, Dt, Longitude, Latitude, Pressure, Temperature,
    RightAscension, Declination, HourAngle, Zenith, Azimuth);

    return snew;
end

function sunpos_Algorithm5(s::sunpos)
    # implementation of Algorithm5
    
    # Inputs
    UT = s.UT;
    Day = s.Day;
    Month = s.Month;
    Year = s.Year;
    Dt = s.Dt;
    Longitude = s.Longitude;
    Latitude = s.Latitude;
    Pressure = s.Pressure;
    Temperature = s.Temperature;

    # Algorithm5
    if s.Month <= 2
        mt = Month + 12;
        yt = Year - 1;
    else
        mt = Month;
        yt = Year;
    end
    
    t = convert(Float64 , floor(365.25*(yt-2000)) + floor(30.6001*(mt+1)) - floor(0.01*(yt)) + Day ) + 0.0416667*UT - 21958.0;
    te = t + 1.1574e-5*Dt;

    wte = 0.0172019715*te;

    s1 = sin(wte);
    c1 = cos(wte);
    s2 = 2.0*s1*c1;
    c2 = (c1+s1)*(c1-s1);
    s3 = s2*c1 + c2*s1;
    c3 = c2*c1 - s2*s1;

    L = 1.7527901 + 1.7202792159e-2*te + 3.33024e-2*s1 - 2.0582e-3*c1 + 3.512e-4*s2 - 4.07e-5*c2 + 5.2e-6*s3 - 9e-7*c3 
        -8.23e-5*s1*sin(2.92e-5*te) + 1.27e-5*sin(1.49e-3*te - 2.337) + 1.21e-5*sin(4.31e-3*te + 3.065) 
        + 2.33e-5*sin(1.076e-2*te - 1.533) + 3.49e-5*sin(1.575e-2*te - 2.358) + 2.67e-5*sin(2.152e-2*te + 0.074) 
        + 1.28e-5*sin(3.152e-2*te + 1.547) + 3.14e-5*sin(2.1277e-1*te - 0.488);

    nu = 9.282e-4*te - 0.8;
    Dlam = 8.34e-5*sin(nu);
    lambda = L + PI + Dlam;

    epsi = 4.089567e-1 - 6.19e-9*te + 4.46e-5*cos(nu);

    sl = sin(lambda);
    cl = cos(lambda);
    se = sin(epsi);
    ce = sqrt(1-se*se);

    RightAscension = atan(sl*ce, cl);
    if (RightAscension < 0.0) 
        RightAscension += PI2;
    end

    Declination = asin(sl*se);

    HourAngle = 1.7528311 + 6.300388099*t + Longitude - RightAscension + 0.92*Dlam;
    HourAngle = rem(HourAngle + PI, PI2) - PI;

    # Outputs
    sp = sin(Latitude);
    cp = sqrt((1-sp*sp));
    sd = sin(Declination);
    cd = sqrt(1-sd*sd);
    sH = sin(HourAngle);
    cH = cos(HourAngle);
    se0 = sp*sd + cp*cd*cH;
    ep = asin(se0) - 4.26e-5*sqrt(1.0-se0*se0);

    if (ep > 0.0)
        De = (0.08422*Pressure) / ((273.0+Temperature)*tan(ep + 0.003138/(ep + 0.08919)));
    else
        De = 0.0;
    end

    Zenith = PIM - ep - De;
    Azimuth = atan(sH, cH*sp - sd*cp/cd);

    snew = sunpos(UT, Day, Month, Year, Dt, Longitude, Latitude, Pressure, Temperature,
    RightAscension, Declination, HourAngle, Zenith, Azimuth);

    return snew;
end

# Main Function
"""
    sun_position_el(JD::Number)

Compute the Sun position represented in the Local Horizon Reference System at the
Julian Day `JD`, Longitude `Longitude`, Latitude `Latitude`, Atmospheric pressure `Pressure`, 
Ambient Temperature `Temperature`, Output Flag `flag` and Algorithm `algorithm`.

Inputs:

JD: Julian Day;
Longitude and Latitude of the observer, in degrees, WSG84;
Pressure in atm;
Temperature in Celsius degrees.


Outputs:

Equatorial system: flag -e
    RightAscension, in radians;
    Declination, in radians;

Local Coordinates: flag -l
    HourAngle, in radians;
    Zenith, in radians;
    Azimuth, in radians;

flag -a: all outputs

"""
function sun_position_el(JD::Number, Longitude::Number=0.0, Latitude::Number=0.0, 
    Pressure::Number=1.0, Temperature::Number=20.0, flag::Char='l', algorithm::Char='5')

    s = sunpos(JD, true, Longitude*PI/180, Latitude*PI/180, Pressure, Temperature);

    # Switch Algorithms
    if algorithm == '1'
        s_new = sunpos_Algorithm1(s); # Use Algorithm1
    elseif algorithm == '2'
        s_new = sunpos_Algorithm2(s); # Use Algorithm2
    elseif algorithm == '3'
        s_new = sunpos_Algorithm3(s); # Use Algorithm3
    elseif algorithm == '4'
        s_new = sunpos_Algorithm4(s); # Use Algorithm4
    elseif algorithm == '5'
        s_new = sunpos_Algorithm5(s); # Use Algorithm5
    end

    # Switch Output
    if flag == 'e' # Equtorial System
        return s_new.RightAscension, s_new.Declination; 
    elseif flag == 'l' # Local Coordinates
        return s_new.HourAngle, s_new.Zenith, s_new.Azimuth;
    elseif flag == 'a' # All 
        return s_new.RightAscension, s_new.Declination, 
        s_new.HourAngle, s_new.Zenith, s_new.Azimuth;
    end
end