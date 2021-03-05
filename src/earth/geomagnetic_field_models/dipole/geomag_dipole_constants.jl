#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Constants for the dipole model for the Earth geomagnetic field.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] http://wdc.kugi.kyoto-u.ac.jp/poles/polesexp.html
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# Location of the South pole (which lies in the North hemisphere) along the
# years [°]. This data was obtained from [1].
const _dipole_pole = [1900    +78.7   -68.8
                      1905    +78.7   -68.7
                      1910    +78.7   -68.7
                      1915    +78.6   -68.6
                      1920    +78.6   -68.4
                      1925    +78.6   -68.3
                      1930    +78.6   -68.3
                      1935    +78.6   -68.4
                      1940    +78.5   -68.5
                      1945    +78.5   -68.5
                      1950    +78.5   -68.8
                      1955    +78.5   -69.2
                      1960    +78.6   -69.5
                      1965    +78.6   -69.9
                      1970    +78.7   -70.2
                      1975    +78.8   -70.5
                      1980    +78.9   -70.8
                      1985    +79.0   -70.9
                      1990    +79.2   -71.1
                      1995    +79.4   -71.4
                      2000    +79.6   -71.6
                      2005    +79.8   -71.8
                      2010    +80.1   -72.2
                      2011    +80.1   -72.3
                      2012    +80.2   -72.4
                      2013    +80.3   -72.5
                      2014    +80.3   -72.5
                      2015    +80.4   -72.6
                      2016    +80.4   -72.7
                      2017    +80.5   -72.8
                      2018    +80.5   -73.0
                      2019    +80.6   -73.1
                      2020    +80.6   -73.2]

# Dipole momentum along the years [10²² A.m²]. This data was obtained from [1].
const _dipole_mag = [1900 8.32
                     1905 8.30
                     1910 8.27
                     1915 8.24
                     1920 8.20
                     1925 8.16
                     1930 8.13
                     1935 8.11
                     1940 8.09
                     1945 8.08
                     1950 8.06
                     1955 8.05
                     1960 8.03
                     1965 8.00
                     1970 7.97
                     1975 7.94
                     1980 7.91
                     1985 7.87
                     1990 7.84
                     1995 7.81
                     2000 7.79
                     2005 7.77
                     2010 7.75
                     2011 7.74
                     2012 7.74
                     2013 7.73
                     2014 7.73
                     2015 7.72
                     2016 7.72
                     2017 7.72
                     2018 7.71
                     2019 7.71
                     2020 7.70]

# Interpolations.
const _itp_dipole_lat = extrapolate(
    interpolate( (_dipole_pole[:,1],), _dipole_pole[:,2], Gridded(Linear())),
    Flat())
const _itp_dipole_lon = extrapolate(
    interpolate( (_dipole_pole[:,1],), _dipole_pole[:,3], Gridded(Linear())),
    Flat())
const _itp_dipole_mag = extrapolate(
    interpolate( (_dipole_mag[:,1],), _dipole_mag[:,2], Gridded(Linear())),
    Flat())
