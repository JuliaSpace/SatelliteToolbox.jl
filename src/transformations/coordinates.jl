#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    Coordinate transformations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] ESA Navipedia: http://www.navipedia.net/
#
#   [3] mu-blox ag (1999). Datum Transformations of GPS Positions. Application
#       Note.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-03-27: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   The functions related to the GMST were moved to another file to improve code
#   organization.
#
# 2016-03-10: Ronan Arraes Jardim Chagsa <ronan.arraes@inpe.br>
#   Add code to convert from LLA (WGS-84) to ECEF reference frame.
#
# 2015-11-23: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export ECEFtoLLA, LLAtoECEF

"""
### function ECEFtoLLA(r_e::Vector)

Convert the vector `r_e` in the ECEF reference frame into LLA (WGS-84).

##### Args

* r_e: Position vector represented in the ECEF reference frame.

##### Returns

* Latitude [rad].
* Longitude [rad].
* Altitude [m].

##### Remarks

Based on algorithm in [3].

"""

function ECEFtoLLA(r_e::Vector)
    # Auxiliary variables.
    X = r_e[1]
    Y = r_e[2]
    Z = r_e[3]

    p = sqrt(X^2 + Y^2)
    θ = atan2(Z*a_wgs84, p*b_wgs84)

    # Compute LLA.
    lon = atan2(Y, X)
    lat = atan2(Z + el_wgs84^2*b_wgs84*sin(θ)^3,
                p -  e_wgs84^2*a_wgs84*cos(θ)^3)
    N   = a_wgs84/sqrt(1 - e_wgs84^2*sin(lat)^2 )
    h   = p/cos(lat) - N

    (lat, lon, h)
end

"""
### function LLAtoECEF(lat::Number, lon::Number, h::Number)

Convert the latitude `lat`, longitude `lon`, and altitude `h` (WGS-84) into the
ECEF reference frame.

##### Args

* lat: Latitude in WGS-84 [rad].
* lon: Longitude in WGS-84 [rad].
* h: Altitude in WGS-84 [rad].

##### Returns

The converted vector represented in the ECEF reference frame.

##### Remarks

Based on algorithm in [3].

"""

function LLAtoECEF(lat::Number, lon::Number, h::Number)
    # Radius of curvature [m].
    N = a_wgs84/sqrt(1 - e_wgs84^2*sin(lat)^2 )

    # Compute the position in ECEF frame.
    [ (                     N + h)*cos(lat)*cos(lon);
      (                     N + h)*cos(lat)*sin(lon);
      ( (b_wgs84/a_wgs84)^2*N + h)*sin(lat);]
end
