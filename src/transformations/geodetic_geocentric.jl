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
#    Coordinate transformations related with the geodetic and geocentric
#    coordinates.
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
#   [4] ISO TC 20/SC 14 N (2011). Geomagnetic Reference Models.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-05-13: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   The `LLA` in function names was changed to `Geodetic`.
#
# 2018-05-07: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Add function to convert from geodetic to geocentric latitude.
#
# 2018-03-27: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   The functions related to the GMST were moved to another file to improve code
#   organization.
#
# 2016-03-10: Ronan Arraes Jardim Chagsa <ronan.arraes@inpe.br>
#   Add code to convert from Geodetic (WGS-84) to ECEF reference frame.
#
# 2015-11-23: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export ECEFtoGeodetic, GeodetictoECEF

export GeodetictoGeocentric

"""
### function ECEFtoGeodetic(r_e::Vector)

Convert the vector `r_e` represented in the Earth-Centered, Earth-Fixed (ECEF)
reference frame into Geodetic coordinates (WGS-84).

##### Args

* r_e: Position vector represented in the ECEF reference frame.

##### Returns

* Latitude [rad].
* Longitude [rad].
* Altitude [m].

##### Remarks

Based on algorithm in [3].

"""

function ECEFtoGeodetic(r_e::Vector)
    # Auxiliary variables.
    X = r_e[1]
    Y = r_e[2]
    Z = r_e[3]

    p = sqrt(X^2 + Y^2)
    θ = atan2(Z*a_wgs84, p*b_wgs84)

    # Compute Geodetic.
    lon = atan2(Y, X)
    lat = atan2(Z + el_wgs84^2*b_wgs84*sin(θ)^3,
                p -  e_wgs84^2*a_wgs84*cos(θ)^3)
    N   = a_wgs84/sqrt(1 - e_wgs84^2*sin(lat)^2 )
    h   = p/cos(lat) - N

    (lat, lon, h)
end

"""
### function GeodetictoECEF(lat::Number, lon::Number, h::Number)

Convert the latitude `lat`, longitude `lon`, and altitude `h` (WGS-84) into the
Earth-Centered, Earth-Fixed (ECEF) reference frame.

##### Args

* lat: Latitude in WGS-84 [rad].
* lon: Longitude in WGS-84 [rad].
* h: Altitude in WGS-84 [rad].

##### Returns

The converted vector represented in the ECEF reference frame.

##### Remarks

Based on algorithm in [3].

"""

function GeodetictoECEF(lat::Number, lon::Number, h::Number)
    # Radius of curvature [m].
    N = a_wgs84/sqrt(1 - e_wgs84^2*sin(lat)^2 )

    # Compute the position in ECEF frame.
    [ (                     N + h)*cos(lat)*cos(lon);
      (                     N + h)*cos(lat)*sin(lon);
      ( (b_wgs84/a_wgs84)^2*N + h)*sin(lat);]
end

"""
### function GeodetictoGeocentric(ϕ_gd::Number, h::Number)

Compute the geocentric latitude and radius from the geodetic latitude `ϕ_gd` and
height above the reference ellipsoid `h` (WGS-84). Notice that the longitude is
the same in both geocentric and geodetic coordinates.

##### Args

* ϕ_gd: Geodetic latitude (-π/2,+π/2) [rad].
* h: Height above the reference ellipsoid (WGS-84) [m].

##### Returns

* Geocentric latitude [rad].
* Radius from the center of the Earth [m].

##### Remarks

Based on algorithm in [4, p. 3].

"""

function GeodetictoGeocentric(ϕ_gd::Number, h::Number)
    # Auxiliary variables to decrease computational burden.
    sin_ϕ_gd  = sin(ϕ_gd)
    sin2_ϕ_gd = sin_ϕ_gd^2
    e2_wgs84  = e_wgs84^2

    # Radius of curvature in the prime vertical [m].
    N    = a_wgs84/sqrt(1 - e2_wgs84*sin2_ϕ_gd )

    # Compute the geocentric latitude and radius from the Earth center.
    ρ    = (N + h)*cos(ϕ_gd)
    z    = (N * (1 - e2_wgs84) + h)*sin_ϕ_gd
    r    = sqrt(ρ^2 + z^2)
    ϕ_gc = asin(z/r)

    (ϕ_gc, r)
end
