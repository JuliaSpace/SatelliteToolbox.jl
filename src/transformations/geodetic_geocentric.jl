# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
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
#   [5] Borkowski, K. M (1987). Transformation of geocentric to geodetic
#       coordinates without approximations. Astrophysics and Space Science, vol.
#       139, pp. 1-4.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export ecef_to_geodetic, geodetic_to_ecef
export geocentric_to_geodetic, GeodetictoGeocentric

"""
    ecef_to_geodetic(r_e::AbstractVector)

Convert the vector `r_e` [m] represented in the Earth-Centered, Earth-Fixed
(ECEF) reference frame into Geodetic coordinates (WGS-84).

# Returns

* Latitude [rad].
* Longitude [rad].
* Altitude [m].

# Remarks

Based on algorithm in [3].

"""
function ecef_to_geodetic(r_e::AbstractVector)
    # Auxiliary variables.
    X = r_e[1]
    Y = r_e[2]
    Z = r_e[3]

    # Auxiliary variables.
    p = √(X^2 + Y^2)
    θ = atan(Z * a_wgs84, p * b_wgs84)
    e_wgs84² = e_wgs84^2
    sin_θ, cos_θ = sincos(θ)

    # Compute Geodetic.
    lon = atan(Y, X)
    lat = atan(Z + el_wgs84^2 * b_wgs84 * sin_θ^3,
               p -   e_wgs84² * a_wgs84 * cos_θ^3)

    sin_lat, cos_lat = sincos(lat)

    N = a_wgs84/√(1 - e_wgs84² * sin_lat^2)

    # Avoid singularity if we are near the poles (~ 1 deg according to [1,
    # p.172]). Note that `cosd(1) = -0.01745240643728351`.
    if !(-0.01745240643728351 < cos_lat < 0.01745240643728351)
        h = p/cos_lat - N
    else
        h = Z/sin_lat - N*(1 - e_wgs84²)
    end

    return lat, lon, h
end

"""
    geodetic_to_ecef(lat::Number, lon::Number, h::Number)

Convert the latitude `lat` [rad], longitude `lon` [rad], and altitude `h` \\[m]
(WGS-84) into a vector represented on the Earth-Centered, Earth-Fixed (ECEF)
reference frame.

# Remarks

Based on algorithm in [3].

"""
function geodetic_to_ecef(lat::Number, lon::Number, h::Number)
    # Auxiliary variables.
    sin_lat, cos_lat = sincos(lat)
    sin_lon, cos_lon = sincos(lon)

    # Radius of curvature [m].
    N = a_wgs84/√(1 - e_wgs84^2 * sin_lat^2)

    # Compute the position in ECEF frame.
    return SVector(
        (                        N + h) * cos_lat * cos_lon,
        (                        N + h) * cos_lat * sin_lon,
        ((b_wgs84 / a_wgs84)^2 * N + h) * sin_lat
    )
end

"""
    geocentric_to_geodetic(ϕ_gc::Number, r::Number)

Compute the geodetic latitude and altitude (WGS-84) from the geocentric latitude
`ϕ_gc` (-π/2, π/2) [rad] and radius `r` [m]. Notice that the longitude is the
same in both geocentric and geodetic coordinates.

# Returns

* Geodetic latitude [rad].
* Altitude above the reference ellipsoid (WGS-84) [m].

# Remarks

Based on algorithm in [5].

"""
function geocentric_to_geodetic(ϕ_gc::Number, r::Number)
    # Obtain the `z` component and the equatorial component `re`.
    sin_ϕ_gc, cos_ϕ_gc = sincos(ϕ_gc)
    re = r * cos_ϕ_gc
    z  = r * sin_ϕ_gc

    # Auxiliary variables.
    a  = a_wgs84
    a² = a^2
    b  = b_wgs84
    b² = b^2

    # Compute the parameters.
    E  = (b * z - (a² - b²)) / (a * re)
    E² = E^2
    F  = (b * z + (a² - b²)) / (a * re)
    P  = (4 / 3) * (E * F + 1)
    Q  = 2 * (E² - F^2)
    D  = P^3 + Q^2

    if D ≥ 0
        aux = √D
        v = (aux - Q)^(1/3) - (aux + Q)^(1/3)
    else
        v = 2 * √(-P) * cos( acos(Q / (-P)^(3 / 2)) / 3)
    end

    # We must select the correct root depending on where we are (North or South
    # hemisphere) to avoid complex results.
    G = (sign(z) * √(E² + v) + E)/2
    t = sign(z) * √(G^2 + (F - v * G) / (2 * G - E)) - G

    # Compute the geodetic latitude and altitude.
    ϕ_gd = atan(a * (1 - t^2), 2b * t)
    sin_ϕ_gd, cos_ϕ_gd = sincos(ϕ_gd)
    h = (re - a * t) * cos_ϕ_gd + (z - b) * sin_ϕ_gd

    return ϕ_gd, h
end

"""
    GeodetictoGeocentric(ϕ_gd::Number, h::Number)

Compute the geocentric latitude and radius from the geodetic latitude `ϕ_gd`
(-π/2, π/2) [rad] and height above the reference ellipsoid `h` \\[m] (WGS-84).
Notice that the longitude is the same in both geocentric and geodetic
coordinates.

# Returns

* Geocentric latitude [rad].
* Radius from the center of the Earth [m].

# Remarks

Based on algorithm in [4, p. 3].

"""
function GeodetictoGeocentric(ϕ_gd::Number, h::Number)
    # Auxiliary variables to decrease computational burden.
    sin_ϕ_gd, cos_ϕ_gd = sincos(ϕ_gd)
    sin²_ϕ_gd = sin_ϕ_gd^2
    e²_wgs84  = e_wgs84^2

    # Radius of curvature in the prime vertical [m].
    N    = a_wgs84/√(1 - e²_wgs84 * sin²_ϕ_gd )

    # Compute the geocentric latitude and radius from the Earth center.
    ρ    = (N + h) * cos_ϕ_gd
    z    = (N * (1 - e²_wgs84) + h) * sin_ϕ_gd
    r    = √(ρ^2 + z^2)
    ϕ_gc = asin(z / r)

    return ϕ_gc, r
end
