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
# ==============================================================================
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
export geocentric_to_geodetic, geodetic_to_geocentric

"""
    ecef_to_geodetic(r_e::AbstractVector; ellipsoid=wgs84_ellipsoid)

Convert the vector `r_e` [m] represented in the Earth-Centered, Earth-Fixed
(ECEF) reference frame into Geodetic coordinates for a custom target ellipsoid
(defaults to WGS-84).

!!! info
    The algorithm is based in **[1]**.

# Returns

- Latitude [rad].
- Longitude [rad].
- Altitude [m].

# Reference

- **[1]**: mu-blox ag (1999). Datum Transformations of GPS Positions.
    Application Note.
"""
function ecef_to_geodetic(r_e::AbstractVector; ellipsoid = wgs84_ellipsoid)
    # Auxiliary variables.
    x = r_e[1]
    y = r_e[2]
    z = r_e[3]

    # Auxiliary variables.
    a = ellipsoid.a
    b = ellipsoid.b
    e² = ellipsoid.e²
    el² = ellipsoid.el²
    p = √(x^2 + y^2)
    θ = atan(z * a, p * b)
    sin_θ, cos_θ = sincos(θ)

    # Compute Geodetic.
    lon = atan(y, x)
    lat = atan(
        z + el² * b * sin_θ^3,
        p -  e² * a * cos_θ^3
    )

    sin_lat, cos_lat = sincos(lat)

    N = a/√(1 - e² * sin_lat^2)

    # Avoid singularity if we are near the poles (~ 1 deg according to [1,
    # p.172]). Note that `cosd(1) = -0.01745240643728351`.
    if !(-0.01745240643728351 < cos_lat < 0.01745240643728351)
        h = p / cos_lat - N
    else
        h = z / sin_lat - N * (1 - e²)
    end

    return lat, lon, h
end

"""
    geodetic_to_ecef(lat::Number, lon::Number, h::Number; ellipsoid = wgs84_ellipsoid)

Convert the latitude `lat` [rad], longitude `lon` [rad], and altitude `h` \\[m] above the
reference ellipsoid (defaults to WGS-84) into a vector represented on the Earth-Centered,
Earth-Fixed (ECEF) reference frame.

!!! info
    The algorithm is based in **[1]**.

# Reference

- **[1]**: mu-blox ag (1999). Datum Transformations of GPS Positions.
    Application Note.
"""
function geodetic_to_ecef(lat::Number, lon::Number, h::Number; ellipsoid = wgs84_ellipsoid)
    # Auxiliary variables.
    sin_lat, cos_lat = sincos(lat)
    sin_lon, cos_lon = sincos(lon)
    a = ellipsoid.a
    b = ellipsoid.b
    e² = ellipsoid.e²

    # Radius of curvature [m].
    N = a/√(1 - e² * sin_lat^2)

    # Compute the position in ECEF frame.
    return SVector(
        (            N + h) * cos_lat * cos_lon,
        (            N + h) * cos_lat * sin_lon,
        ((b / a)^2 * N + h) * sin_lat
    )
end

"""
    geocentric_to_geodetic(ϕ_gc::Number, r::Number; ellipsoid = wgs84_ellipsoid)

Compute the geodetic latitude and altitude above the reference ellipsoid (defaults to
WGS-84) from the geocentric latitude `ϕ_gc` (-π/2, π/2) [rad] and radius `r` [m].
Notice that the longitude is the same in both geocentric and geodetic coordinates.

!!! info
    The algorithm is based in **[1]**.

# Returns

* Geodetic latitude [rad].
* Altitude above the reference ellipsoid (defaults to WGS-84) [m].

# References

- **[1]** Borkowski, K. M (1987). Transformation of geocentric to geodetic
    coordinates without approximations. Astrophysics and Space Science, vol.
    139, pp. 1-4.
"""
function geocentric_to_geodetic(ϕ_gc::Number, r::Number; ellipsoid = wgs84_ellipsoid)
    # Obtain the `z` component and the equatorial component `re`.
    sin_ϕ_gc, cos_ϕ_gc = sincos(ϕ_gc)
    re = r * cos_ϕ_gc
    z  = r * sin_ϕ_gc
    sign_z = z == 0 ? +1 : sign(z)

    # Auxiliary variables.
    a  = ellipsoid.a
    a² = a^2
    b  = sign_z * ellipsoid.b
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
        v = (aux - Q)^(1/3) - (Q + aux)^(1/3)
    else
        aux = √(-P)
        v = 2 * aux * cos(acos(Q / (P * aux)) / 3)
    end

    G = (√(E² + v) + E)/2

    # NOTE: Reference [5] appears to have an error in Eq. (13), where we must
    # have G^2 instead of G inside the square root. The correct version can be
    # seen here:
    #
    #   https://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm
    #
    t = √(G^2 + (F - v * G) / (2 * G - E)) - G

    # Compute the geodetic latitude and altitude.
    ϕ_gd = atan(a * (1 - t^2)/(2b * t))
    sin_ϕ_gd, cos_ϕ_gd = sincos(ϕ_gd)
    h = (re - a * t) * cos_ϕ_gd + (z - b) * sin_ϕ_gd

    return ϕ_gd, h
end

"""
    geodetic_to_geocentric(ϕ_gd::Number, h::Number; ellipsoid = wgs84_ellipsoid)

Compute the geocentric latitude and radius from the geodetic latitude `ϕ_gd`
(-π/2, π/2) [rad] and height above the reference ellipsoid `h` \\[m] (defaults to WGS-84).
Notice that the longitude is the same in both geocentric and geodetic
coordinates.

!!! info
    The algorithm is based in **[1]**(p. 3).

# Returns

* Geocentric latitude [rad].
* Radius from the center of the Earth [m].

# References

- **[1]** ISO TC 20/SC 14 N (2011). Geomagnetic Reference Models.
"""
function geodetic_to_geocentric(ϕ_gd::Number, h::Number; ellipsoid = wgs84_ellipsoid)
    # Auxiliary variables to decrease computational burden.
    sin_ϕ_gd, cos_ϕ_gd = sincos(ϕ_gd)
    sin²_ϕ_gd = sin_ϕ_gd^2
    a = ellipsoid.a
    e² = ellipsoid.e²

    # Radius of curvature in the prime vertical [m].
    N = a / √(1 - e² * sin²_ϕ_gd )

    # Compute the geocentric latitude and radius from the Earth center.
    ρ    = (N + h) * cos_ϕ_gd
    z    = (N * (1 - e²) + h) * sin_ϕ_gd
    r    = √(ρ^2 + z^2)
    ϕ_gc = asin(z / r)

    return ϕ_gc, r
end
