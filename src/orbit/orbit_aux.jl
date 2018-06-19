#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Propagate the satellite orbit.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A., McClain, W. D (2013). Fundamentals of Astrodynamics
#       and Applications. Microcosm Press.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export satellite_orbit_compute_f

"""
    function satellite_orbit_compute_f(a::Number, e::Number, i::Number, m::Number, tol::Number = 1e-10)

Compute the true anomaly given the mean anomaly `m`, the semi-major axis `a`,
the eccentricity `e`, and the inclination `i`.

# Args

* `a`: Semi-major axis [m].
* `e`: Eccentricity.
* `i`: Inclination [rad].
* `m`: Mean anomaly [rad].
* `tol`: (OPTIONAL) Tolerance for the Newton-Raphson method (**Default**: 1e-10).

# Returns

The true anomaly [rad].

# Remarks

This functions uses the Newton-Raphson method to compute the true anomaly.

"""
function satellite_orbit_compute_f(a::Number,
                                   e::Number,
                                   i::Number,
                                   m::Number,
                                   tol::Number = 1e-10)

    # Update the eccentric anomaly using the Newton-Raphson method.
    # =============================================================

    # Initial guess.
    u = m

    sin_u = sin(u)
    cos_u = cos(u)

    # Newton-Raphson iterations.
    while ( abs(u - e*sin_u - m) > tol )
        u = u - (u - e*sin_u - m)/(1-e*cos_u)

        sin_u = sin(u)
        cos_u = cos(u)
    end

    # Transform u to the interval [0, 2*pi].
    u = mod(u, 2*pi)

    # Compute the true anomaly.
    half_f = atan( sqrt( (1+e)/(1-e) )*(sin_u/(1+cos_u))  )

    # f/2 and u/2 are in the same quadrant.
    if ( (u/2 > pi/2) && (u/2 <= 1.5*pi) )
        half_f += pi
    elseif ( u/2 > 1.5*pi )
        half_f = 2*pi - half_f
    end

    # Compute the true anomaly in the interval [0, 2*pi].
    f = mod(2*half_f, 2*pi)
end
