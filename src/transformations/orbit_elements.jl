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
#   Conversions related to the orbit elements.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Schwarz, R (2014). Memorandum No. 2: Cartesian State Vectors to
#       Keplerian Orbit Elements. Available at www.rene-schwarz.com
#
#       https://downloads.rene-schwarz.com/dc/category/18
#       (Accessed on 2017-08-09).
#
#   [2] Vallado, D. A., McClain, W. D (2013). Fundamentals of Astrodynamics
#       and Applications. Microcosm Press.
#
#   [3] Kuga, H. K., Carrara, V., Rao, K. R (2005). Introdução à Mecânica
#       Orbital. 2ª ed. Instituto Nacional de Pesquisas Espaciais.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2017-08-09: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export kepler_to_rv
export rv_to_kepler

"""
### function kepler_to_rv(a::Number, e::Number, i::Number, Ω::Number, ω::Number, f::Number)

Convert the Keplerian elements (`a`, `e`, `i`, `Ω`, `ω`, and `f`) to a Cartesian
representation (position vector `r` and velocity vector `v`)

##### Args

* a: Semi-major axis [m].
* e: Excentricity.
* i: Inclination [rad].
* Ω: Right ascension of the ascending node [rad].
* ω: Argument of perigee [rad].
* f: True anomaly [rad].

##### Returns

* The position vector represented in the inertial reference frame [m].
* The velocity vector represented in the inertial reference frame [m].

##### References

This algorithm was adapted from [1] and [3, p. 37-38].

"""

################################################################################
#                                 TEST RESULTS
################################################################################
#
# This function was tested using the example in [1, pp. 119-120].
#
# Scenario
# ========
#
#    p    = 11067.790 km
#    e    = 0.83285
#    i    = 87.87°
#    RAAN = 227.89°
#    w    = 53.38°
#    f    = 92.335°
#    a    = p/(1-e^2)
#
#
#                    +----------+----------+-----------+
#                    | X [km/s] | Y [km/s] | Z [km/s]  |
# +------------------+----------+----------+-----------+
# | Algorithm in [1] | 4.902276 | 5.533124 | -1.975709 |
# +------------------+----------+----------+-----------+
# | This function    | 4.902279 | 5.533140 | -1.975710 |
# +------------------+----------+----------+-----------+
# | Difference       | 0.000003 | 0.000016 |  0.000001 |
# +------------------+----------+----------+-----------+
#
# References
# ==========
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th
# ed. Microcosm Press. Hawthorn, CA, USA.
#
# TODO: The position part of the algorithm must be tested.
#
################################################################################

function kepler_to_rv(a::Number,
                      e::Number,
                      i::Number,
                      Ω::Number,
                      ω::Number,
                      f::Number)
    # Check eccentricity.
    if !(0 <= e < 1)
        throw(ArgumentError("Eccentricity must be in the interval [0,1)."))
    end

    # Auxiliary variables.
    sin_f     = sin(f)
    cos_f     = cos(f)

    # Compute the geocentric distance.
    r = a*(1-e^2)/(1+e*cos_f)

    # Compute the position vector in the orbit plane, defined as:
    #   - The X axis points towards the perigee;
    #   - The Z axis is perpendicular to the orbital plane (right-hand);
    #   - The Y axis completes a right-hand coordinate system.
    r_o = SVector{3}(r*cos_f, r*sin_f, 0)

    # Compute the velocity vector in the orbit plane without perturbations.
    n = angvel(a, e, i, :J0)
    v_o = n*a/sqrt(1-e^2)*SVector{3}(-sin_f, e+cos_f, 0)

    # Compute the matrix that rotates the orbit reference frame into the
    # inertial reference frame.
    Dio = angle2dcm(-ω, -i, -Ω, :ZXZ)

    # Compute the position and velocity represented in the inertial frame.
    r_i = Dio*r_o
    v_i = Dio*v_o

    (r_i, v_i)
end

"""
### function rv_to_kepler(r::Vector, v::Vector)

Convert a Cartesian representation (position vector `r` and velocity vector `v`)
to the Keplerian elements.

##### Args

* r: Position vector in an inertial reference frame [m].
* v: Velocity vector in an inertial reference frame [m].

##### Returns

An instance of the structure `Orbit` with the Keplerian elements [SI units].

##### References

The algorithm was adapted from [1].

"""

################################################################################
#                                 TEST RESULTS
################################################################################
#
# This function was tested using the examples in [2, pp. 114 - 116].
#
# Example 2-5
# ===========
#
# Cartesian representation:
#
#     r = 6524.835    I + 6862.875    J + 6448.296    K km
#     v =    4.901327 I +    5.533756 J -    1.976341 K km

# Results of the conversion to Keplerian elements using this function compared
# to those in [2].
#
#     ╔═════════════════╦══════════════╦═════════════════════╗
#     ║    Parameter    ║   Vallado    ║ SatelliteToolbox.jl ║
#     ╠═════════════════╬══════════════╬═════════════════════╣
#     ║ Semi-major axis ║ 36127.343 km ║ 36127.349 km        ║
#     ║ Eccentricity    ║ 0.832853     ║ 0.832853            ║
#     ║ Inclination     ║ 87.870°      ║ 87.869°             ║
#     ║ RAAN            ║ 227.898°     ║ 227.898             ║
#     ║ Arg. of Perigee ║ 53.38        ║ 53.38°              ║
#     ║ True Anomaly    ║ 92.335°      ║ 92.335°             ║
#     ╚═════════════════╩══════════════╩═════════════════════╝
#
################################################################################

function rv_to_kepler(r_i::AbstractVector, v_i::AbstractVector)
    # Position and velocity vector norms.
    r2 = r_i'*r_i
    v2 = v_i'*v_i

    r  = sqrt(r2)
    v  = sqrt(v2)

    # Angular momentum vector.
    h_i = cross( r_i, v_i )
    h   = norm(h_i)

    # Vector that points to the right ascension of the ascending node (RAAN).
    n_i = cross( [0;0;1], h_i )
    n   = norm(n_i)

    # Eccentricity vector.
    e_i = ( (v2 - m0/r)*r_i - dot(r_i, v_i)*v_i )/m0

    # Orbit energy.
    ξ = v2/2 - m0/r

    # Eccentricity
    # ============

    ecc = norm(e_i)

    # Semi-major axis
    # ===============

    if abs(ecc) <= 1.0-1e-6
        a = -m0/(2*ξ)
    else
        error("Could not convert the provided Cartesian values to Kepler elements.\n" *
              "The computed eccentricity was not between 0 and 1");
    end

    # Inclination
    # ===========
    i = acos(h_i[3]/h)

    # Right Ascension of the Ascending Node.
    # ======================================
    Ω = acos(n_i[1]/n)
    (n_i[2] < 0) && (Ω = 2*pi - Ω)

    # Argument of Perigee
    # ===================
    ω = acos(n_i'*e_i/(n*ecc))
    (e_i[3] < 0) && (ω = 2*pi - ω)

    # True anomaly
    # ============
    v = acos(e_i'*r_i/(ecc*r))
    (r_i'*v_i < 0) && (v = 2*pi - v)

    # Return the Keplerian elements.
    # ==============================

    Orbit(0.0,a,ecc,i,Ω,ω,v)
end

"""
### function rv_to_kepler(x::Number, y::Number, z::Number, vx::Number, vy::Number, vz::Number)

Convert a Cartesian representation (position vector `[x;y;z]` and velocity
vector `[vx;vy;vz]`) to the Keplerian elements.

##### Args

* x: X component of the position vector in an inertial reference frame [m].
* y: Y component of the position vector in an inertial reference frame [m].
* z: Z component of the position vector in an inertial reference frame [m].
* vx: X component of the velocity vector in an inertial reference frame [m/s].
* vy: Y component of the velocity vector in an inertial reference frame [m/s].
* vz: Z component of the velocity vector in an inertial reference frame [m/s].

##### Returns

* The Keplerian elements in this order:
    - Semi-major axis [km].
    - Eccentricity.
    - Inclination [rad].
    - Right ascension of the ascending node [rad].
    - Argument of perigee [rad].
    - True anomaly [rad].

"""

function rv_to_kepler(x::Number,  y::Number,  z::Number,
                      vx::Number, vy::Number, vz::Number)
    # Create the position and velocity vectors.
    r_i = SVector{3}( x, y, z)
    v_i = SVector{3}(vx,vy,vz)

    # Compute the Keplerian orbit elements.
    rv_to_kepler(r_i,v_i)
end
