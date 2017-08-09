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
#    Conversions related to the orbit elements.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#     [1] Schwarz, R (2014). Memorandum No. 2: Cartesian State Vectors to
#         Keplerian Orbit Elements. Available at www.rene-schwarz.com
#
#         https://downloads.rene-schwarz.com/dc/category/18
#         (Accessed on 2017-08-09).
#
#    [2] Vallado, D. A., McClain, W. D (2013). Fundamentals of Astrodynamics
#        and Applications. Microcosm Press.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2017-08-09: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#     Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export rv_to_kepler

"""
### function rv_to_kepler(r::Vector{Float64}, v::Vector{Float64})

Convert a Cartesian representation (position and velocity) to the Keplerian
elements.

##### Args

* r_i: Position vector in an inertial reference frame [m].
* v_i: Velocity vector in an inertial reference frame [m].

##### Returns

* The Keplerian elements in this order:
    - Semi-major axis [km].
    - Eccentricity.
    - Inclination [rad].
    - Right ascension of the ascending node [rad].
    - Argument of perigee [rad].
    - True anomaly [rad].

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
#     ╔═════════════════╦══════════════╦═══════════════╗
#     ║    Parameter    ║   Vallado    ║ SatToolbox.jl ║
#     ╠═════════════════╬══════════════╬═══════════════╣
#     ║ Semi-major axis ║ 36127.343 km ║ 36127.349 km  ║
#     ║ Eccentricity    ║ 0.832853     ║ 0.832853      ║
#     ║ Inclination     ║ 87.870°      ║ 87.869°       ║
#     ║ RAAN            ║ 227.898°     ║ 227.898       ║
#     ║ Arg. of Perigee ║ 53.38        ║ 53.38°        ║
#     ║ True Anomaly    ║ 92.335°      ║ 92.335°       ║
#     ╚═════════════════╩══════════════╩═══════════════╝
#
################################################################################

function rv_to_kepler(r_i::Vector{Float64}, v_i::Vector{Float64})
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

    (a,ecc,i,Ω,ω,v)
end

"""
### function rv_to_kepler(x::Real, y::Real, z::Real, vx::Real, vy::Real, vz::Real)

Convert a Cartesian representation (position and velocity) to the Keplerian
elements.

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

function rv_to_kepler(x::Real, y::Real, z::Real, vx::Real, vy::Real, vz::Real)
    # Create the position and velocity vectors.
    r_i = [ x; y; z]
    v_i = [vx;vy;vz]

    # Compute the Keplerian orbit elements.
    rv_to_kepler(r_i,v_i)
end
