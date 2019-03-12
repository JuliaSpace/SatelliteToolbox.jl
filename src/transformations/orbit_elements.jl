#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export change_oe_frame
export kepler_to_rv
export rv_to_kepler
export kepler_to_sv, sv_to_kepler

"""
    function change_oe_frame(a::Number, e::Number, i::Number, Ω::Number, ω::Number, f::Number, conv_args...)
    function change_oe_frame(oe::Orbit, conv_args...)

Change the reference frame of orbit elements. The orbit elements can be
specified by `a`, `e`, `i`, `Ω`, `ω`, and `f`, or the structure `oe` (see
`Orbit`).

The conversion arguments `conv_args` are **the same** arguments that one should
pass to the function `rECItoECI` to convert between the desired frames. For more
information, see the documentation of the function `rECItoECI`.

# Args

* `a`: Semi-major axis [m].
* `e`: Excentricity.
* `i`: Inclination [rad].
* `Ω`: Right-ascension of the ascending node [rad].
* `ω`: Argument of perigee [rad].
* `f`: True anomaly [rad].
* `conv_args...`: Conversion arguments, which are the same arguments that one
                  would pass to the function `rECItoECI` to convert between the
                  desired frames.

* `oe`: An instance of the structure `Orbit` with the orbit elements that will
        be converted [SI units].

# Returns

An instance of the structure `Orbit` with the Keplerian elements [SI units]
converted to the new frame.

# Examples

```julia-repl
julia> eop = get_iers_eop(:IAU1980);

julia> teme_epoch = DatetoJD(2016,6,1,11,0,0);

julia> tod_epoch  = DatetoJD(2016,1,1,0,0,0);

julia> oe_teme    = Orbit(0,
                          7130.982e3,
                          0.001111,
                          98.405*pi/180,
                          227.336*pi/180,
                          90*pi/180,
                          320*pi/180)

                 Orbit
  =====================================
                  t =        0.0
    Semi-major axis =     7130.9820 km
       Eccentricity =        0.001111
        Inclination =       98.4050 ˚
               RAAN =      227.3360 ˚
    Arg. of Perigee =       90.0000 ˚
       True Anomaly =      320.0000 ˚

julia> oe_j2000 = change_oe_frame(oe_teme, TEME(), J2000(), teme_epoch, eop)

                 Orbit
  ======================================
                  t =        0.0
    Semi-major axis =     7130.9820 km
       Eccentricity =        0.001111
        Inclination =       98.3365 ˚
               RAAN =      227.1345 ˚
    Arg. of Perigee =       90.0604 ˚
       True Anomaly =      320.0000 ˚

julia> oe_tod   = change_oe_frame(oe_teme, TEME(), teme_epoch, TOD(), tod_epoch, eop)

                 Orbit
  ======================================
                  t =        0.0
    Semi-major axis =     7130.9820 km
       Eccentricity =        0.001111
        Inclination =       98.4037 ˚
               RAAN =      227.3306 ˚
    Arg. of Perigee =       90.0014 ˚
       True Anomaly =      320.0000 ˚
```
"""
function change_oe_frame(a::Number,
                         e::Number,
                         i::Number,
                         Ω::Number,
                         ω::Number,
                         f::Number,
                         conv_args...)

    # The approach is to transform the orbit elements to Cartesian
    # representation, convert the frame, and then convert back to orbit
    # elements.
    #
    # NOTE: In my benchmarks, the operation with DCMs are faster than with
    # quaternions after the DCM representation was changed to SMatrix.

    r_o, v_o    = kepler_to_rv(a, e, i, Ω, ω, f)
    D_ECIf_ECIo = rECItoECI(DCM, conv_args...)
    r_f         = D_ECIf_ECIo*r_o
    v_f         = D_ECIf_ECIo*v_o

    rv_to_kepler(r_f, v_f)
end

function change_oe_frame(oe::Orbit, conv_args...)
    change_oe_frame(oe.a, oe.e, oe.i, oe.Ω, oe.ω, oe.f, conv_args...)
end

"""
    function kepler_to_rv(a::Number, e::Number, i::Number, Ω::Number, ω::Number, f::Number)
    function kepler_to_rv(o::Orbit)

Convert the Keplerian elements (`a`, `e`, `i`, `Ω`, `ω`, and `f`) to a Cartesian
representation (position vector `r` and velocity vector `v`). The Keplerian
elements can also be passed inside an instance of the `Orbit` structure.

# Args

* `a`: Semi-major axis [m].
* `e`: Excentricity.
* `i`: Inclination [rad].
* `Ω`: Right ascension of the ascending node [rad].
* `ω`: Argument of perigee [rad].
* `f`: True anomaly [rad].

# Returns

* The position vector represented in the inertial reference frame [m].
* The velocity vector represented in the inertial reference frame [m].

# References

This algorithm was adapted from [1] and [3, p. 37-38].

"""
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
    sin_f, cos_f = sincos(f)

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
    Dio = angle_to_dcm(-ω, -i, -Ω, :ZXZ)

    # Compute the position and velocity represented in the inertial frame.
    r_i = Dio*r_o
    v_i = Dio*v_o

    (r_i, v_i)
end

kepler_to_rv(o::Orbit) = kepler_to_rv(o.a, o.e, o.i, o.Ω, o.ω, o.f)

"""
    function rv_to_kepler(r::Vector, v::Vector)

Convert a Cartesian representation (position vector `r` [m] and velocity vector
`v` [m/s²]) to the Keplerian elements.

# Returns

An instance of the structure `Orbit` with the Keplerian elements [SI units].

# Remarks

The special cases are treated as follows:

* **Circular and equatorial**: the right ascension of the ascending node and the
  argument of perigee are set to 0. Hence, the true anomaly is equal to the true
  longitude.
* **Elliptical and equatorial**: the right ascension of the ascending node is
  set to 0. Hence, the argument of perigee is equal to the longitude of
  periapsis.
* **Circular and inclined**: the argument of perigee is set to 0. Hence, the
  true anomaly is equal to the argument of latitude.

# References

The algorithm was adapted from [1].

"""
function rv_to_kepler(r_i::AbstractVector, v_i::AbstractVector)
    # Check inputs.
    length(r_i) != 3 && error("The vector r_i must have 3 dimensions.")
    length(v_i) != 3 && error("The vector v_i must have 3 dimensions.")

    @inbounds begin

        # Position and velocity vector norms and auxiliary dot products.
        r2 = dot(r_i,r_i)
        v2 = dot(v_i,v_i)

        r  = sqrt(r2)
        v  = sqrt(v2)

        rv = dot(r_i,v_i)

        # The type of `r` will be the type of the orbit elements.
        T   = typeof(r)
        Tm0 = T(m0)

        # Angular momentum vector.
        h_i = cross( r_i, v_i )
        h   = norm(h_i)

        # Vector that points to the right ascension of the ascending node (RAAN).
        n_i = SVector{3}(0,0,1) × h_i
        n   = norm(n_i)

        # Eccentricity vector.
        e_i = ( (v2 - Tm0/r)*r_i - rv*v_i )/Tm0

        # Orbit energy.
        ξ = v2/2 - Tm0/r

        # Eccentricity
        # ============

        ecc = norm(e_i)

        # Semi-major axis
        # ===============

        if abs(ecc) <= 1.0-1e-6
            a = -Tm0/(2ξ)
        else
            error("Could not convert the provided Cartesian values to Kepler elements.\n" *
                  "The computed eccentricity was not between 0 and 1");
        end

        # Inclination
        # ===========

        cos_i = h_i[3]/h
        cos_i = abs(cos_i) > 1 ? sign(cos_i) : cos_i
        i     = acos(cos_i)

        # Check the type of the orbit to account for special cases
        # ======================================================================

        # Equatorial
        # ----------------------------------------------------------------------

        if abs(n) <= 1e-6

            # Right Ascension of the Ascending Node.
            # ======================================

            Ω = T(0)

            # Equatorial and elliptical
            # ------------------------------------------------------------------

            if abs(ecc) > 1e-6

                # Argument of Perigee
                # ===================

                cos_ω = e_i[1]/ecc
                cos_ω = abs(cos_ω) > 1 ? sign(cos_ω) : cos_ω
                ω     = acos(cos_ω)

                (e_i[3] < 0) && (ω = T(2π) - ω)

                # True anomaly
                # ============

                cos_v = dot(e_i,r_i)/(ecc*r)
                cos_v = abs(cos_v) > 1 ? sign(cos_v) : cos_v
                v     = acos(cos_v)

                (rv < 0) && (v = T(2π) - v)

            # Equatorial and circular
            # ------------------------------------------------------------------

            else
                # Argument of Perigee
                # ===================

                ω = T(0)

                # True anomaly
                # ============

                cos_v = r_i[1]/r
                cos_v = abs(cos_v) > 1 ? sign(cos_v) : cos_v
                v     = acos(cos_v)

                (r_i[2] < 0) && (v = T(2π) - v)
            end

        # Inclined
        # ----------------------------------------------------------------------
        else

            # Right Ascension of the Ascending Node.
            # ======================================

            cos_Ω = n_i[1]/n
            cos_Ω = abs(cos_Ω) > 1 ? sign(cos_Ω) : cos_Ω
            Ω     = acos(cos_Ω)

            (n_i[2] < 0) && (Ω = T(2π) - Ω)

            # Circular and inclined
            # ------------------------------------------------------------------

            if abs(ecc) < 1e-6

                # Argument of Perigee
                # ===================

                ω = T(0)

                # True anomaly
                # ============

                cos_v = dot(n_i,r_i)/(n*r)
                cos_v = abs(cos_v) > 1 ? sign(cos_v) : cos_v
                v     = acos(cos_v)

                (r_i[3] < 0) && (v = T(2π) - v)
            else

                # Argument of Perigee
                # ===================

                cos_ω = dot(n_i,e_i)/(n*ecc)
                cos_ω = abs(cos_ω) > 1 ? sign(cos_ω) : cos_ω
                ω     = acos(cos_ω)

                (e_i[3] < 0) && (ω = T(2π) - ω)

                # True anomaly
                # ============

                cos_v = dot(e_i,r_i)/(ecc*r)
                cos_v = abs(cos_v) > 1 ? sign(cos_v) : cos_v
                v     = acos(cos_v)

                (rv < 0) && (v = T(2π) - v)
            end
        end
    end

    # Return the Keplerian elements.
    # ==============================

    Orbit(0,a,ecc,i,Ω,ω,v)
end

"""
    function rv_to_kepler(x::Number, y::Number, z::Number, vx::Number, vy::Number, vz::Number)

Convert a Cartesian representation (position vector `[x;y;z]` [m] and velocity
vector `[vx;vy;vz]` [m/s²]) to the Keplerian elements.

# Returns

An instance of the structure `Orbit` with the Keplerian elements [SI units].

"""
function rv_to_kepler(x::Number,  y::Number,  z::Number,
                      vx::Number, vy::Number, vz::Number)
    # Create the position and velocity vectors.
    r_i = SVector{3}( x, y, z)
    v_i = SVector{3}(vx,vy,vz)

    # Compute the Keplerian orbit elements.
    rv_to_kepler(r_i,v_i)
end


"""
    function kepler_to_sv(orb::Orbit)

Convert the Keplerian elements in the structure `orb` to a state vector.

"""
function kepler_to_sv(orb::Orbit)
    r_i, v_i = kepler_to_rv(orb)
    satsv(orb.t, r_i, v_i)
end

"""
    function sv_to_kepler(sv::SatelliteStateVector)

Convert the state vector `sv` to Keplerian elements represented by an instance
of the structure `Orbit`.

"""
sv_to_kepler(sv::SatelliteStateVector) = rv_to_kepler(sv.r, sv.v)
