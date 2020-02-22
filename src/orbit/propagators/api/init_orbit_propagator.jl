#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    SatelliteToolbox orbit propagator API: init_orbit_propagator
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Hoots, F. R., Roehrich, R. L (1980). Models for Propagation of NORAD
#       Elements Set. Spacetrack Report No. 3.
#
#   [2] Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S (2006). Revisiting
#       Spacetrack Report #3: Rev1. AIAA.
#
#   [3] SGP4 Source code of STRF: https://github.com/cbassa/strf
#       The SGP4 C code available on STRF was converted by Paul. S. Crawford and
#       Andrew R. Brooks.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export init_orbit_propagator

"""
    init_orbit_propagator(T, epoch::Number, a_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number, ...)
    init_orbit_propagator(T, orb_0::Orbit, ...)

Initialize the orbit propagator `T` using the initial mean orbital elements. The
propagator type `T` can be:

* `Val{:J2}`: J2 orbit propagator;
* `Val{:J4}`: J4 orbit propagator; or
* `Val{:twobody}`: Two-body orbit propagator.

The mean orbital elements can be passed individually of using an instance of the
structure `Orbit`.

# Args

* `epoch`: Initial orbit epoch [Julian Day].
* `a_0`: Initial semi-major axis [m].
* `e_0`: Initial eccentricity.
* `i_0`: Initial inclination [rad].
* `Ω_0`: Initial right ascension of the ascending node [rad].
* `ω_0`: Initial argument of perigee [rad].
* `f_0`: Initial true anomaly [rad].
* `n_0`: Initial angular velocity [rad/s].
* `M_0`: Initial mean anomaly [rad].
* `orb_0`: Instance of the structure `Orbit` with the initial mean orbital
           elements [SI].

## Additional optional arguments for the J2 orbit propagator

The initialization function of the J2 orbit propagator can receive the
following optional parameters:

* `dn_o2`: (OPTIONAL) First time derivative of mean motion divided by 2
           \\[rad/s²] (**Default** = 0).
* `ddn_o6`: (OPTIONAL) Second time derivative of mean motion divided by 6
            \\[rad/s³] (**Default** = 0).
* `j2_gc`: (OPTIONAL) J2 orbit propagator gravitational constants
           (**Default** = `j2_gc_egm08`).

## Additional optional arguments for the J4 orbit propagator

The initialization function of the J4 orbit propagator can receive the
following optional parameters:

* `dn_o2`: (OPTIONAL) First time derivative of mean motion divided by 2
           \\[rad/s²] (**Default** = 0).
* `ddn_o6`: (OPTIONAL) Second time derivative of mean motion divided by 6
            \\[rad/s³] (**Default** = 0).
* `j4_gc`: (OPTIONAL) J4 orbit propagator gravitational constants
           (**Default** = `j4_gc_egm08`).

## Additional optional arguments for the two body orbit propagator

The initialization function of the two body orbit propagator can receive the
following optional parameter:

* `μ`: (OPTIONAL) Standard gravitational parameter of the central body
       \\[m^3/s^2] (**Default** = `m0`).

# Returns

A new instance of the orbit propagator structure that stores the information of
the orbit propagator.

# Remarks

If the orbit is defined in terms of the angular velocity (mean motion) instead
of the semi-major axis, then it is possible to use the function `angvel_to_a` to
convert.

"""
function init_orbit_propagator(::Type{Val{:J2}}, epoch::Number, a_0::Number,
                               e_0::Number, i_0::Number, Ω_0::Number,
                               ω_0::Number, f_0::Number, dn_o2::Number = 0,
                               ddn_o6::Number = 0,
                               j2_gc::J2_GravCte{T} = j2_gc_egm08) where T

    # Create the new Two Body propagator structure.
    j2d = j2_init(j2_gc, epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, dn_o2, ddn_o6)

    # Create the `Orbit` structure.
    orb_0 = Orbit(epoch, j2d.al_0*j2_gc.R0, e_0, i_0, Ω_0, ω_0, j2d.f_k)

    # Create and return the orbit propagator structure.
    OrbitPropagatorJ2(orb_0, j2d)
end

function init_orbit_propagator(::Type{Val{:J4}}, epoch::Number, a_0::Number,
                               e_0::Number, i_0::Number, Ω_0::Number,
                               ω_0::Number, f_0::Number, dn_o2::Number = 0,
                               ddn_o6::Number = 0,
                               j4_gc::J4_GravCte{T} = j4_gc_egm08) where T

    # Create the new J4 propagator structure.
    j4d = j4_init(j4_gc, epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, dn_o2, ddn_o6)

    # Create the `Orbit` structure.
    orb_0 = Orbit(epoch, j4d.al_0*j4_gc.R0, e_0, i_0, Ω_0, ω_0, j4d.f_k)

    # Create and return the orbit propagator structure.
    OrbitPropagatorJ4(orb_0, j4d)
end

function init_orbit_propagator(::Type{Val{:twobody}}, epoch::Number,
                               a_0::Number, e_0::Number, i_0::Number,
                               Ω_0::Number, ω_0::Number, f_0::Number,
                               μ::T = m0) where T

    # Create the new Two Body propagator structure.
    tbd = twobody_init(epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, μ)

    # Create the `Orbit` structure.
    orb_0 = Orbit(epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0)

    # Create and return the orbit propagator structure.
    OrbitPropagatorTwoBody(orb_0, tbd)
end

init_orbit_propagator(::Type{Val{:J2}}, orb_0::Orbit, dn_o2::Number = 0,
                      ddn_o6::Number = 0, j2_gc::J2_GravCte = j2_gc_egm08) =
    init_orbit_propagator(Val{:J2}, orb_0.t, orb_0.a, orb_0.e, orb_0.i, orb_0.Ω,
                          orb_0.ω, orb_0.f, dn_o2, ddn_o6, j2_gc)

init_orbit_propagator(::Type{Val{:J4}}, orb_0::Orbit, dn_o2::Number = 0,
                      ddn_o6::Number = 0, j4_gc::J4_GravCte = j4_gc_egm08) =
    init_orbit_propagator(Val{:J4}, orb_0.t, orb_0.a, orb_0.e, orb_0.i, orb_0.Ω,
                          orb_0.ω, orb_0.f, dn_o2, ddn_o6, j4_gc)

init_orbit_propagator(::Type{Val{:twobody}}, orb_0::Orbit, μ::Number = m0) =
    init_orbit_propagator(Val{:twobody}, orb_0.t, orb_0.a, orb_0.e, orb_0.i,
                          orb_0.Ω, orb_0.ω, orb_0.f, μ)

"""
    init_orbit_propagator(T, tle::TLE, ...)

Initialize the orbit propagator `T` using the TLE `tle`. The propagator type `T`
can be:

* `Val{:J2}`: J2 orbit propagator;
* `Val{:J4}`: J4 orbit propagator;
* `Val{:twobody}`: Two-body orbit propagator; or
* `Val{:sgp4}`: SGP4 orbit propagator.

## Additional optional arguments for the J2 orbit propagator

The initialization function of the J2 orbit propagator can receive the
following optional parameter:

* `j2_gc`: (OPTIONAL) J2 orbit propagator gravitational constants
           (**Default** = `j2_gc_egm08`).

## Additional optional arguments for the J4 orbit propagator

The initialization function of the J4 orbit propagator can receive the
following optional parameter:

* `j4_gc`: (OPTIONAL) J4 orbit propagator gravitational constants
           (**Default** = `j4_gc_egm08`).

## Additional optional arguments for the two body orbit propagator

The initialization function of the two body orbit propagator can receive the
following optional parameter:

* `μ`: (OPTIONAL) Standard gravitational parameter of the central body
       \\[m^3/s^2] (**Default** = `m0`).

## Additional optional arguments for the SGP4 orbit propagator

The initialization function of the SGP4 orbit propagator can receive the
following optional parameter:

* `sgp4_gc`: (OPTIONAL) Gravitational constants (**Default** = `sgp4_gc_wgs84`).

# Returns

A new instance of the orbit propagator structure that stores the information of
the orbit propagator.

# Remarks

The SGP4 implementation includes also the deep space perturbations, which was
originally called SDP4 algorithm. Modern approaches, such as [2] and [3],
identifies if the selected orbit must be propagated using the deep space
perturbations and automatically applied them. This is sometimes called SGDP4
algorithm.

"""
function init_orbit_propagator(::Type{Val{:J2}},
                               tle::TLE,
                               j2_gc::J2_GravCte = j2_gc_egm08)

    # Unpack the gravitational constants to improve code readability.
    @unpack_J2_GravCte j2_gc

    # Constants.
    revday2radsec = 2π/86400    # Revolutions per day to radians per second.
    d2r           =  π/180      # Degrees to radians.

    # Obtain the data from the TLE.
    n_0 = tle.n*revday2radsec
    e_0 = tle.e
    i_0 = tle.i*d2r
    Ω_0 = tle.Ω*d2r
    ω_0 = tle.ω*d2r
    M_0 = tle.M*d2r

    # Auxiliary variables.
    e_0²    = e_0^2
    cos_i_0 = cos(i_0)

    # Recover the original mean motion (nll_0) and semi-major axis (all_0) from
    # the input elements and the same algorithm of SGP4.
    aux   = 1/2*J2*(3*cos_i_0^2-1)/(1-e_0²)^(3/2)
    a_1   = (μm/n_0)^(2/3)
    δ_1   = 3/2*aux/a_1^2
    a_0   = a_1*(1 - 1/3*δ_1 - δ_1^2 - 134/81*δ_1^3)
    δ_0   = 3/2*aux/a_0^2
    nll_0 = n_0/(1 + δ_0)
    a_0   = (μm/nll_0)^(2/3)*R0
    f_0   = M_to_f(e_0, M_0)

    init_orbit_propagator(Val{:J2},
                          tle.epoch,
                          a_0,
                          e_0,
                          i_0,
                          Ω_0,
                          ω_0,
                          f_0,
                          tle.dn_o2*2π/(86400)^2,
                          tle.ddn_o6*2π/(86400)^3,
                          j2_gc_egm08)
end

function init_orbit_propagator(::Type{Val{:J4}},
                               tle::TLE,
                               j4_gc::J4_GravCte = j4_gc_egm08)

    # Unpack the gravitational constants to improve code readability.
    @unpack_J2_GravCte j4_gc

    revday2radsec = 2π/86400    # Revolutions per day to radians per second.
    # Constants.
    d2r           =  π/180      # Degrees to radians.

    # Obtain the data from the TLE.
    n_0 = tle.n*revday2radsec
    e_0 = tle.e
    i_0 = tle.i*d2r
    Ω_0 = tle.Ω*d2r
    ω_0 = tle.ω*d2r
    M_0 = tle.M*d2r

    # Auxiliary variables.
    e_0²    = e_0^2
    cos_i_0 = cos(i_0)

    # Recover the original mean motion (nll_0) and semi-major axis (all_0) from
    # the input elements and the same algorithm of SGP4.
    aux   = 1/2*J2*(3*cos_i_0^2-1)/(1-e_0²)^(3/2)
    a_1   = (μm/n_0)^(2/3)
    δ_1   = 3/2*aux/a_1^2
    a_0   = a_1*(1 - 1/3*δ_1 - δ_1^2 - 134/81*δ_1^3)
    δ_0   = 3/2*aux/a_0^2
    nll_0 = n_0/(1 + δ_0)
    a_0   = (μm/nll_0)^(2/3)*R0
    f_0   = M_to_f(e_0, M_0)

    init_orbit_propagator(Val{:J4}, tle.epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0,
                          tle.dn_o2*2π/(86400)^2, tle.ddn_o6*2π/(86400)^3,
                          j4_gc_egm08)
end

function init_orbit_propagator(::Type{Val{:twobody}},
                               tle::TLE,
                               μ::Number = m0)

    # Constants.
    revday2radsec = 2π/86400    # Revolutions per day to radians per second.
    d2r           =  π/180      # Degrees to radians.

    # Obtain the data from the TLE.
    n_0 = tle.n*revday2radsec
    e_0 = tle.e
    i_0 = tle.i*d2r
    Ω_0 = tle.Ω*d2r
    ω_0 = tle.ω*d2r
    M_0 = tle.M*d2r

    # Since we do not have in the Two Body algorithms the coefficient `J2`, then
    # the semi-major axis will be recovered considering a Keplerian orbit.

    a_0 = (μ/n_0^2)^(1/3)
    f_0 = M_to_f(e_0, M_0)

    init_orbit_propagator(Val{:twobody},
                          tle.epoch,
                          a_0,
                          e_0,
                          i_0,
                          Ω_0,
                          ω_0,
                          f_0,
                          μ)
end

function init_orbit_propagator(::Type{Val{:sgp4}}, tle::TLE,
                               sgp4_gc::SGP4_GravCte{T} = sgp4_gc_wgs84) where T

    # Constants.
    d2r           =  π/180      # Degrees to radians.
    revday2radmin = 2π/1440     # Revolutions per day to radians per minute.

    # Obtain the data from the TLE.
    epoch = tle.epoch
    n_0   = tle.n*revday2radmin
    e_0   = tle.e
    i_0   = tle.i*d2r
    Ω_0   = tle.Ω*d2r
    ω_0   = tle.ω*d2r
    M_0   = tle.M*d2r
    bstar = tle.bstar

    # Create the new SGP4 structure.
    sgp4d = sgp4_init(sgp4_gc, tle.epoch, n_0, e_0, i_0, Ω_0, ω_0, M_0,
                      tle.bstar)

    # Create the `Orbit` structure.
    orb_0 = Orbit(epoch, 1000sgp4d.a_k*sgp4_gc.R0, e_0, i_0, Ω_0, ω_0,
                  M_to_f(e_0, M_0))

    # Create and return the orbit propagator structure.
    OrbitPropagatorSGP4(orb_0, sgp4_gc, sgp4d)
end

