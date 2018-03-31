#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    SatToolbox orbit propagator API for J2 orbit propagator algorithm.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-03-31: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export OrbitPropagatorJ2
export init_orbit_propagator, step!, propagate!

################################################################################
#                             Types and Structures
################################################################################

mutable struct OrbitPropagatorJ2
    orb::Orbit

    # J2 orbit propagator related fields.
    j2d::J2_Structure
end

################################################################################
#                                  Functions
################################################################################

"""
### function init_orbit_propagator(::Type{Val{:J2}}, t_0::Number, n_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, M_0::Number, dn_o2::Number = 0, ddn_o6::Number = 0, j2_gc::J2_GravCte = j2_gc_wgs84)

Initialize the J2 orbit propagator using the initial orbit specified by the
elements `t_0, `n_0, `e_0`, `i_0`, `Ω_0`, `ω_0`, and `M_0`, and the
gravitational parameters `j2_gc` (see `J2_GravCte`).

##### Args

* t_0: Initial orbit epoch [s].
* n_0: Initial angular velocity [rad/s].
* e_0: Initial eccentricity.
* i_0: Initial inclination [rad].
* Ω_0: Initial right ascension of the ascending node [rad].
* ω_0: Initial argument of perigee [rad].
* M_0: Initial mean anomaly [rad].
* dn_o2: (OPTIONAL) First time derivative of mean motion divided by 2 [rad/s²]
         (**DEFAULT** = 0).
* ddn_o6: (OPTIONAL) Second time derivative of mean motion divided by 6 [rad/s³]
          (**DEFAULT** = 0).
* j2_gc: (OPTIONAL) J2 orbit propagator gravitational constants (**DEFAULT** =
         `j2_gc_wgs84`).

##### Returns

A new instance of the structure `OrbitPropagatorJ2` that stores the information
of the orbit propagator.

"""

function init_orbit_propagator(::Type{Val{:J2}},
                               t_0::Number,
                               n_0::Number,
                               e_0::Number,
                               i_0::Number,
                               Ω_0::Number,
                               ω_0::Number,
                               M_0::Number,
                               dn_o2::Number = 0,
                               ddn_o6::Number = 0,
                               j2_gc::J2_GravCte = j2_gc_wgs84)
    # Create the new Two Body propagator structure.
    j2d = j2_init(j2_gc, t_0, n_0, e_0, i_0, Ω_0, ω_0, M_0, dn_o2, ddn_o6)

    # Create the `Orbit` structure.
    orb_0 = Orbit(t_0, j2d.a_0*j2_gc.R0, e_0, i_0, Ω_0, ω_0, j2d.f_k)

    # Create and return the orbit propagator structure.
    OrbitPropagatorJ2(orb_0, j2d)
end

"""
### function init_orbit_propagator(::Type{Val{:J2}}, orb_0::Orbit, dn_o2::Number = 0, ddn_o6::Number = 0, j2_gc::J2_GravCte = j2_gc_wgs84)

Initialize the J2 orbit propagator using the initial orbit specified in `orb_0`,
and the gravitational parameters in the structure `j2_gc`.

##### Args

* orb_0: Initial orbital elements (see `Orbit`).
* dn_o2: (OPTIONAL) First time derivative of mean motion divided by 2 [rad/s²]
         (**DEFAULT** = 0).
* ddn_o6: (OPTIONAL) Second time derivative of mean motion divided by 6 [rad/s³]
          (**DEFAULT** = 0).
* j2_gc: (OPTIONAL) J2 orbit propagator gravitational constants (**DEFAULT** =
         `j2_gc_wgs84`).

##### Returns

A new instance of the structure `OrbitPropagatorJ2` that stores the information
of the orbit propagator.

"""

function init_orbit_propagator(::Type{Val{:J2}},
                               orb_0::Orbit,
                               dn_o2::Number = 0,
                               ddn_o6::Number = 0,
                               j2_gc::J2_GravCte = j2_gc_wgs84)
    init_orbit_propagator(Val{:J2},
                          orb_0.t,
                          angvel(orb_0, :J2),
                          orb_0.e,
                          orb_0.i,
                          orb_0.Ω,
                          orb_0.ω,
                          f_to_M(orb_0.e, orb_0.f),
                          dn_o2,
                          ddn_o6,
                          j2_gc_wgs84)
end

"""
### function init_orbit_propagator(::Type{Val{:J2}}, tle::TLE, j2_gc::J2_GravCte = j2_gc_wgs84)

Initialize the J2 orbit propagator using the initial orbit specified in the TLE
`tle`. The orbit epoch `t0` will be defined as the number of seconds since the
beginning of the year (see `TLE.epoch_day`).

##### Args

* tle: TLE that will be used to initialize the propagator.
* j2_gc: (OPTIONAL) J2 orbit propagator gravitational constants (**DEFAULT** =
         `j2_gc_wgs84`).

##### Returns

A new instance of the structure `OrbitPropagatorJ2` that stores the information
of the orbit propagator.

"""

function init_orbit_propagator(::Type{Val{:J2}},
                               tle::TLE,
                               j2_gc::J2_GravCte = j2_gc_wgs84)
    init_orbit_propagator(Val{:J2},
                          tle.epoch_day*24*60*60,
                          tle.n*2*pi/(24*60*60),
                          tle.e,
                          tle.i*pi/180,
                          tle.Ω*pi/180,
                          tle.ω*pi/180,
                          tle.M*pi/180,
                          tle.dn_o2*2*pi/(24*60*60)^2,
                          tle.ddn_o6*2*pi/(24*60*60)^3,
                          j2_gc_wgs84)
end

"""
### function step!(orbp::OrbitPropagatorJ2, Δt::Number)

Propagate the orbit in `orbp` by `Δt` s using the J2 orbit propagator algorithm.
The new parameters will be written in `orbp`.

##### Args

* orbp: Propagator structure (see `OrbitPropagatorJ2`).
* Δt: Step time [s].

##### Returns

* The Keplerian elements represented in the inertial frame after the step (see
  `Orbit`) [SI units].
* The position vector represented in the inertial frame after the step [m].
* The velocity vector represented in the inertial frame after the step [m].

##### Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. If the orbit parameters are obtained
from a TLE, then the inertial frame will be TEME. Notice, however, that the
perturbation theory requires an inertial frame with true equator.

"""

function step!(orbp::OrbitPropagatorJ2, Δt::Number)
    # Auxiliary variables.
    orb = orbp.orb
    j2d = orbp.j2d

    # Propagate the orbit.
    (r_i, v_i) = j2!(j2d, orb.t + Δt)

    # Update the elements in the `orb` structure.
    orb.t += Δt
    orb.a  = j2d.a_k*j2d.j2_gc.R0
    orb.e  = j2d.e_k
    orb.i  = j2d.i_k
    orb.Ω  = j2d.Ω_k
    orb.ω  = j2d.ω_k
    orb.f  = j2d.f_k

    # Return the information about the step.
    (copy(orbp.orb), r_i, v_i)
end

"""
### function propagate!(orbp::OrbitPropagatorJ2, t::Vector)

Propagate the orbit in `orbp` using the time instants defined in the vector `t`
using the J2 orbit propagator. The structure `orbp` will contain the elements at
the last propagation instant.

##### Args

* orbp: Propagator structure (see `OrbitPropagatorJ2`).
* t: Time instants from orbit epoch in which the orbit will be propagated
     [s].

##### Returns

* An array with the mean Keplerian elements represented in inertial frame in
  each time instant (see `Orbit`) [SI units].
* An array with the position vector represented in inertial frame in each time
  instant [m].
* An array with the velocity vector represented in inertial frame in each time
  instant [m].

##### Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. If the orbit parameters are obtained
from a TLE, then the inertial frame will be TEME. Notice, however, that the
perturbation theory requires an inertial frame with true equator.

"""

function propagate!(orbp::OrbitPropagatorJ2, t::Vector)
    # Auxiliary variables.
    orb = orbp.orb
    j2d = orbp.j2d

    # Output.
    result_orb = Array{Orbit}(0)
    result_r   = Array{Vector}(0)
    result_v   = Array{Vector}(0)

    for k in t
        # Propagate the orbit.
        (r_i_k, v_i_k) = j2!(j2d, j2d.t_0 + k)

        # Update the elements in the `orb` structure.
        orb.t = j2d.t_0 + k
        orb.a = j2d.a_k*j2d.j2_gc.R0
        orb.e = j2d.e_k
        orb.i = j2d.i_k
        orb.Ω = j2d.Ω_k
        orb.ω = j2d.ω_k
        orb.f = j2d.f_k

        push!(result_orb, copy(orb))
        push!(result_r,   r_i_k)
        push!(result_v,   v_i_k)
    end

    (result_orb, result_r, result_v)
end
