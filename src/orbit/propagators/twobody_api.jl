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
#    SatelliteToolbox orbit propagator API for Two Body orbit propagator
#    algorithm.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-04-08: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Restrict types in the structures, which led to a huge performance gain.
#
# 2018-03-30: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export step!, propagate!

################################################################################
#                                  Functions
################################################################################

"""
    function init_orbit_propagator(::Type{Val{:twobody}}, t_0::Number, n_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, M_0::Number, μ::T = m0) where T

Initialize the Two Body orbit propagator using the initial orbit specified by
the elements `t_0, `n_0, `e_0`, `i_0`, `Ω_0`, `ω_0`, and `M_0`, and the standard
gravitational parameters of the central body `μ`.

##### Args

* `t_0`: Initial orbit epoch [s].
* `n_0`: Initial angular velocity [rad/s].
* `e_0`: Initial eccentricity.
* `i_0`: Initial inclination [rad].
* `Ω_0`: Initial right ascension of the ascending node [rad].
* `ω_0`: Initial argument of perigee [rad].
* `M_0`: Initial mean anomaly [rad].
* `μ`: (OPTIONAL) Standard gravitational parameter of the central body [m^3/s^2]
     (**Default** = `m0`).

##### Returns

A new instance of the structure `OrbitPropagatorTwoBody` that stores the
information of the orbit propagator.

"""
function init_orbit_propagator(::Type{Val{:twobody}},
                               t_0::Number,
                               n_0::Number,
                               e_0::Number,
                               i_0::Number,
                               Ω_0::Number,
                               ω_0::Number,
                               M_0::Number,
                               μ::T = m0) where T
    # Create the new Two Body propagator structure.
    tbd = twobody_init(t_0, n_0, e_0, i_0, Ω_0, ω_0, M_0, μ)

    # Create the `Orbit` structure.
    orb_0 = Orbit{T,T,T,T,T,T,T}(t_0, tbd.a, e_0, i_0, Ω_0, ω_0, tbd.f_k)

    # Create and return the orbit propagator strucutre.
    OrbitPropagatorTwoBody(orb_0, tbd)
end

"""
    function init_orbit_propagator(::Type{Val{:twobody}}, orb_0::Orbit, μ::Number = m0)

Initialize the Two Body orbit propagator using the initial orbit specified in
`orb_0`, and the standard gravitational parameters of the central body `μ`.

##### Args

* `orb_0`: Initial orbital elements (see `Orbit`).
* `μ`: (OPTIONAL) Standard gravitational parameter of the central body [m^3/s^2]
       (**Default** = `m0`).

##### Returns

A new instance of the structure `OrbitPropagatorTwoBody` that stores the
information of the orbit propagator.

"""
function init_orbit_propagator(::Type{Val{:twobody}},
                               orb_0::Orbit,
                               μ::Number = m0)
    init_orbit_propagator(Val{:twobody},
                          orb_0.t,
                          angvel(orb_0, :J0),
                          orb_0.e,
                          orb_0.i,
                          orb_0.Ω,
                          orb_0.ω,
                          f_to_M(orb_0.e, orb_0.f),
                          μ)
end

"""
    function init_orbit_propagator(::Type{Val{:twobody}}, tle::TLE, μ::Number = m0)

Initialize the Two Body orbit propagator using the initial orbit specified in
the TLE `tle`. The orbit epoch `t0` will be defined as the number of seconds
since the beginning of the year (see `TLE.epoch_day`).

##### Args

* `tle`: TLE that will be used to initialize the propagator.
* `μ`: (OPTIONAL) Standard gravitational parameter of the central body [m^3/s^2]
       (**Default** = `m0`).

##### Returns

A new instance of the structure `OrbitPropagatorTwoBody` that stores the
information of the orbit propagator.

"""
function init_orbit_propagator(::Type{Val{:twobody}},
                               tle::TLE,
                               μ::Number = m0)
    init_orbit_propagator(Val{:twobody},
                          tle.epoch_day*24*60*60,
                          tle.n*2*pi/(24*60*60),
                          tle.e,
                          tle.i*pi/180,
                          tle.Ω*pi/180,
                          tle.ω*pi/180,
                          tle.M*pi/180,
                          μ)
end

"""
    function step!(orbp::OrbitPropagatorTwoBody, Δt::Number)

Propagate the orbit in `orbp` by `Δt` s using the Two Body orbit propagator
algorithm. The new parameters will be written in `orbp`.

##### Args

* `orbp`: Propagator structure (see `OrbitPropagatorTwoBody`).
* `Δt`: Step time [s].

##### Returns

* The Keplerian elements represented in the inertial frame after the step (see
  `Orbit`) [SI units].
* The position vector represented in the inertial frame after the step [m].
* The velocity vector represented in the inertial frame after the step [m].

###### Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. If the orbit parameters are obtained
from a TLE, then the inertial frame will be TEME.

"""
function step!(orbp::OrbitPropagatorTwoBody, Δt::Number)
    # Auxiliary variables.
    orb = orbp.orb
    tbd = orbp.tbd

    # Propagate the orbit.
    (r_i, v_i) = twobody!(tbd, orb.t + Δt)

    # Update the elements in the `orb` structure.
    orb.t += Δt
    orb.f  = tbd.f_k

    # Return the information about the step.
    (copy(orbp.orb), r_i, v_i)
end

"""
    function propagate!(orbp::OrbitPropagatorTwoBody, t::Vector)

Propagate the orbit in `orbp` using the time instants defined in the vector `t`
using the Two Body orbit propagator. The structure `orbp` will contain the
elements at the last propagation instant.

##### Args

* `orbp`: Propagator structure (see `OrbitPropagatorTwoBody`).
* `t`: Time instants from orbit epoch in which the orbit will be propagated [s].

##### Returns

* An array with the mean Keplerian elements represented in inertial frame in
  each time instant (see `Orbit`) [SI units].
* An array with the position vector represented in inertial frame in each time
  instant [m].
* An array with the velocity vector represented in inertial frame in each time
  instant [m].

###### Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. If the orbit parameters are obtained
from a TLE, then the inertial frame will be TEME.

"""
function propagate!(orbp::OrbitPropagatorTwoBody, t::Vector)
    # Auxiliary variables.
    orb = orbp.orb
    tbd = orbp.tbd

    # Output.
    result_orb = Array{Orbit}(0)
    result_r   = Array{Vector}(0)
    result_v   = Array{Vector}(0)

    for k in t
        # Propagate the orbit.
        (r_i_k, v_i_k) = twobody!(tbd, tbd.t_0 + k)

        # Update the elements in the `orb` structure.
        orb.t = tbd.t_0 + k
        orb.f = tbd.f_k

        push!(result_orb, copy(orb))
        push!(result_r,   r_i_k)
        push!(result_v,   v_i_k)
    end

    (result_orb, result_r, result_v)
end
