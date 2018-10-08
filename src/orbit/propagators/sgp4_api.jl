#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    SatelliteToolbox orbit propagator API for SGP4 algorithm.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export init_orbit_propagator, step!, propagate!

################################################################################
#                                    Macros
################################################################################

"""
### macro update_orb!(orbp, t)

Macro that updates the parameters of `orbp.orb` using `t` and `orbp.sgp4d`.

# Args

* `orbp`: Orbit propagator structure.
* `t`: New orbit epoch [s].

"""
macro update_orb!(orbp, t)
    quote
        local orb   = $(esc(orbp)).orb
        local sgp4d = $(esc(orbp)).sgp4d
        local R0    = $(esc(orbp)).sgp4_gc.R0

        orb.t = sgp4d.epoch + $(esc(t))/86400
        orb.a = sgp4d.a_k*R0*1000
        orb.e = sgp4d.e_k
        orb.i = sgp4d.i_k
        orb.Ω = sgp4d.Ω_k
        orb.ω = sgp4d.ω_k
        orb.f = M_to_f(sgp4d.e_k, sgp4d.M_k)
    end
end

################################################################################
#                                  Functions
################################################################################

"""
    function init_orbit_propagator(::Type{Val{:sgp4}}, epoch::Number, n_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, M_0::Number, bstar::Number, sgp4_gc::SGP4_GravCte{T} = sgp4_gc_wgs84) where T

Initialize the SGP4 orbit propagator using the initial orbit specified by the
elements `epoch, `n_0, `e_0`, `i_0`, `Ω_0`, `ω_0`, and `M_0`, the B* parameter
`bstar`, and the gravitational constants in the structure `sgp4_gc`.

Notice that the orbit elements **must be** represented in TEME frame.

# Args

* `epoch`: Initial orbit epoch [Julian Day].
* `n_0`: Initial angular velocity [rad/s].
* `e_0`: Initial eccentricity.
* `i_0`: Initial inclination [rad]
* `Ω_0`: Initial right ascension of the ascending node [rad].
* `ω_0`: Initial argument of perigee [rad].
* `M_0`: Initial mean anomaly [rad].
* `bstar`: Initial B* parameter of the SGP4.
* `sgp4_gc`: (OPTIONAL) Gravitational constants (**Default** = `sgp4_gc_wgs84`).

# Returns

A new instance of the structure `OrbitPropagatorSGP4` that stores the
information of the orbit propagator.

# Remarks

This implementation includes also the deep space perturbations, which was
originally called SDP4 algorithm. Modern approaches, such as [2] and [3],
identifies if the selected orbit must be propagated using the deep space
perturbations and automatically applied them. This is sometimes called SGDP4
algorithm.

"""
function init_orbit_propagator(::Type{Val{:sgp4}},
                               epoch::Number,
                               n_0::Number,
                               e_0::Number,
                               i_0::Number,
                               Ω_0::Number,
                               ω_0::Number,
                               M_0::Number,
                               bstar::Number,
                               sgp4_gc::SGP4_GravCte{T} = sgp4_gc_wgs84) where T
    # Create the new SGP4 structure.
    sgp4d = sgp4_init(sgp4_gc,
                      epoch,
                      n_0*60,
                      e_0,
                      i_0,
                      Ω_0,
                      ω_0,
                      M_0,
                      bstar)

    # Create the `Orbit` structure.
    orb_0 = Orbit{T,T,T,T,T,T,T}(epoch,
                                 sgp4d.a_k*sgp4_gc_wgs84.R0,
                                 e_0,
                                 i_0,
                                 Ω_0,
                                 ω_0,
                                 M_to_f(e_0, M_0))

    # Create and return the orbit propagator structure.
    OrbitPropagatorSGP4(orb_0, sgp4_gc, sgp4d)
end

"""
    function init_orbit_propagator(::Type{Val{:sgp4}}, orb_0::Orbit, bstar::Number = 0.0, sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84)

Initialize the SGP4 orbit propagator using the initial orbit specified in
`orb_0`, the B* parameter `bstar`, and the gravitational constants in the
structure `sgp4_gc`.

Notice that the orbit elements specified in `orb_0` **must be** represented in
TEME frame.

# Args

* `orb_0`: Initial orbital elements (see `Orbit`).
* `bstar`: B* parameter of the SGP4.
* `sgp4_gc`: (OPTIONAL) Gravitational constants (**Default** = `sgp4_gc_wgs84`).

# Returns

A new instance of the structure `OrbitPropagatorSGP4` that stores the
information of the orbit propagator.

# Remarks

This implementation includes also the deep space perturbations, which was
originally called SDP4 algorithm. Modern approaches, such as [2] and [3],
identifies if the selected orbit must be propagated using the deep space
perturbations and automatically applied them. This is sometimes called SGDP4
algorithm.

"""
function init_orbit_propagator(::Type{Val{:sgp4}},
                               orb_0::Orbit,
                               bstar::Number = 0.0,
                               sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84)
    init_orbit_propagator(Val{:sgp4},
                          orb_0.t,
                          angvel(orb_0, :J2),
                          orb_0.e,
                          orb_0.i,
                          orb_0.Ω,
                          orb_0.ω,
                          f_to_M(orb_0.e, orb_0.f),
                          bstar,
                          sgp4_gc)
end

"""
    function init_orbit_propagator(::Type{Val{:sgp4}}, tle::TLE, sgp4_gc::SGP4_Structure = sgp4_gc_wgs84)

Initialize the SGP4 orbit propagator using the initial orbit specified in the
TLE `tle`. The orbit epoch `t0` will be defined as the number of seconds since
the beginning of the year (see `TLE.epoch_day`).

# Args

* `tle`: TLE that will be used to initialize the propagator.
* `sgp4_gc`: (OPTIONAL) Gravitational constants (**Default** = `sgp4_gc_wgs84`).

# Returns

A new instance of the structure `OrbitPropagatorSGP4` that stores the
information of the orbit propagator.

# Remarks

This implementation includes also the deep space perturbations, which was
originally called SDP4 algorithm. Modern approaches, such as [2] and [3],
identifies if the selected orbit must be propagated using the deep space
perturbations and automatically applied them. This is sometimes called SGDP4
algorithm.

"""
function init_orbit_propagator(::Type{Val{:sgp4}},
                               tle::TLE,
                               sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84)

    init_orbit_propagator(Val{:sgp4},
                          tle.epoch,
                          tle.n*2*pi/(24*60*60),
                          tle.e,
                          tle.i*pi/180,
                          tle.Ω*pi/180,
                          tle.ω*pi/180,
                          tle.M*pi/180,
                          tle.bstar,
                          sgp4_gc)
end

"""
    function step!(orbp::OrbitPropagatorSGP4{T}, Δt::Number) where T

Propagate the orbit in `orbp` by `Δt` s using the SGP4 propagator. The new
parameters will be written in `orbp`.

# Args

* `orbp`: Propagator structure (see `OrbitPropagatorSGP4`).
* `Δt`: Step time [s].

# Returns

* The mean Keplerian elements represented in TEME frame after the step (see
  `Orbit`) [SI units].
* The position vector represented in TEME frame after the step [m].
* The velocity vector represented in TEME frame after the step [m].

"""
function step!(orbp::OrbitPropagatorSGP4{T}, Δt::Number) where T
    # Auxiliary variables.
    orb   = orbp.orb
    sgp4d = orbp.sgp4d

    # Propagate the orbit.
    (r_teme, v_teme) = sgp4!(sgp4d, sgp4d.Δt + Δt/60)

    # Convert km to m.
    r_teme *= 1000
    v_teme *= 1000

    # Update the elements in the `orb` structure.
    @update_orb!(orbp, sgp4d.Δt*60)

    # Return the information about the step.
    (copy(orbp.orb), r_teme, v_teme)
end

"""
    function propagate!(orbp::OrbitPropagatorSGP4{T}, t::AbstractVector) where T

Propagate the orbit in `orbp` using the time instants defined in the vector `t`
using the SGP4 propagator. The structure `orbp` will contain the orbit elements
at the last propagation instant.

# Args

* `orbp`: Propagator structure (see `OrbitPropagatorSGP4`).
* `t`: Time instants from orbit epoch in which the orbit will be propagated [s].

# Returns

* An array with the mean Keplerian elements represented in TEME frame in each
  time instant (see `Orbit`) [SI units].
* An array with the position vector represented in TEME frame in each time
  instant [m].
* An array with the velocity vector represented in TEME frame in each time
  instant [m].

"""
function propagate!(orbp::OrbitPropagatorSGP4{T}, t::AbstractVector) where T
    # Auxiliary variables.
    orb   = orbp.orb
    sgp4d = orbp.sgp4d

    # Output.
    result_orb = Array{Orbit{T,T,T,T,T,T,T}}(undef, 0)
    result_r   = Array{Vector{T}}(undef, 0)
    result_v   = Array{Vector{T}}(undef, 0)

    for k in t
        # Propagate the orbit.
        (r_teme_k, v_teme_k) = sgp4!(sgp4d, k/60)

        # Convert km to m.
        r_teme_k *= 1000
        v_teme_k *= 1000

        # Update the elements in the `orb` structure.
        @update_orb!(orbp, k)

        push!(result_orb, copy(orb))
        push!(result_r,   r_teme_k)
        push!(result_v,   v_teme_k)
    end

    (result_orb, result_r, result_v)
end
