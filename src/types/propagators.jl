# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Structures related to the orbit propagators.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export OrbitPropagator

"""
    OrbitPropagator{Tepoch, T}

Abstract type of the orbit propagator. Every propagator structure must be a
subtype of this type and must implement the following API functions:

    propagate!(orbp, t::Number)
    propagate!(orbp, t::AbstractVector)
    propagate_to_epoch!(orbp, JD::Number)
    propagate_to_epoch!(orbp, JD::AbstractVector)
    step!(orbp, Δt::Number)

!!! note
    `Tepoch` is the type of the epoch representation and `T` is the type of the
    internal propagator variables.
"""
abstract type OrbitPropagator{Tepoch, T} end

#                          Two Body Orbit Propagator
# ==============================================================================

export TwoBodyPropagator, OrbitPropagatorTwoBody

"""
    TwoBodyPropagator{Tepoch, T}

Low level Two Body orbit propagator structure.
"""
@with_kw mutable struct TwoBodyPropagator{Tepoch, T}
    # Initial mean orbital elements
    # ==========================================================================
    epoch::Tepoch # .................... Epoch of the initial mean elements [JD]
    a_0::T        # ............... Initial mean normalized semi-major axis [m]
    n_0::T        # ................................ Initial mean motion [rad/s]
    e_0::T        # ............................. Initial mean eccentricity [  ]
    i_0::T        # ............................. Initial mean inclination [rad]
    Ω_0::T        # .................................... Initial mean RAAN [rad]
    ω_0::T        # ..................... Initial mean argument of perigee [rad]
    f_0::T        # ............................ Initial mean true anomaly [rad]
    M_0::T        # ............................ Initial mean mean anomaly [rad]
    μ::T          # . Std. gravitational parameter of the central body [m^3/s^2]

    # Current mean orbital elements
    # ==========================================================================
    Δt::T   # .............. Timespan from the epoch of the initial elements [s]
    f_k::T  # .................................. Current mean true anomaly [rad]
    M_k::T  # .................................. Current mean mean anomaly [rad]
end

"""
    OrbitPropagatorTwoBody{Tepoch, T} <: OrbitPropagator{Tepoch, T}

Structure that holds the information related to the Two Body orbit propagator.

# Fields

- `tbd::TwoBodyPropagator`: Structure that stores the Two Body orbit propagator
    data (see [`TwoBodyPropagator`](@ref)).
"""
struct OrbitPropagatorTwoBody{Tepoch, T} <: OrbitPropagator{Tepoch, T}
    tbd::TwoBodyPropagator{Tepoch, T}
end

#                             J2 Orbit Propagator
# ==============================================================================

export J2PropagatorConstants, J2Propagator, OrbitPropagatorJ2

"""
    J2PropagatorConstants{T}

Constants for the J2 orbit propagator.

# Fields

- `R0::T`: Earth equatorial radius [m].
- `μm::T`: √(GM / R0^3) [er/s]^(3/2).
- `J2::T`: The second gravitational zonal harmonic of the Earth.
"""
@with_kw struct J2PropagatorConstants{T}
    R0::T
    μm::T
    J2::T
end

"""
    J2Propagator{Tepoch, T}

Low level J2 orbit propagator structure.
"""
@with_kw mutable struct J2Propagator{Tepoch, T}
    # Initial mean orbital elements
    # ==========================================================================

    epoch::Tepoch # .................... Epoch of the initial mean elements [JD]
    al_0::T       # ............... Initial mean normalized semi-major axis [er]
    n_0::T        # ................................ Initial mean motion [rad/s]
    e_0::T        # ............................. Initial mean eccentricity [  ]
    i_0::T        # ............................. Initial mean inclination [rad]
    Ω_0::T        # .................................... Initial mean RAAN [rad]
    ω_0::T        # ..................... Initial mean argument of perigee [rad]
    f_0::T        # ............................ Initial mean true anomaly [rad]
    M_0::T        # ............................ Initial mean mean anomaly [rad]
    dn_o2::T      # .............. First time derivative of mean motion [rad/s²]
    ddn_o6::T     # ............. Second time derivative of mean motion [rad/s³]

    # J2 orbit propagator gravitational constants.
    j2c::J2PropagatorConstants{T}

    # Current mean orbital elements
    # ==========================================================================

    Δt::T   # .............. Timespan from the epoch of the initial elements [s]
    al_k::T # ..................... Current mean normalized semi-major axis [er]
    e_k::T  # ................................... Current mean eccentricity [  ]
    i_k::T  # ................................... Current mean inclination [rad]
    Ω_k::T  # .......................................... Current mean RAAN [rad]
    ω_k::T  # ........................... Current mean argument of perigee [rad]
    f_k::T  # .................................. Current mean true anomaly [rad]
    M_k::T  # .................................. Current mean mean anomaly [rad]

    # Auxiliary variables
    # ==========================================================================

    δa::T   # ........................... Semi-major axis time derivative [er/s]
    δe::T   # ............................... Eccentricity time derivative [1/s]
    δΩ::T   # ..................................... RAAN time derivative [rad/s]
    δω::T   # ...................... Argument of perigee time derivative [rad/s]
    n̄::T    # ....................................... "Mean" mean motion [rad/s]
end

"""
    OrbitPropagatorJ2{Tepoch, T} <: OrbitPropagator{Tepoch, T}

Structure that holds the information related to the J2 orbit propagator.

# Fields

- `j2d`: Structure that stores the J2 orbit propagator data (see
    [`J2Propagator`](@ref)).
"""
struct OrbitPropagatorJ2{Tepoch, T} <: OrbitPropagator{Tepoch, T}
    j2d::J2Propagator{Tepoch, T}
end

#                        J2 osculating orbit propagator
# ==============================================================================

export J2oscPropagator, OrbitPropagatorJ2osc

"""
    J2oscPropagator{Tepoch, T}

Low level J2 osculating orbit propagator structure.
"""
@with_kw mutable struct J2oscPropagator{Tepoch, T}
    # J2 orbit propagator to propagate the mean elements.
    j2d::J2Propagator{Tepoch, T}

    # Propagation time from epoch.
    Δt::T

    # Current osculating Keplerian elements
    # ==========================================================================
    a_k::T # ................................... Osculating semi-major axis [er]
    e_k::T # ....................................... Osculating eccentricity [ ]
    i_k::T # ...................................... Osculating inclination [rad]
    Ω_k::T # ............................................. Osculating RAAN [rad]
    ω_k::T # .............................. Osculating argument of perigee [rad]
    f_k::T # ..................................... Osculating true anomaly [rad]
    M_k::T # ..................................... Osculating mean anomaly [rad]
end

"""
    OrbitPropagatorJ2osc{Tepoch, T} <: OrbitPropagator{Tepoch, T}

Structure that holds the information related to the J2 osculating orbit
propagator.

# Fields

- `j2oscd`: Structure that stores the J2 osculating orbit propagator data (see
    [`J2oscPropagator`](@ref)).
"""
struct OrbitPropagatorJ2osc{Tepoch, T} <: OrbitPropagator{Tepoch, T}
    j2oscd::J2oscPropagator{Tepoch, T}
end

#                             J4 orbit propagator
# ==============================================================================

export J4PropagatorConstants, J4Propagator, OrbitPropagatorJ4

"""
    J4PropagatorConstants{T}

Gravitational constants for J4 orbit propagator.

# Fields

- `R0::T`: Earth equatorial radius [m].
- `μm::T`: √(GM / R0^3) [er/s]^(3/2).
- `J2::T`: The second gravitational zonal harmonic of the Earth.
- `J4::T`: The fourth gravitational zonal harmonic of the Earth.
"""
@with_kw struct J4PropagatorConstants{T}
    R0::T
    μm::T
    J2::T
    J4::T
end

"""
    J4Propagator{Tepoch, T}

Low level J4 orbit propagator structure.
"""
@with_kw mutable struct J4Propagator{Tepoch, T}
    # Initial mean orbital elements
    # ==========================================================================

    epoch::Tepoch # .................... Epoch of the initial mean elements [JD]
    al_0::T       # ............... Initial mean normalized semi-major axis [er]
    n_0::T        # ................................ Initial mean motion [rad/s]
    e_0::T        # ............................. Initial mean eccentricity [  ]
    i_0::T        # ............................. Initial mean inclination [rad]
    Ω_0::T        # .................................... Initial mean RAAN [rad]
    ω_0::T        # ..................... Initial mean argument of perigee [rad]
    f_0::T        # ............................ Initial mean true anomaly [rad]
    M_0::T        # ............................ Initial mean mean anomaly [rad]
    dn_o2::T      # .............. First time derivative of mean motion [rad/s²]
    ddn_o6::T     # ............. Second time derivative of mean motion [rad/s³]

    # J4 orbit propagator gravitational constants.
    j4c::J4PropagatorConstants{T}

    # Current mean orbital elements
    # ==========================================================================

    Δt::T   # .............. Timespan from the epoch of the initial elements [s]
    al_k::T # ..................... Current mean normalized semi-major axis [er]
    e_k::T  # ................................... Current mean eccentricity [  ]
    i_k::T  # ................................... Current mean inclination [rad]
    Ω_k::T  # .......................................... Current mean RAAN [rad]
    ω_k::T  # ........................... Current mean argument of perigee [rad]
    f_k::T  # .................................. Current mean true anomaly [rad]
    M_k::T  # .................................. Current mean mean anomaly [rad]

    # Auxiliary variables
    # ==========================================================================

    δa::T   # ........................... Semi-major axis time derivative [er/s]
    δe::T   # ............................... Eccentricity time derivative [1/s]
    δΩ::T   # ..................................... RAAN time derivative [rad/s]
    δω::T   # ...................... Argument of perigee time derivative [rad/s]
    n̄::T    # ....................................... "Mean" mean motion [rad/s]
end

"""
    OrbitPropagatorJ4{Tepoch, T} <: OrbitPropagator{Tepoch, T}

Structure that holds the information related to the J4 orbit propagator.

# Fields

- `j4d`: Structure that stores the J4 orbit propagator data (see
    [`J4Propagator`](@ref)).
"""
struct OrbitPropagatorJ4{Tepoch, T} <: OrbitPropagator{Tepoch, T}
    j4d::J4Propagator{Tepoch, T}
end

#                                     SGP4
# ==============================================================================

export OrbitPropagatorSGP4

"""
    OrbitPropagatorSGP4{Tepoch, T} <: OrbitPropagator{Tepoch, T}

Structure that holds the information related to the SGP4 propagator.

# Fields

- `sgp4d`: Structure that stores the SGP4 data (see [`Sgp4Propagator`](@ref)).
"""
struct OrbitPropagatorSGP4{Tepoch, T} <: OrbitPropagator{Tepoch, T}
    sgp4d::Sgp4Propagator{Tepoch, T}
end
