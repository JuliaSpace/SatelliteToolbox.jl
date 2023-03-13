# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#    API implementation for J2 osculating orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

get_epoch(orbp::OrbitPropagatorJ2osc) = orbp.j2oscd.j2d.epoch

function get_mean_elements(orbp::OrbitPropagatorJ2osc)
    j2d = orbp.j2oscd.j2d
    orb = KeplerianElements(
        j2d.epoch + j2d.Δt / 86400,
        j2d.al_k * j2d.j2c.R0,
        j2d.e_k,
        j2d.i_k,
        j2d.Ω_k,
        j2d.ω_k,
        j2d.f_k
    )

    return orb
end

"""
    init_orbit_propagator(Val(:J2osc), epoch::Number, a_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...)
    init_orbit_propagator(Val(:J2osc), orb_0::Orbit, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...)

Initialize the J2 osculating orbit propagator.

# Args

- `epoch::Number`: Initial orbit epoch [Julian Day].
- `a_0::Number`: Initial mean semi-major axis [m].
- `e_0::Number`: Initial mean eccentricity.
- `i_0::Number`: Initial mean inclination [rad].
- `Ω_0::Number`: Initial mean right ascension of the ascending node [rad].
- `ω_0::Number`: Initial mean argument of perigee [rad].
- `f_0::Number`: Initial mean true anomaly [rad].
- `dn_o2::Number`: (OPTIONAL) First time derivative of mean motion divided by 2
    \\[rad/s²] (**Default** = 0).
- `ddn_o6::Number`: (OPTIONAL) Second time derivative of mean motion divided by 6
    \\[rad/s³] (**Default** = 0).
- `orb_0::Orbit`: Object of type [`Orbit`](@ref) with the initial mean orbital
    elements [SI].

# Keywords

- `j2c::J2PropagatorConstants`: J2 orbit propagator constants.
    (**Default** = `j2c_egm08`).
"""
function init_orbit_propagator(
    ::Val{:J2osc},
    epoch::Number,
    a_0::Number,
    e_0::Number,
    i_0::Number,
    Ω_0::Number,
    ω_0::Number,
    f_0::Number,
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j2c::J2PropagatorConstants{T} = j2c_egm08
) where T
    # Create the new Two Body propagator structure.
    j2oscd = j2osc_init(
        epoch,
        a_0,
        e_0,
        i_0,
        Ω_0,
        ω_0,
        f_0,
        dn_o2,
        ddn_o6;
        j2c = j2c
    )

    # Create and return the orbit propagator structure.
    return OrbitPropagatorJ2osc(j2oscd)
end

function init_orbit_propagator(
    ::Val{:J2osc},
    orb_0::Orbit,
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j2c::J2PropagatorConstants = j2c_egm08
)
    # Convert the orbit representation to Keplerian elements.
    k_0 = convert(KeplerianElements, orb_0)

    return init_orbit_propagator(
        Val(:J2osc),
        k_0.t,
        k_0.a,
        k_0.e,
        k_0.i,
        k_0.Ω,
        k_0.ω,
        k_0.f,
        dn_o2,
        ddn_o6;
        j2c = j2c
    )
end

function propagate!(orbp::OrbitPropagatorJ2osc, t::Number)
    # Auxiliary variables.
    j2oscd = orbp.j2oscd

    # Propagate the orbit.
    r_i, v_i = j2osc!(j2oscd, t)

    # Return.
    return r_i, v_i
end

function step!(orbp::OrbitPropagatorJ2osc, Δt::Number)
    # Auxiliary variables.
    j2oscd = orbp.j2oscd

    # Propagate the orbit.
    r_i, v_i = j2osc!(j2oscd, j2oscd.Δt + Δt)

    # Return the information about the step.
    return r_i, v_i
end
