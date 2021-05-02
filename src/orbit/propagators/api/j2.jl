# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#    API implementation for J2 orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

get_epoch(orbp::OrbitPropagatorJ2) = orbp.j2d.epoch

"""
    init_orbit_propagator(Val(:J2), epoch::Number, a_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number, dn_o2::Number = 0, ddn_o6::Number = 0; j2_gc::J2_GravCte{T} = j2_gc_egm08) where T
    init_orbit_propagator(Val(:J2), orb_0::Orbit, dn_o2::Number = 0, ddn_o6::Number = 0; j2_gc::J2_GravCte = j2_gc_egm08)

Initialize the J2 orbit propagator.

# Args

* `epoch`: Initial orbit epoch [Julian Day].
* `a_0`: Initial mean semi-major axis [m].
* `e_0`: Initial mean eccentricity.
* `i_0`: Initial mean inclination [rad].
* `Ω_0`: Initial mean right ascension of the ascending node [rad].
* `ω_0`: Initial mean argument of perigee [rad].
* `f_0`: Initial mean true anomaly [rad].
* `dn_o2`: (OPTIONAL) First time derivative of mean motion divided by 2
           \\[rad/s²] (**Default** = 0).
* `ddn_o6`: (OPTIONAL) Second time derivative of mean motion divided by 6
            \\[rad/s³] (**Default** = 0).
* `orb_0`: Instance of the structure `KeplerianElements` with the initial mean
           orbital elements [SI].

## Keywords

* `j2_gc`: (OPTIONAL) J2 orbit propagator gravitational constants
           (**Default** = `j2_gc_egm08`).

"""
function init_orbit_propagator(::Val{:J2}, epoch::Number, a_0::Number,
                               e_0::Number, i_0::Number, Ω_0::Number,
                               ω_0::Number, f_0::Number, dn_o2::Number = 0,
                               ddn_o6::Number = 0;
                               j2_gc::J2_GravCte{T} = j2_gc_egm08) where T

    # Create the new Two Body propagator structure.
    j2d = j2_init(epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, dn_o2, ddn_o6;
                  j2_gc = j2_gc)

    # Create and return the orbit propagator structure.
    return OrbitPropagatorJ2(j2d)
end

function init_orbit_propagator(::Val{:J2},
                               orb_0::Orbit,
                               dn_o2::Number = 0,
                               ddn_o6::Number = 0;
                               j2_gc::J2_GravCte = j2_gc_egm08)
    # Convert the orbit representation to Keplerian elements.
    k_0 = convert(KeplerianElements, orb_0)

    return init_orbit_propagator(Val(:J2), k_0.t, k_0.a, k_0.e, k_0.i, k_0.Ω,
                                 k_0.ω, k_0.f, dn_o2, ddn_o6; j2_gc = j2_gc)
end

function propagate!(orbp::OrbitPropagatorJ2{T}, t::Number) where T
    # Auxiliary variables.
    j2d = orbp.j2d

    # Propagate the orbit.
    r_i, v_i = j2!(j2d, t)

    # Return.
    return r_i, v_i
end

function step!(orbp::OrbitPropagatorJ2, Δt::Number)
    # Auxiliary variables.
    j2d = orbp.j2d

    # Propagate the orbit.
    r_i, v_i = j2!(j2d, j2d.Δt + Δt)

    # Return the information about the step.
    return r_i, v_i
end
