# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#    API implementation for J2 osculating orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

get_epoch(orbp::OrbitPropagatorJ2osc) = orbp.j2d.epoch

"""
    init_orbit_propagator(Val(:J2osc), epoch::Number, a_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number, dn_o2::Number = 0, ddn_o6::Number = 0; j2_gc::J2_GravCte{T} = j2_gc_egm08) where T
    init_orbit_propagator(Val(:J2osc), orb_0::Orbit, dn_o2::Number = 0, ddn_o6::Number = 0; j2_gc::J2_GravCte = j2_gc_egm08)

Initialize the J2 osculating orbit propagator.

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
function init_orbit_propagator(::Val{:J2osc}, epoch::Number, a_0::Number,
                               e_0::Number, i_0::Number, Ω_0::Number,
                               ω_0::Number, f_0::Number, dn_o2::Number = 0,
                               ddn_o6::Number = 0;
                               j2_gc::J2_GravCte{T} = j2_gc_egm08) where T

    # Create the new Two Body propagator structure.
    j2oscd = j2osc_init(epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, dn_o2, ddn_o6;
                        j2_gc = j2_gc)

    # Create the `KeplerianElements` structure.
    j2d = j2oscd.j2d
    orb_0 = KeplerianElements(j2d.epoch,
                              j2d.al_0*j2_gc.R0,
                              j2d.e_0,
                              j2d.i_0,
                              j2d.Ω_0,
                              j2d.ω_0,
                              j2d.f_0)

    # Create and return the orbit propagator structure.
    return OrbitPropagatorJ2osc(orb_0, j2oscd)
end

function init_orbit_propagator(::Val{:J2osc},
                               orb_0::Orbit,
                               dn_o2::Number = 0,
                               ddn_o6::Number = 0;
                               j2_gc::J2_GravCte = j2_gc_egm08)
    # Convert the orbit representation to Keplerian elements.
    k_0 = convert(KeplerianElements, orb_0)

    return init_orbit_propagator(Val(:J2osc), k_0.t, k_0.a, k_0.e, k_0.i, k_0.Ω,
                                 k_0.ω, k_0.f, dn_o2, ddn_o6; j2_gc = j2_gc)
end

function propagate!(orbp::OrbitPropagatorJ2osc{T}, t::Number) where T
    # Auxiliary variables.
    j2oscd = orbp.j2oscd

    # Propagate the orbit.
    r_i, v_i = j2osc!(j2oscd, t)

    # Update the elements in the `orb` structure, which are the mean elements.
    j2d = j2oscd.j2d
    orbp.orb = KeplerianElements(j2d.epoch + t/86400,
                                 j2d.al_k*j2d.j2_gc.R0,
                                 j2d.e_k,
                                 j2d.i_k,
                                 j2d.Ω_k,
                                 j2d.ω_k,
                                 j2d.f_k)

    # Create the osculating element object.
    orb = KeplerianElements(j2oscd.j2d.epoch + t/86400,
                            j2oscd.a_k,
                            j2oscd.e_k,
                            j2oscd.i_k,
                            j2oscd.Ω_k,
                            j2oscd.ω_k,
                            j2oscd.f_k)

    # Return.
    return orb, r_i, v_i
end

function step!(orbp::OrbitPropagatorJ2osc, Δt::Number)
    # Auxiliary variables.
    j2oscd = orbp.j2oscd

    # Propagate the orbit.
    r_i, v_i = j2osc!(j2oscd, j2oscd.Δt + Δt)

    # Update the elements in the `orb` structure.
    j2d = j2oscd.j2d
    orbp.orb = KeplerianElements(orb.t + Δt/86400,
                                 j2d.al_k*j2d.j2_gc.R0,
                                 j2d.e_k,
                                 j2d.i_k,
                                 j2d.Ω_k,
                                 j2d.ω_k,
                                 j2d.f_k)

    # Create the osculating element object.
    orb = KeplerianElements(j2oscd.j2d.epoch + t/86400,
                            j2oscd.a_k,
                            j2oscd.e_k,
                            j2oscd.i_k,
                            j2oscd.Ω_k,
                            j2oscd.ω_k,
                            j2oscd.f_k)

    # Return the information about the step.
    return orb, r_i, v_i
end
