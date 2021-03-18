# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#    API implementation for SGP4 orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

get_epoch(orbp::OrbitPropagatorSGP4)    = orbp.sgp4d.epoch

"""
    init_orbit_propagator(Val(:sgp4), tle::TLE, sgp4_gc::SGP4_GravCte{T} = sgp4_gc_wgs84) where T

Initialize the SGP4 orbit propagator using the TLE `tle`.

## Keywords

* `sgp4_gc`: (OPTIONAL) Gravitational constants. (**Default** = `sgp4_gc_wgs84`)

"""
function init_orbit_propagator(::Val{:sgp4}, tle::TLE;
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

    # Create the `KeplerianElements` structure.
    orb_0 = KeplerianElements(epoch,
                              1000sgp4d.a_k*sgp4_gc.R0,
                              e_0,
                              i_0,
                              Ω_0,
                              ω_0,
                              M_to_f(e_0, M_0))

    # Create and return the orbit propagator structure.
    return OrbitPropagatorSGP4(orb_0, sgp4_gc, sgp4d)
end

function propagate!(orbp::OrbitPropagatorSGP4{T}, t::Number) where T
    # Auxiliary variables.
    sgp4d   = orbp.sgp4d
    sgp4_gc = orbp.sgp4_gc

    # Propagate the orbit.
    r_teme, v_teme = sgp4!(sgp4d, t/60)

    # Convert km to m.
    r_teme *= 1000
    v_teme *= 1000

    # Update the elements in the `orb` structure.
    orbp.orb = KeplerianElements(sgp4d.epoch + t/86400,
                                 sgp4d.a_k*sgp4_gc.R0*1000,
                                 sgp4d.e_k,
                                 sgp4d.i_k,
                                 sgp4d.Ω_k,
                                 sgp4d.ω_k,
                                 M_to_f(sgp4d.e_k, sgp4d.M_k))

    # Return.
    return orbp.orb, r_teme, v_teme
end

function step!(orbp::OrbitPropagatorSGP4{T}, Δt::Number) where T
    # Auxiliary variables.
    sgp4d   = orbp.sgp4d
    sgp4_gc = orbp.sgp4_gc

    # Propagate the orbit.
    r_teme, v_teme = sgp4!(sgp4d, sgp4d.Δt + Δt/60)

    # Convert km to m.
    r_teme *= 1000
    v_teme *= 1000

    # Update the elements in the `orb` structure.
    orbp.orb = KeplerianElements(sgp4d.epoch + sgp4d.Δt*60/86400,
                                 sgp4d.a_k*sgp4_gc.R0*1000,
                                 sgp4d.e_k,
                                 sgp4d.i_k,
                                 sgp4d.Ω_k,
                                 sgp4d.ω_k,
                                 M_to_f(sgp4d.e_k, sgp4d.M_k))

    # Return the information about the step.
    return orbp.orb, r_teme, v_teme
end
