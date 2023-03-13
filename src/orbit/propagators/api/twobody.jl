# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#    API implementation for two body orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

get_epoch(orbp::OrbitPropagatorTwoBody) = orbp.tbd.epoch

function get_mean_elements(orbp::OrbitPropagatorTwoBody)
    tbd = orbp.tbd
    orb = KeplerianElements(
        tbd.epoch + tbd.Δt / 86400,
        tbd.a_0,
        tbd.e_0,
        tbd.i_0,
        tbd.Ω_0,
        tbd.ω_0,
        tbd.f_k
    )

    return orb
end

"""
    init_orbit_propagator(::Val{:twobody}, epoch::Number, a_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number; kwargs...) where T

Initialize the two body orbit propagator.

# Args

- `epoch::Number`: Initial orbit epoch [Julian Day].
- `a_0::Number`: Initial mean semi-major axis [m].
- `e_0::Number`: Initial mean eccentricity.
- `i_0::Number`: Initial mean inclination [rad].
- `Ω_0::Number`: Initial mean right ascension of the ascending node [rad].
- `ω_0::Number`: Initial mean argument of perigee [rad].
- `f_0::Number`: Initial mean true anomaly [rad].
- `orb_0::Orbit`: Object of type [`Orbit`](@ref) with the initial mean orbital
    elements [SI].

## Keywords

- `μ::Number`: Standard gravitational parameter of the central body [m^3/s^2].
    (**Default** = `m0`).
"""
function init_orbit_propagator(
    ::Val{:twobody},
    epoch::Number,
    a_0::Number,
    e_0::Number,
    i_0::Number,
    Ω_0::Number,
    ω_0::Number,
    f_0::Number;
    μ::T = m0
) where T
    # Create the new Two Body propagator structure.
    tbd = twobody_init(epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0; μ = μ)

    # Create and return the orbit propagator structure.
    return OrbitPropagatorTwoBody(tbd)
end

function init_orbit_propagator(::Val{:twobody}, orb_0::Orbit; μ::Number = m0)
    # Convert the orbit representation to Keplerian elements.
    k_0 = convert(KeplerianElements, orb_0)

    return init_orbit_propagator(
        Val(:twobody),
        k_0.t,
        k_0.a,
        k_0.e,
        k_0.i,
        k_0.Ω,
        k_0.ω,
        k_0.f;
        μ = μ
    )
end

function propagate!(orbp::OrbitPropagatorTwoBody, t::Number)
    # Auxiliary variables.
    tbd = orbp.tbd

    # Propagate the orbit.
    r_i, v_i = twobody!(tbd, t)

    return r_i, v_i
end

function step!(orbp::OrbitPropagatorTwoBody, Δt::Number)
    # Auxiliary variables.
    tbd = orbp.tbd

    # Propagate the orbit.
    r_i, v_i = twobody!(tbd, tbd.Δt + Δt)

    # Return the information about the step.
    return r_i, v_i
end
