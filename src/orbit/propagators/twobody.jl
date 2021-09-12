# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Two Body orbit propagator.
#
#   This algorithm considers a perfect Keplerian orbit. In other words, no
#   perturbation is considered during the propagation and the Earth is modeled
#   as a perfect sphere.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export TwoBody_Structure
export twobody_init, twobody!

################################################################################
#                                  Functions
################################################################################

"""
    twobody_init(epoch::Number, a_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number; μ::T = m0) where T

Initialize the data structure of two body orbit propagator algorithm.

!!! note
    The type used in the propagation will be the same as used to define the
    gravitational constant `μ`.

# Args

- `epoch::Number`: Epoch of the initial mean orbital elements [Julian Day].
- `a_0::Number`: Initial mean semi-major axis [m].
- `e_0::Number`: Initial mean eccentricity.
- `i_0::Number`: Initial mean inclination [rad].
- `Ω_0::Number`: Initial mean right ascension of the ascending node [rad].
- `ω_0::Number`: Initial mean argument of perigee [rad].
- `f_0::Number`: Initial mean true anomaly [rad].

# Keywords

- `μ::T`: Standard gravitational parameter of the central body [m^3/s^2].
    (**Default** = `m0`)

# Returns

The structure [`TwoBody_Structure`](@ref) with the initialized parameters.
"""
function twobody_init(
    epoch::Number,
    a_0::Number,
    e_0::Number,
    i_0::Number,
    Ω_0::Number,
    ω_0::Number,
    f_0::Number;
    μ::T = m0
) where T
    # The propagator is only defined for 0 <= e < 1.
    if !(0 <= e_0 < 1)
        throw(ArgumentError("The two body propagator only supports eccentricities in the interval [0,1)"))
    end

    # Compute the mean motion using the semi-major axis.
    n_0 = sqrt(μ / T(a_0)^3)

    # Compute the initial mean anomaly.
    M_0 = f_to_M(T(e_0), T(f_0))

    # Create and return the Two Body orbit propagator structure.
    TwoBody_Structure{T}(
        epoch = epoch,
        a_0   = a_0,
        n_0   = n_0,
        e_0   = e_0,
        i_0   = i_0,
        Ω_0   = Ω_0,
        ω_0   = ω_0,
        M_0   = M_0,
        f_0   = f_0,
        μ     = μ,
        Δt    = 0,
        f_k   = f_0,
        M_k   = M_0
    )
end

"""
    twobody!(tbd::TwoBody_Structure, t::Number)

Propagate the orbit defined in `tbd` (see [`TwoBody_Structure`](@ref)) until the
time `t` [s].

!!! note
    The internal values in `tbd` will be modified.

# Returns

- The position vector represented in the inertial frame at time `t` [m].
- The velocity vector represented in the inertial frame at time `t` [m/s]

# Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters.
"""
function twobody!(tbd::TwoBody_Structure{T}, t::Number) where T
    # Time elapsed since epoch.
    Δt     = T(t)
    tbd.Δt = Δt

    # Update the mean anomaly.
    tbd.M_k = tbd.M_0 + tbd.n_0 * Δt

    # Convert the mean anomaly to true anomaly.
    tbd.f_k = M_to_f(tbd.e_0, tbd.M_k)

    # Compute the position and velocity vectors given the orbital elements.
    r_i_k, v_i_k = kepler_to_rv(tbd.a_0, tbd.e_0, tbd.i_0, tbd.Ω_0, tbd.ω_0, tbd.f_k)

    # Return the position and velocity vector represented in the inertial
    # reference frame.
    return r_i_k, v_i_k
end
