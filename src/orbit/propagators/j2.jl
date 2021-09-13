# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   J2 orbit propagator algorithm.
#
#   This algorithm propagates the orbit considering the perturbed two-body
#   equations as presented in [1, p. 690-692]. It uses the first-order
#   approximation of Kepler's problem, considering the effects of secular
#   gravitational and drag perturbations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] Wertz, J. R (1978). Spacecraft attitude determination and control.
#       Kluwer Academic Publishers, Dordrecht, The Netherlands.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export j2_gc_egm08, j2_gc_egm96, j2_gc_jgm02, j2_gc_jgm03
export j2_gc_egm08_f32, j2_gc_egm96_f32, j2_gc_jgm02_f32, j2_gc_jgm03_f32
export j2_init, j2!

################################################################################
#                                     TODO
################################################################################
#
# 1. Analyze the reference frame representation of the inputs for this
#    algorithm.
#
#   The SGP4 algorithm expects that the input parameters are represented in the
#   TEME (true equator, mean equinox) reference frame. This J2 orbit propagator
#   model requires that the input parameters are consistent with the
#   gravitational perturbation theory in which the `J2` coefficient was
#   computed. Looking at [1, p. 642], it appears that the perturbations are
#   considering a frame in which the Z-axis is aligned with the CIP (Celestial
#   Intermediate Pole, or the Earth rotation axis). Hence, the J2 parameter is
#   defined based on the PEF. Since no rotations or adaptations are programmed,
#   then the input parameters for this propagator should be represented in any
#   reference frame with a true Equator, because of the symmetry.
#
#   This needs to be further analyzed and confirmed.
#
################################################################################

################################################################################
#                                  Constants
################################################################################

# These constants were obtained from the GFC files. Remember that:
#
#   J_n = -C_n,0 * sqrt(2n+1)
#

# EGM-08 Gravitational constants.
const j2_gc_egm08 = J2_GravCte(
    6378137.0,
    sqrt(3.986004415e14/6378137.0^3),
    0.0010826261738522227
)

const j2_gc_egm08_f32 = J2_GravCte{Float32}(
    6378137.0,
    sqrt(3.986004415e14/6378137.0^3),
    0.0010826261738522227
)

# EGM-96 Gravitational constants.
const j2_gc_egm96 = J2_GravCte(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826266835531513
)

const j2_gc_egm96_f32 = J2_GravCte{Float32}(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826266835531513
)

# JGM-02 Gravitational constants.
const j2_gc_jgm02 = J2_GravCte(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826269256388149
)

const j2_gc_jgm02_f32 = J2_GravCte{Float32}(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826269256388149
)

# JGM-03 Gravitational constants.
const j2_gc_jgm03 = J2_GravCte(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826360229829945
)

const j2_gc_jgm03_f32 = J2_GravCte{Float32}(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826360229829945
)

################################################################################
#                                  Functions
################################################################################

"""
    j2_init(epoch::Tepoch, a_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number, dn_o2::Number = 0, ddn_o6::Number = 0; j2_gc::J2_GravCte{T} = j2_gc_egm08) where {Tepoch, T}

Initialize the data structure of J2 orbit propagator algorithm.

!!! note
    The type used in the propagation will be the same as used to define the
    gravitational constants in the structure `j2_gc`.

# Args

- `epoch::Number`: Epoch of the initial mean orbital elements [Julian Day].
- `a_0::Number`: Initial mean semi-major axis [m].
- `e_0::Number`: Initial mean eccentricity.
- `i_0::Number`: Initial mean inclination [rad].
- `Ω_0::Number`: Initial mean right ascension of the ascending node [rad].
- `ω_0::Number`: Initial mean argument of perigee [rad].
- `f_0::Number`: Initial mean true anomaly [rad].
- `dn_o2::Number`: First time derivative of the mean motion divided by two
    [rad/s^2].
- `ddn_o6::Number`: Second time derivative of the mean motion divided by six
    [rad/s^3].

# Keywords

- `j2_gc::J2_GravCte`: J2 orbit propagator gravitational constants (see
    [`J2_GravCte`](@ref)). (**Default** = `j2_gc_egm08`)

# Returns

The structure [`J2_Structure`](@ref) with the initialized parameters.
"""
function j2_init(
    epoch::Tepoch,
    a_0::Number,
    e_0::Number,
    i_0::Number,
    Ω_0::Number,
    ω_0::Number,
    f_0::Number,
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j2_gc::J2_GravCte{T} = j2_gc_egm08
) where {Tepoch, T}
    # Unpack the gravitational constants to improve code readability.
    @unpack_J2_GravCte j2_gc

    # Convert all inputs to the correct type.
    a_0    = T(a_0)
    e_0    = T(e_0)
    i_0    = T(i_0)
    Ω_0    = T(Ω_0)
    ω_0    = T(ω_0)
    f_0    = T(f_0)
    dn_o2  = T(dn_o2)
    ddn_o6 = T(ddn_o6)

    # Initial values.
    al_0 = a_0 / R0              # ............. Normalized semi-major axis [er]
    n_0  = μm / al_0^(T(3 / 2))  # ............. Unperturbed mean motion [rad/s]
    p_0  = al_0 * (1 - e_0^2)    # ...................... Semi-latus rectum [er]
    M_0  = f_to_M(e_0, f_0)      # .................. Initial mean anomaly [rad]

    # Auxiliary variables.
    dn   = 2dn_o2            # ..... Time-derivative of the mean motion [rad/s²]
    p_0² = p_0^2
    e_0² = e_0^2

    sin_i_0, cos_i_0 = sincos(i_0)
    sin_i_0² = sin_i_0^2

    # First-order time-derivative of the orbital elements.
    #
    # See [1, p 692].

    δa   = +T(2 / 3) * al_0 * dn / n_0
    δe   = +T(2 / 3) * (1 - e_0) * dn / n_0
    δΩ   = -T(3 / 2) * n_0 * J2 / p_0² * cos_i_0
    δω   = +T(3 / 4) * n_0 * J2 / p_0² * (4 - 5sin_i_0²)
    δM_0 = +T(3 / 4) * n_0 * J2 / p_0² * sqrt(1 - e_0²) * (2 - 3sin_i_0²)

    # Create the output structure with the data.
    return J2_Structure{Tepoch, T}(
        epoch  = epoch,
        al_0   = al_0,
        n_0    = n_0,
        e_0    = e_0,
        i_0    = i_0,
        Ω_0    = Ω_0,
        ω_0    = ω_0,
        f_0    = f_0,
        M_0    = M_0,
        dn_o2  = dn_o2,
        ddn_o6 = ddn_o6,
        j2_gc  = j2_gc,
        Δt     = 0,
        al_k   = al_0,
        e_k    = e_0,
        i_k    = i_0,
        Ω_k    = Ω_0,
        ω_k    = ω_0,
        f_k    = f_0,
        M_k    = M_0,
        δa     = δa,
        δe     = δe,
        δΩ     = δΩ,
        δω     = δω,
        δM_0   = δM_0,
    )
end

"""
    j2!(j2d::J2_Structure{Tepoch, T}, t::Number) where {Tepoch, T}

Propagate the orbit defined in `j2d` (see [`J2_Structure`](@ref)) until the time
`t` [s].

!!! note
    The internal values in `j2d` will be modified.

# Returns

- The position vector represented in the inertial frame at time `t` [m].
- The velocity vector represented in the inertial frame at time `t` [m/s]

# Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. Notice that the perturbation theory
requires an inertial frame with true equator.
"""
function j2!(j2d::J2_Structure{Tepoch, T}, t::Number) where {Tepoch, T}
    # Unpack the variables.
    @unpack al_0, n_0, e_0, i_0, Ω_0, ω_0, f_0, M_0, dn_o2, ddn_o6 = j2d
    @unpack δa, δe, δΩ, δω, δM_0, j2_gc = j2d
    @unpack R0, μm, J2 = j2_gc

    # Time elapsed since epoch.
    Δt = T(t)

    # Propagate the orbital elements.
    al_k = al_0 - δa * Δt
    e_k  = e_0 - δe * Δt
    i_k  = i_0
    Ω_k  = mod(Ω_0 + δΩ * Δt, T(2π))
    ω_k  = mod(ω_0 + δω * Δt, T(2π))

    # In [1, p. 692], it is mentioned that the mean anomaly must be updated
    # considering the equation:
    #                           .           ..
    #                          n_0          n_0
    #   M_k = M_0 + n_0 ⋅ Δt + ---- ⋅ Δt² + ---- ⋅ Δt³
    #                           2            6
    #
    # However, the time derivative of M_0 for perturbed orbits considering terms
    # up to J2 is not 0. This is clear in [1, p. 652].
    #
    # Hence, we will use the equation:
    #                                     .           ..
    #                .                    n_0          n_0
    #   M_k = M_0 + M_0 ⋅ Δt + n_0 ⋅ Δt + ---- ⋅ Δt² + ---- ⋅ Δt³
    #                                      2            6
    #

    M_k = mod(@evalpoly(Δt, M_0, (n_0 + δM_0), dn_o2, ddn_o6), T(2π))
    f_k = M_to_f(e_k, M_k)

    # Make sure that eccentricity is not lower than 0.
    e_k < 0 && (e_k = T(0))

    # Compute the position and velocity vectors given the orbital elements.
    r_i_k, v_i_k = kepler_to_rv(al_k * T(R0), e_k, i_k, Ω_k, ω_k, f_k)

    # Update the J2 orbit propagator structure.
    @pack! j2d = Δt, al_k, e_k, i_k, Ω_k, ω_k, M_k, f_k

    # Return the position and velocity vector represented in the inertial
    # reference frame.
    return r_i_k, v_i_k
end
