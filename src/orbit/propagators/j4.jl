# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   J4 orbit propagator algorithm.
#
#   This algorithm propagates the orbit considering the secular perturbations of
#   central body zonal harmonics as presented in [1, p. 647-654] and in [2].
#   Notice that only the terms J2, J2², and J4 are considered, i.e. the J6 is
#   assumed to be 0. The effect of the drag is also taken into account. This
#   can be used as a propagator of mean elements for mission analysis in which
#   the satellite orbit is maintained.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] Merson, R. H (1961). The motion of a satellite in an axi-symmetric
#       gravitational field. Geophysical journal international, Vol. 4(1),
#       p. 17-52.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
export j4_gc_egm08, j4_gc_egm96, j4_gc_jgm02, j4_gc_jgm03
export j4_gc_egm08_f32, j4_gc_egm96_f32, j4_gc_jgm02_f32, j4_gc_jgm03_f32
export j4_init, j4!

################################################################################
#                                  Constants
################################################################################

# These constants were obtained from the GFC files. Remember that:
#
#   J_n = -C_n,0 * sqrt(2n+1)
#

# EGM-08 Gravitational constants.
const j4_gc_egm08 = J4_GravCte(
    6378137.0,
    sqrt(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227,
    -1.6198975999169731e-6
)

const j4_gc_egm08_f32 = J4_GravCte{Float32}(
    6378137.0,
    sqrt(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227,
    -1.6198975999169731e-6
)

# EGM-96 Gravitational constants.
const j4_gc_egm96 = J4_GravCte(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826266835531513,
    -1.619621591367e-6
)

const j4_gc_egm96_f32 = J4_GravCte{Float32}(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826266835531513,
    -1.619621591367e-6
)

# JGM-02 Gravitational constants.
const j4_gc_jgm02 = J4_GravCte(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826269256388149,
    -1.62042999e-6
)

const j4_gc_jgm02_f32 = J4_GravCte{Float32}(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826269256388149,
    -1.62042999e-6
)

# JGM-03 Gravitational constants.
const j4_gc_jgm03 = J4_GravCte(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826360229829945,
    -1.619331205071e-6
)

const j4_gc_jgm03_f32 = J4_GravCte{Float32}(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826360229829945,
    -1.619331205071e-6
)

################################################################################
#                                  Functions
################################################################################

"""
    j4_init(epoch::Number, a_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number, dn_o2::Number = 0, ddn_o6::Number = 0; j4_gc::J4_GravCte{T} = j4_gc_egm08) where T

Initialize the data structure of J4 orbit propagator algorithm.

!!! note
    The type used in the propagation will be the same as used to define the
    gravitational constants in the structure `j4_gc`.

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

* `j4_gc::J4_GravCtr`: J4 orbit propagator gravitational constants (see
    [`J4_GravCte`](@ref)). (**Default** = `j4_gc_egm08`)

# Returns

The structure [`J4_Structure`](@ref) with the initialized parameters.
"""
function j4_init(
    epoch::Tepoch,
    a_0::Number,
    e_0::Number,
    i_0::Number,
    Ω_0::Number,
    ω_0::Number,
    f_0::Number,
    dn_o2::Number = 0,
    ddn_o6::Number = 0;
    j4_gc::J4_GravCte{T} = j4_gc_egm08
) where {Tepoch, T}
    # Unpack the gravitational constants to improve code readability.
    @unpack_J4_GravCte j4_gc

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
    al_0 = a_0 / R0               # ............ Normalized semi-major axis [er]
    n_0  = μm / al_0^(T(3 / 2))   # ............ Unperturbed mean motion [rad/s]
    p_0  = al_0 * (1 - T(e_0)^2)  # ..................... Semi-latus rectum [er]
    M_0  = f_to_M(e_0, f_0)       # ................. Initial mean anomaly [rad]

    # Auxiliary variables.
    dn               = 2 * dn_o2  # .... Time-deriv. of the mean motion [rad/s²]
    e_0²             = e_0^2
    e_0⁴             = e_0^4
    sin_i_0, cos_i_0 = sincos(i_0)
    sin_i_0          = sin(i_0)
    sin_i_0²         = sin_i_0^2
    sin_i_0⁴         = sin_i_0^4
    aux              = (1 - e_0²)
    saux             = sqrt(aux)
    p_0²             = p_0^2
    p_0⁴             = p_0^4
    J2²              = J2^2

    # First-order time-derivative of the orbital elements.
    #
    # See [1, p 692].

    δa   = -T( 2 /  3) * al_0 * dn / n_0
    δe   = -T( 2 /  3) * (1 - e_0) * dn / n_0
    δΩ   = -T( 3 /  2) * n_0 * J2 /p_0² * cos_i_0 +
            T( 3 / 32) * n_0 * J2²/ p_0⁴ * cos_i_0 * (12 -  4e_0² - (80 +  5e_0²) * sin_i_0²) +
            T(15 / 32) * n_0 * J4 / p_0⁴ * cos_i_0 * ( 8 + 12e_0² - (14 + 21e_0²) * sin_i_0²)

    δω   = T( 3 /   4) * n_0 * J2  / p_0² * (4 - 5sin_i_0²) +
           T( 9 / 384) * n_0 * J2² / p_0⁴ * (     56e_0² + (760 -  36e_0²) * sin_i_0² - (890 +  45e_0²) * sin_i_0⁴) -
           T(15 / 128) * n_0 * J4  / p_0⁴ * (64 + 72e_0² - (248 + 252e_0²) * sin_i_0² + (196 + 189e_0²) * sin_i_0⁴)

    δM_0 = T( 3 /   4) * n_0 * J2 / p_0² * saux * (2 - 3sin_i_0²) +
           T( 3 / 512) * n_0 * J2²/ p_0⁴ / saux * (
               (         320e_0² - 280e_0⁴) +
               ( 1600 - 1568e_0² + 328e_0⁴) * sin_i_0² +
               (-2096 + 1072e_0² +  79e_0⁴) * sin_i_0⁴
           ) - T(45 / 128) * n_0 * J4 / p_0⁴ * saux * e_0² * (-8 + 40sin_i_0 - 35sin_i_0²)

    # Create the output structure with the data.
    J4_Structure{Tepoch, T}(
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
        j4_gc  = j4_gc,
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
    j4!(j4d::J4_Structure{T}, t::Number) where T

Propagate the orbit defined in `j4d` (see [`J4_Structure`](@ref)) until the time
`t` [s].

!!! note
    The internal values in `j4d` will be modified.

# Returns

- The position vector represented in the inertial frame at time `t` [m].
- The velocity vector represented in the inertial frame at time `t` [m/s]

# Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. Notice that the perturbation theory
requires an inertial frame with true equator.
"""
function j4!(j4d::J4_Structure{Tepoch, T}, t::Number) where {Tepoch, T}
    # Unpack the variables.
    @unpack epoch, al_0, n_0, e_0, i_0, Ω_0, ω_0, f_0, M_0 = j4d
    @unpack dn_o2, ddn_o6, j4_gc, δa, δe, δΩ, δω, δM_0 = j4d
    @unpack R0, μm, J2, J4 = j4_gc

    # Time elapsed since epoch.
    Δt = T(t)

    # Propagate the orbital elements.
    al_k = al_0 + δa * Δt
    e_k  = e_0  + δe * Δt
    i_k  = i_0
    Ω_k  = mod(Ω_0 + δΩ * Δt, T(2π))
    ω_k  = mod(ω_0 + δω * Δt, T(2π))
    M_k  = mod(M_0 + (δM_0 + n_0) * Δt, T(2π))
    f_k  = M_to_f(e_k, M_k)

    # Make sure that eccentricity is not lower than 0.
    e_k < 0 && (e_k = T(0))

    # Compute the position and velocity vectors given the orbital elements.
    r_i_k, v_i_k = kepler_to_rv(al_k * R0, e_k, i_k, Ω_k, ω_k, f_k)

    # Update the J4 orbit propagator structure.
    @pack! j4d = Δt, al_k, e_k, i_k, Ω_k, ω_k, M_k, f_k

    # Return the position and velocity vector represented in the inertial
    # reference frame.
    return r_i_k, v_i_k
end
