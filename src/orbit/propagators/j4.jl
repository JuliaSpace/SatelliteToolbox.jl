# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   J4 orbit propagator algorithm.
#
#   This algorithm propagates the orbit considering the secular perturbations of
#   central body zonal harmonics as presented in [1, p. 647-654, 692-692], which
#   is Kozai's method but neglecting long-periodic and short-periodic
#   perturbations.
#
#   The terms J2, J2², and J4 are considered, i.e. the J6 is assumed to be 0.
#   The effect of the drag is also taken into account. This can be used as a
#   propagator of mean elements for mission analysis in which the satellite
#   orbit is maintained.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] Hoots, F. R., Roehrich, R. L (1980). Models for Propagation of NORAD
#       Elements Set. Spacetrack Report No. 3.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export j4c_egm08, j4c_egm96, j4c_jgm02, j4c_jgm03
export j4c_egm08_f32, j4c_egm96_f32, j4c_jgm02_f32, j4c_jgm03_f32
export j4_init, j4!

################################################################################
#                                  Constants
################################################################################

# These constants were obtained from the GFC files. Remember that:
#
#   J_n = -C_n,0 * sqrt(2n+1)
#

# EGM-08 gravitational constants.
const j4c_egm08 = J4PropagatorConstants(
    6378137.0,
    sqrt(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227,
    -1.6198975999169731e-6
)

const j4c_egm08_f32 = J4PropagatorConstants{Float32}(
    6378137.0,
    sqrt(3.986004415e14 / 6378137.0^3),
    0.0010826261738522227,
    -1.6198975999169731e-6
)

# EGM-96 gravitational constants.
const j4c_egm96 = J4PropagatorConstants(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826266835531513,
    -1.619621591367e-6
)

const j4c_egm96_f32 = J4PropagatorConstants{Float32}(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826266835531513,
    -1.619621591367e-6
)

# JGM-02 gravitational constants.
const j4c_jgm02 = J4PropagatorConstants(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826269256388149,
    -1.62042999e-6
)

const j4c_jgm02_f32 = J4PropagatorConstants{Float32}(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826269256388149,
    -1.62042999e-6
)

# JGM-03 gravitational constants.
const j4c_jgm03 = J4PropagatorConstants(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826360229829945,
    -1.619331205071e-6
)

const j4c_jgm03_f32 = J4PropagatorConstants{Float32}(
    6378136.3,
    sqrt(3.986004415e14/6378136.3^3),
    0.0010826360229829945,
    -1.619331205071e-6
)

################################################################################
#                                  Functions
################################################################################

"""
    j4_init(epoch::Number, a_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number, dn_o2::Number = 0, ddn_o6::Number = 0; kwargs...)

Initialize the data structure of J4 orbit propagator algorithm.

!!! note
    The type used in the propagation will be the same as used to define the
    constants in the structure `j4c`.

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

- `j4c::J4_GravCtr`: J4 orbit propagator constants (see
    [`J4PropagatorConstants`](@ref)). (**Default** = `j4c_egm08`)

# Returns

The structure [`J4Propagator`](@ref) with the initialized parameters.
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
    j4c::J4PropagatorConstants{T} = j4c_egm08
) where {Tepoch, T}
    # Unpack the gravitational constants to improve code readability.
    @unpack R0, μm, J2, J4 = j4c

    # Initial values.
    al_0 = T(a_0) / R0            # ............ Normalized semi-major axis [er]
    n_0  = μm / al_0^(T(3 / 2))   # ............ Unperturbed mean motion [rad/s]
    p_0  = al_0 * (1 - T(e_0)^2)  # ..................... Semi-latus rectum [er]
    M_0  = f_to_M(T(e_0), T(f_0)) # ................. Initial mean anomaly [rad]

    # Auxiliary variables.
    dn               = 2 * dn_o2  # .... Time-deriv. of the mean motion [rad/s²]
    e_0²             = T(e_0)^2
    sin_i_0, cos_i_0 = sincos(T(i_0))
    sin_i_0²         = sin_i_0^2
    sin_i_0⁴         = sin_i_0^4
    aux              = (1 - e_0²)
    saux             = sqrt(aux)
    p_0²             = p_0^2
    p_0⁴             = p_0^4
    J2²              = J2^2

    # We need to compute the "mean" mean motion that is used to calculate the
    # first-order time derivative of the orbital elements.
    #
    # NOTE: Description of J4 propagator in [1, p. 648-653].
    #
    # Using the equations in [1], we could not match the results from STK. After
    # analyzing the perturbation equations, it turns out that the
    # time-derivative depends on the mean motion instead of the unperturbed mean
    # motion. We can see this by looking at the algorithm in Kozai's method in
    # [1, p. 693].
    #
    # Notice that using the full expression here, with the J2² and J4 terms,
    # yields a solution with much higher error compared with STK result.
    āl = al_0 * (1 - T(3 / 4) * J2 / p_0² * saux * (2 - 3sin_i_0²))
    p̄  = āl   * aux
    n̄  = n_0  * (1 + T(3 / 4) * J2 / p̄^2 * saux * (2 - 3sin_i_0²))

    # First-order time-derivative of the orbital elements.
    δa = -T( 2 / 3  ) * al_0 * dn / n_0
    δe = -T( 2 / 3  ) * (1 - e_0) * dn / n_0

    # TODO: Check J4 perturbation term sign.
    #
    # We needed to flip the J4 perturbation term sign to obtain values that
    # match those of STK. However, this modification does not seem right if we
    # observe the RAAN secular perturbation term in SGP4 orbit propagator
    # [2, p. 16]. For more information, see:
    #
    #   https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/91
    #
    δΩ = -T( 3 / 2  ) * n̄ * J2  / p_0² * cos_i_0 +
          T( 3 / 32 ) * n̄ * J2² / p_0⁴ * cos_i_0 * (-36 -  4e_0² + 48saux + (40 - 5e_0² - 72saux) * sin_i_0²) +
          T(15 / 32 ) * n̄ * J4  / p_0⁴ * cos_i_0 * (  8 + 12e_0² - (14 + 21e_0²) * sin_i_0²)

    δω = +T( 3 / 4  ) * n̄ * J2  / p_0² * (4 - 5sin_i_0²) +
          T( 9 / 384) * n̄ * J2² / p_0⁴ * (192 + 56e_0² - 192saux + (-172 + 288saux) * sin_i_0² + e_0² * sin_i_0⁴) -
          T(15 / 128) * n̄ * J4  / p_0⁴ * (64 + 72e_0² - (248 + 252e_0²) * sin_i_0² + (196 + 189e_0²) * sin_i_0⁴)

    # Create the output structure with the data.
    J4Propagator{Tepoch, T}(
        epoch  = epoch,
        al_0   = al_0,
        n_0    = n_0,
        e_0    = T(e_0),
        i_0    = T(i_0),
        Ω_0    = T(Ω_0),
        ω_0    = T(ω_0),
        f_0    = T(f_0),
        M_0    = T(M_0),
        dn_o2  = T(dn_o2),
        ddn_o6 = T(ddn_o6),
        j4c    = j4c,
        Δt     = 0,
        al_k   = al_0,
        e_k    = T(e_0),
        i_k    = T(i_0),
        Ω_k    = T(Ω_0),
        ω_k    = T(ω_0),
        f_k    = T(f_0),
        M_k    = T(M_0),
        δa     = δa,
        δe     = δe,
        δΩ     = δΩ,
        δω     = δω,
        n̄      = n̄
    )
end

"""
    j4!(j4d::J4Propagator{T}, t::Number) where T

Propagate the orbit defined in `j4d` (see [`J4Propagator`](@ref)) until the time
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
function j4!(j4d::J4Propagator{Tepoch, T}, t::Number) where {Tepoch, T}
    # Unpack the variables.
    @unpack epoch, al_0, n_0, e_0, i_0, Ω_0, ω_0, f_0, M_0 = j4d
    @unpack dn_o2, ddn_o6, j4c, δa, δe, δΩ, δω, n̄ = j4d
    @unpack R0, μm, J2, J4 = j4c

    # Time elapsed since epoch.
    Δt = T(t)

    # Propagate the orbital elements.
    al_k = al_0 + δa * Δt
    e_k  = e_0  + δe * Δt
    i_k  = i_0
    Ω_k  = mod(Ω_0 + δΩ * Δt, T(2π))
    ω_k  = mod(ω_0 + δω * Δt, T(2π))
    M_k  = mod(M_0 + n̄  * Δt, T(2π))
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
