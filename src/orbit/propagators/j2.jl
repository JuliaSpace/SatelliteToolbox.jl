#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   J2 orbit propagator algorithm.
#
#       This algorithm propagates the orbit considering the perturbed two-body
#       equations as presented in [1, p. 690-692]. It uses the first-order
#       approximation of Kepler's problem, considering the effects of secular
#       gravitational and drag perturbations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] Wertz, J. R (1978). Spacecraft attitude determination and control.
#       Kluwer Academic Publishers, Dordrecht, The Netherlands.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export j2_gc_egm08
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
#   This need to be further analyzed and confirmed.
#
################################################################################

################################################################################
#                                  Overloads
################################################################################

# Copy for J2_Structure
Base.copy(m::J2_Structure) = J2_Structure([ getfield(m, k) for k = 1:length(fieldnames(m)) ]...)

# Deepcopy for J2_Structure.
Base.deepcopy(m::J2_Structure) = J2_Structure([ deepcopy(getfield(m, k)) for k = 1:length(fieldnames(m)) ]...)

################################################################################
#                                  Constants
################################################################################

# EGM-08 Gravitational constants.
#
# These constants were obtained from the GFC file of EGM-08. Remember that:
#
#   J_n = -C_n,0 * sqrt(2n+1)
#
j2_gc_egm08 = J2_GravCte(
        6378137.0,
        sqrt(3.986004415e14/6378137.0^3),
         0.0010826261738522227,
       )

################################################################################
#                                  Functions
################################################################################

"""
    function j2_init(j2_gc::J2_GravCte{T}, epoch::Number, a_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number, dn_o2::Number, ddn_o6::Number) where T

Initialize the data structure of J2 orbit propagator algorithm.

# Args

* `j2_gc`: J2 orbit propagator gravitational constants (see `J2_GravCte`).
* `epoch`: Epoch of the orbital elements [Julian Day].
* `a_0`: Initial semi-major axis [m].
* `e_0`: Initial eccentricity.
* `i_0`: Initial inclination [rad].
* `Ω_0`: Initial right ascension of the ascending node [rad].
* `ω_0`: Initial argument of perigee [rad].
* `f_0`: Initial true anomaly [rad].
* `dn_o2`: First time derivative of the mean motion divided by two [rad/s^2].
* `ddn_o6`: Second time derivative of the mean motion divided by six [rad/s^3].

# Returns

The structure `J2_Structure` with the initialized parameters.

# Remarks

The inputs are the mean orbital elements.

"""
function j2_init(j2_gc::J2_GravCte{T}, epoch::Number, a_0::Number, e_0::Number,
                 i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number,
                 dn_o2::Number, ddn_o6::Number) where T

    # Unpack the gravitational constants to improve code readability.
    @unpack_J2_GravCte j2_gc

    # Initial values.
    al_0 = a_0/R0            # Normalized semi-major axis [er].
    n_0  = μm/al_0^(3/2)     # Unperturbed mean motion [rad/s].
    p_0  = al_0*(1-e_0^2)    # Semi-latus rectum [er].
    M_0  = f_to_M(e_0, f_0)  # Initial mean anomaly [rad].

    # Auxiliary variables.
    dn   = 2dn_o2            # Time-derivative of the mean motion [rad/s²].
    p_0² = p_0^2
    e_0² = e_0^2

    sin_i_0, cos_i_0 = sincos(i_0)
    sin_i_0² = sin_i_0^2

    # First-order time-derivative of the orbital elements.
    #
    # See [1, p 692].

    δa   = +2/3*al_0*dn/n_0
    δe   = +2/3*(1-e_0)*dn/n_0
    δΩ   = -3/2*n_0*J2/p_0²*cos_i_0
    δω   = +3/4*n_0*J2/p_0²*(4 - 5sin_i_0²)
    δM_0 = +3/4*n_0*J2/p_0²*sqrt(1-e_0²)*(2 - 3sin_i_0²)

    # Create the output structure with the data.
    J2_Structure{T}(
        epoch, al_0, n_0, e_0, i_0, Ω_0, ω_0, f_0, M_0, 0, dn_o2, ddn_o6, al_0,
        e_0, i_0, Ω_0, ω_0, f_0, M_0, δa, δe, δΩ, δω, δM_0, j2_gc
    )
end

"""
    function j2!(j2d::J2_Structure{T}, t::Number) where T

Propagate the orbit defined in `j2d` (see `J2_Structure`) until the time `t`
[s]. Notice that the values in `j2d` will be modified.

# Returns

* The position vector represented in the inertial frame at time `t` [m].
* The velocity vector represented in the inertial frame at time `t` [m/s]

###### Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. If the orbit parameters are obtained
from a TLE, then the inertial frame will be TEME. Notice, however, that the
perturbation theory requires an inertial frame with true equator.

"""
function j2!(j2d::J2_Structure{T}, t::Number) where T
    # Unpack the variables.
    @unpack_J2_Structure j2d
    @unpack_J2_GravCte   j2_gc

    # Time elapsed since epoch.
    Δt = t

    # Propagate the orbital elements.
    al_k = al_0 - δa*Δt
    e_k  = e_0 - δe*Δt
    i_k  = i_0
    Ω_k  = mod(Ω_0 + δΩ*Δt, 2π)
    ω_k  = mod(ω_0 + δω*Δt, 2π)

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

    M_k = mod(@evalpoly(Δt, M_0, (n_0 + δM_0), dn_o2, ddn_o6), 2π)
    f_k = M_to_f(e_k, M_k)

    # Make sure that eccentricity is not lower than 0.
    (e_k < 0) && (e_k = T(0))

    # Compute the position and velocity vectors given the orbital elements.
    (r_i_k, v_i_k) = kepler_to_rv(al_k*R0, e_k, i_k, Ω_k, ω_k, f_k)

    # Update the J2 orbit propagator structure.
    j2d.Δt   = Δt
    j2d.al_k = al_k
    j2d.e_k  = e_k
    j2d.i_k  = i_k
    j2d.Ω_k  = Ω_k
    j2d.ω_k  = ω_k
    j2d.M_k  = M_k
    j2d.f_k  = f_k

    # Return the position and velocity vector represented in the inertial
    # reference frame.
    (r_i_k, v_i_k)
end
