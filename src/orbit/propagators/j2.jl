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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export j2_gc_wgs84, j2_gc_wgs72
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

# WGS-84 / EGM-08 Gravitational constants.
j2_gc_wgs84 = J2_GravCte(
        R0,
        sqrt(3.986005e14/6378137.0^3),
         0.00108262998905
       )

# WGS-72 Gravitational constants.
j2_gc_wgs72 = J2_GravCte(
        6378135.0,
        sqrt(3.986008e14/6378135.0^3),
         0.001082616
       )

################################################################################
#                                  Functions
################################################################################

"""
    function j2_init(j2_gc::J2_GravCte{T}, epoch::Number, n_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, M_0::Number, dn_o2::Number, ddn_o6::Number) where T

Initialize the data structure of J2 orbit propagator algorithm.

# Args

* `j2_gc`: J2 orbit propagator gravitational constants (see `J2_GravCte`).
* `epoch`: Epoch of the orbital elements [Julian Day].
* `n_0`: Mean motion at epoch [rad/s].
* `e_0`: Initial eccentricity.
* `i_0`: Initial inclination [rad].
* `Ω_0`: Initial right ascension of the ascending node [rad].
* `ω_0`: Initial argument of perigee [rad].
* `M_0`: Initial mean anomaly [rad].
* `dn_o2`: First time derivative of the mean motion divided by two [rad/s^2].
* `ddn_o6`: Second time derivative of the mean motion divided by six [rad/s^3].

# Returns

The structure `J2_Structure` with the initialized parameters.

"""
function j2_init(j2_gc::J2_GravCte{T},
                 epoch::Number,
                 n_0::Number,
                 e_0::Number,
                 i_0::Number,
                 Ω_0::Number,
                 ω_0::Number,
                 M_0::Number,
                 dn_o2::Number,
                 ddn_o6::Number) where T

    # Unpack the gravitational constants to improve code readability.
    @unpack_J2_GravCte j2_gc

    # Get the semi-major axis using J2 perturbation theory [er].
    #
    # This can only be done using a numerical algorithm to solve the following
    # equation for `a`:
    #
    #          μm        3    J2.μm    sqrt(1-e²).(3cos²(i)-1) + (5cos²(i)-1))
    #   n = --------- + ---.---------.-----------------------------------------
    #        a^(3/2)     4   a^(3/2)            a^2.(1-e^2)^2
    #
    # NOTE: This is necessary because we are specifying the angular velocity
    # instead of the semi-major axis.

    # Auxiliary variables to solve for the semi-major axis.
    f_ei = sqrt(1-e_0^2)*(3cos(i_0)^2-1) + (5cos(i_0)^2-1)
    K    = 3/4*J2*f_ei/(1-e_0^2)^2

    # Initial guess using a non-perturbed orbit.
    a_0 = (μm/n_0)^(2/3)

    # Newton-Raphson algorithm.
    #
    # Notice that we will allow, at most, 20 iterations.
    for k = 1:20
        # Auxiliary variables.
        a_0p3o2  = a_0^(3/2)
        a_0p5o2  = a_0p3o2*a_0 # -> a_0^(5/2)
        a_0p7o2  = a_0p5o2*a_0 # -> a_0^(7/2)
        a_0p11o2 = a_0p7o2*a_0 # -> a_0^(11/2)

        # Compute the residue.
        res = n_0 - μm/a_0p3o2 - K*μm/a_0p7o2

        # Compute the Jacobian of the function.
        df = +3/2*μm/a_0p5o2 + 7/2*K*μm/a_0p11o2

        # Compute the new estimate.
        a_0 = a_0 - res/df

        (abs(res) < 1e-10) && break
    end

    # Auxiliary variables.
    dn  = 2dn_o2
    p_0 = a_0*(1-e_0^2)
    f_0 = M_to_f(e_0, M_0)

    # Constants.
    #
    # See [1, p 692].

    C1 = 2/3*a_0*dn/n_0
    C2 = 2/3*(1-e_0)*dn/n_0
    C3 = 3/2*n_0*J2/p_0^2
    C4 = 3/4*n_0*J2/p_0^2

    # Create the output structure with the data.
    J2_Structure{T}(

        epoch, a_0, n_0, e_0, i_0, Ω_0, ω_0, M_0, 0, dn_o2, ddn_o6, a_0, e_0,
        i_0, Ω_0, ω_0, M_0, n_0, f_0, C1, C2, C3, C4, j2_gc

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

    # Auxiliary variables.
    sin_i_k, cos_i_k = sincos(i_k)

    # Propagate the orbital elements.
    a_k = a_0 - C1*Δt
    e_k = e_0 - C2*Δt
    i_k = i_0
    Ω_k = mod(Ω_0 - C3*cos_i_k*Δt,                    2π)
    ω_k = mod(ω_0 + C4*(4 - 5sin_i_k^2)*Δt,           2π)
    M_k = mod(@evalpoly(Δt, M_0, n_0, dn_o2, ddn_o6), 2π)
    f_k = M_to_f(e_k, M_k)

    # Make sure that eccentricity is not lower than 0.
    (e_k < 0) && (e_k = T(0))

    # Compute the position and velocity vectors given the orbital elements.
    (r_i_k, v_i_k) = kepler_to_rv(a_k*R0, e_k, i_k, Ω_k, ω_k, f_k)

    # Update the J2 orbit propagator structure.
    j2d.Δt  = Δt
    j2d.a_k = a_k
    j2d.e_k = e_k
    j2d.i_k = i_k
    j2d.Ω_k = Ω_k
    j2d.ω_k = ω_k
    j2d.M_k = M_k
    j2d.f_k = f_k

    # Return the position and velocity vector represented in the inertial
    # reference frame.
    (r_i_k, v_i_k)
end
