#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-04-08: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Restrict types in the structures, which led to a huge performance gain.
#
# 2017-08-08: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

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
#                             Types and Structures
################################################################################

# Gravitational constants for J2 orbit propagator.
immutable J2_GravCte{T}
    R0::T   # Earth equatorial radius [m].
    μm::T   # sqrt(GM) [er/s]^(3/2).
    J2::T   # The second gravitational zonal harmonic of the Earth.
end

# Serialization of the arguments in J2_GravCtr.
function getindex(j2_gc::J2_GravCte, ::Colon)
    j2_gc.R0, j2_gc.μm, j2_gc.J2
end

# J2 orbit propagator structure.
type J2_Structure{T}
    # Orbit parameters.
    t_0::T
    a_0::T
    n_0::T
    e_0::T
    i_0::T
    Ω_0::T
    ω_0::T
    M_0::T
    # Drag parameters.
    dn_o2::T   # First time derivative of mean motion [rad/s²].
    ddn_o6::T  # Second time derivative of mean motion [rad/s³].
    # Current parameters.
    a_k::T
    e_k::T
    i_k::T
    Ω_k::T
    ω_k::T
    M_k::T
    n_k::T
    f_k::T
    # Useful constants to decrease the computational burden.
    C1::T
    C2::T
    C3::T
    C4::T
    # J2 orbit propagator gravitational constants.
    j2_gc::J2_GravCte{T}
end

# Copy for J2_Structure
function Base.copy(m::J2_Structure)
    J2_Structure([ getfield(m, k) for k = 1:length(fieldnames(m)) ]...)
end

# Deepcopy for J2_Structure.
function Base.deepcopy(m::J2_Structure)
    J2_Structure([ deepcopy(getfield(m, k)) for k = 1:length(fieldnames(m)) ]...)
end

# Serialization of the arguments in J2_Structure.
function getindex(j2d::J2_Structure, ::Colon)

    j2d.t_0, j2d.a_0, j2d.n_0, j2d.e_0, j2d.i_0, j2d.Ω_0, j2d.ω_0, j2d.M_0,
    j2d.dn_o2, j2d.ddn_o6, j2d.a_k, j2d.e_k, j2d.i_k, j2d.Ω_k, j2d.ω_k, j2d.M_k,
    j2d.n_k, j2d.f_k, j2d.C1, j2d.C2, j2d.C3, j2d.C4, j2d.j2_gc

end

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
### function j2_init(j2_gc::J2_GravCte{T}, t_0::Number, n_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, M_0::Number, dn_o2::Number, ddn_o6::Number) where T

Initialize the data structure of J2 orbit propagator algorithm.

##### Args

* j2_gc: J2 orbit propagator gravitational constants (see `J2_GravCte`).
* t_0: Epoch of the orbital elements [s].
* n_0: Mean motion at epoch [rad/s].
* e_0: Initial eccentricity.
* i_0: Initial inclination [rad].
* Ω_0: Initial right ascension of the ascending node [rad].
* ω_0: Initial argument of perigee [rad].
* M_0: Initial mean anomaly [rad].
* dn_o2: First time derivative of the mean motion divided by two [rad/s^2].
* ddn_o6: Second time derivative of the mean motion divided by six [rad/s^3].

##### Returns

The structure `J2_Structure` with the initialized parameters.

"""
function j2_init(j2_gc::J2_GravCte{T},
                 t_0::Number,
                 n_0::Number,
                 e_0::Number,
                 i_0::Number,
                 Ω_0::Number,
                 ω_0::Number,
                 M_0::Number,
                 dn_o2::Number,
                 ddn_o6::Number) where T
    # Unpack the gravitational constants to improve code readability.
    R0, μm, J2 = j2_gc[:]

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

    # Declare the function in which the zero must be find.
    f(a)  = n_0 - μm/a^(3/2) - K*μm/a^(7/2)
    df(a) = +3/2*μm/a^(5/2) + 7/2*K*μm/a^(11/2)

    # Initial guess using a non-perturbed orbit.
    a_0 = (μm/n_0)^(2/3)

    # Newton-Raphson algorithm.
    #
    # Notice that we will allow, at most, 20 iterations.
    for k = 1:20
        res = f(a_0)
        a_0 = a_0 - res/df(a_0)

        (abs(res) < 1e-10) && break
    end

    # Auxiliary variables.
    dn  = dn_o2*2
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

        t_0, a_0, n_0, e_0, i_0, Ω_0, ω_0, M_0, dn_o2, ddn_o6, a_0, e_0, i_0,
        Ω_0, ω_0, M_0, n_0, f_0, C1, C2, C3, C4, j2_gc

    )
end

"""
### function j2!(j2d::J2_Structure{T}, t::Number) where T

Propagate the orbit defined in `j2d` until the time `t`. Notice that the values
in `j2d` will be modified.

##### Args

* j2d: J2 orbit propagator structure (see `J2_Structure`).
* t: Time in which the elements will be propagated [s].

##### Returns

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
    t_0, a_0, n_0, e_0, i_0, Ω_0, ω_0, M_0, dn_o2, ddn_o6, a_k, e_k, i_k, Ω_k,
    ω_k, M_k, n_k, f_k, C1, C2, C3, C4, j2_gc = j2d[:]

    R0, μm, J2 = j2_gc[:]

    # Time elapsed since epoch.
    Δt = t - t_0

    # Propagate the orbital elements.
    a_k = a_0 - C1*Δt
    e_k = e_0 - C2*Δt
    i_k = i_0
    Ω_k = mod(Ω_0 - C3*cos(i_k)*Δt,                    2*pi)
    ω_k = mod(ω_0 + C4*(4 - 5sin(i_k)^2)*Δt,           2*pi)
    M_k = mod(M_0 + n_0*Δt + dn_o2*Δt^2 + ddn_o6*Δt^3, 2*pi)
    f_k = M_to_f(e_k, M_k)

    # Make sure that eccentricity is not lower than 0.
    (e_k < 0) && (e_k = T(0))

    # Compute the position and velocity vectors given the orbital elements.
    (r_i_k, v_i_k) = kepler_to_rv(a_k*R0, e_k, i_k, Ω_k, ω_k, f_k)

    # Update the J2 orbit propagator structure.
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
