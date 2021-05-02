# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   J2 osculating orbit propagator algorithm.
#
#       This algorithm propagates the orbit considering the secular and
#       short-period perturbations introduced by the J2 gravitational term. The
#       algorithm is based on Kwok version as indicated in [1, p. 708-710].
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

export j2osc_init, j2osc!

################################################################################
#                                  Functions
################################################################################

"""
    j2osc_init(epoch::Number, a_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number, dn_o2::Number, ddn_o6::Number; j2_gc::J2_GravCte{T} = j2_gc_egm08) where T

Initialize the data structure of J2 osculating orbit propagator algorithm.

# Args

* `epoch`: Epoch of the orbital elements [Julian Day].
* `a_0`: Initial semi-major axis [m].
* `e_0`: Initial eccentricity.
* `i_0`: Initial inclination [rad].
* `Ω_0`: Initial right ascension of the ascending node [rad].
* `ω_0`: Initial argument of perigee [rad].
* `f_0`: Initial true anomaly [rad].
* `dn_o2`: First time derivative of the mean motion divided by two [rad/s^2].
* `ddn_o6`: Second time derivative of the mean motion divided by six [rad/s^3].

# Keywords

* `j2_gc`: J2 orbit propagator gravitational constants (see `J2_GravCte`).
           (**Default** = `j2_gc_egm08`)

# Returns

The structure `J2osc_Structure` with the initialized parameters.

# Remarks

The inputs are the mean orbital elements.

"""
function j2osc_init(epoch::Number, a_0::Number, e_0::Number, i_0::Number,
                    Ω_0::Number, ω_0::Number, f_0::Number, dn_o2::Number,
                    ddn_o6::Number; j2_gc::J2_GravCte{T} = j2_gc_egm08) where T

    j2d = j2_init(epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, dn_o2, ddn_o6;
                  j2_gc = j2_gc)

    # Initialize the structure.
    j2oscd = J2osc_Structure{T}(j2d, epoch, 0, 0, 0, 0, 0, 0, 0)

    # Call the propagation one first time to update the osculating elements.
    j2osc!(j2oscd, 0)

    return j2oscd
end

"""
    j2osc!(j2d::J2osc_Structure{T}, t::Number) where T

Propagate the orbit defined in `j2oscd` (see `J2osc_Structure`) until the time
`t` [s]. Notice that the values in `j2oscd` will be modified.

# Returns

* The position vector represented in the inertial frame at time `t` [m].
* The velocity vector represented in the inertial frame at time `t` [m/s]

# Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. Notice, however, that the
perturbation theory requires an inertial frame with true equator.

"""
function j2osc!(j2oscd::J2osc_Structure{T}, t::Number) where T
    # First, we need to propagate the mean elements since they are necessary to
    # compute the short-periodic perturbations.
    j2d = j2oscd.j2d
    j2!(j2d, t)

    @unpack_J2_GravCte j2d.j2_gc

    # Obtain the mean elements at this time instant.
    al_k = j2d.al_k
    e_k  = j2d.e_k
    i_k  = j2d.i_k
    Ω_k  = j2d.Ω_k
    ω_k  = j2d.ω_k
    M_k  = j2d.M_k
    f_k  = j2d.f_k
    p_k  = al_k*R0*(1 - e_k^2)
    u_k  = ω_k + f_k

    # Auxiliary variables to reduce the computational burden.
    KJ2 = J2*R0*R0
    sin_i_k, cos_i_k = sincos(i_k)
    sin_f_k, cos_f_k = sincos(f_k)
    sin_2u_k, cos_2u_k = sincos(2u_k)
    sin_2ω_f_k, cos_2ω_f_k = sincos(2ω_k + f_k)
    sin_2ω_3f_k, cos_2ω_3f_k = sincos(2ω_k + 3f_k)
    sin_i_k² = sin_i_k * sin_i_k
    cos_i_k² = cos_i_k * cos_i_k
    p_k² = p_k * p_k
    e_k² = e_k * e_k
    e_cos_f_k = e_k * cos_f_k
    e_sin_f_k = e_k * sin_f_k

    aux1 = 3cos(2u_k) + 3e_k*cos_2ω_f_k + e_k*cos_2ω_3f_k
    aux2 = sqrt(1 - e_k²)
    aux3 = 3cos_i_k² - 1

    # Compute the short-periodic perturbations considering only the J2
    # gravitational term.
    δisp_k = +KJ2*sin_i_k*cos_i_k/(4p_k²)*aux1

    δpsp_k = +KJ2*sin_i_k²/(2p_k)*aux1

    δΩsp_k = -KJ2*cos_i_k/(4p_k²)*(
        6*(f_k - M_k + e_sin_f_k) - 2sin_2u_k - 3e_k*sin_2ω_f_k - e_k*sin_2ω_3f_k
    )

    δrsp_k = -KJ2/(4p_k)*(
        aux3*( 2*aux2/(1 + e_cos_f_k) + e_cos_f_k/(1 + aux2) ) -
        sin_i_k²*cos_2u_k
    )

    δṙsp_k = +KJ2*sqrt(μm)/(4p_k^2.5)*(
        aux3*e_sin_f_k*( aux2 + (1 + e_cos_f_k)^2/(1 + aux2) ) -
        sin_i_k²*(1 - e_cos_f_k)^2*sin_2u_k
    )

    δusp_k = +KJ2/(8p_k²)*(
        (6 - 30cos_i_k²)*(f_k - M_k) +
        4e_k*sin_f_k*((1 - 6cos_i_k²) - aux3/(1 + aux2)) -
        aux3/(1 + aux2)*e_k²*sin(2f_k) +
        (5cos_i_k² - 2)*2e_k*sin_2ω_f_k +
        (7cos_i_k² - 1)*sin_2u_k + 2*cos_i_k²*e_k*sin_2ω_3f_k
    )

    r_k = p_k/(1 + e_k*cos(f_k))
    ṙ_k = sqrt(μm/p_k)*e_k*cos(f_k)

    r_osc_k = r_k + δrsp_k
    ṙ_osc_k = ṙ_k + δṙsp_k
    p_osc_k = p_k + δpsp_k

    A_k = p_osc_k/r_osc_k - 1
    B_k = sqrt(p_osc_k/μm)*ṙ_osc_k

    e_osc_k = sqrt(A_k^2 + B_k^2)
    a_osc_k = p_osc_k/(1-e_osc_k^2)
    i_osc_k = i_k + δisp_k
    Ω_osc_k = Ω_k + δΩsp_k
    u_osc_k = u_k + δusp_k
    f_osc_k = atan(B_k, A_k)
    ω_osc_k = u_osc_k - f_osc_k
    M_osc_k = f_to_M(e_osc_k, f_osc_k)

    # Compute the position and velocity considering the osculating elements.
    r_i_k, v_i_k = kepler_to_rv(a_osc_k, e_osc_k, i_osc_k, Ω_osc_k, ω_osc_k, f_osc_k)

    # Update the J2 orbit propagator structure.
    j2oscd.Δt  = t
    j2oscd.a_k = a_osc_k
    j2oscd.e_k = e_osc_k
    j2oscd.i_k = i_osc_k
    j2oscd.Ω_k = Ω_osc_k
    j2oscd.ω_k = ω_osc_k
    j2oscd.M_k = M_osc_k
    j2oscd.f_k = f_osc_k

    return r_i_k, v_i_k
end
