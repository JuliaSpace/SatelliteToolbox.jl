#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
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
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#   Microcosm Press, Hawthorn, CA, USA.
#
#   [2] Merson, R. H (1961). The motion of a satellite in an axi-symmetric
#   gravitational field. Geophysical journal international, Vol. 4(1), p. 17-52.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#
#
export j4_gc_egm08, j4_gc_egm96, j4_gc_jgm02, j4_gc_jgm03
export j4_init, j4!

################################################################################
#                                  Overloads
################################################################################

# Copy for J4_Structure
Base.copy(m::J4_Structure) = J4_Structure([ getfield(m, k) for k = 1:length(fieldnames(m)) ]...)

# Deepcopy for J4_Structure.
Base.deepcopy(m::J4_Structure) = J4_Structure([ deepcopy(getfield(m, k)) for k = 1:length(fieldnames(m)) ]...)

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
        sqrt(3.986004415e14/6378137.0^3),
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

# JGM-02 Gravitational constants.
const j4_gc_jgm02 = J4_GravCte(
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

################################################################################
#                                  Functions
################################################################################

"""
    j4_init(j4_gc::J4_GravCte{T}, epoch::Number, n_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, M_0::Number) where T

Initialize the data structure of J4 orbit propagator algorithm.

# Args

* `j4_gc`: J4 orbit propagator gravitational constants (see `J4_GravCte`).
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

The structure `J4_Structure` with the initialized parameters.

# Remarks

The inputs are the mean orbital elements.

"""
function j4_init(j4_gc::J4_GravCte{T}, epoch::Number, a_0::Number, e_0::Number,
                 i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number,
                 dn_o2::Number, ddn_o6::Number) where T

    # Unpack the gravitational constants to improve code readability.
    @unpack_J4_GravCte j4_gc

    # Initial values.
    al_0 = a_0/R0                 # Normalized semi-major axis [er].
    n_0  = μm/al_0^(3/2)          # Unperturbed mean motion [rad/s].
    p_0  = al_0*(1-e_0^2)         # Semi-latus rectum [er].
    M_0  = f_to_M(e_0, f_0)       # Initial mean anomaly [rad].

    # Auxiliary variables.
    dn               = 2dn_o2     # Time-derivative of the mean motion [rad/s²].
    e_0²             = e_0^2
    e_0⁴             = e_0^4
    sin_i_0, cos_i_0 = sincos(i_0)
    sin_i_0          = sin(i_0)
    sin_i_0²         = sin_i_0^2
    sin_i_0⁴         = sin_i_0^4
    aux              = (1-e_0²)
    saux             = sqrt(aux)
    p_0²             = p_0^2
    p_0⁴             = p_0^4

    # First-order time-derivative of the orbital elements.
    #
    # See [1, p 692].

    δa   = -2/3*al_0*dn/n_0
    δe   = -2/3*(1-e_0)*dn/n_0
    δΩ   = -3/2 *n_0*J2  /p_0²*cos_i_0 +
            3/32*n_0*J2^2/p_0⁴*cos_i_0*(12 -  4e_0² - (80 +  5e_0²)*sin_i_0²) +
           15/32*n_0*J4  /p_0⁴*cos_i_0*( 8 + 12e_0² - (14 + 21e_0²)*sin_i_0²)

    δω   =  3/  4*n_0*J2  /p_0²*(4 - 5sin_i_0²) +
            9/384*n_0*J2^2/p_0⁴*(     56e_0² + (760 -  36e_0²)*sin_i_0² - (890 +  45*e_0²)*sin_i_0⁴) -
           15/128*n_0*J4  /p_0⁴*(64 + 72e_0² - (248 + 252e_0²)*sin_i_0² + (196 + 189*e_0²)*sin_i_0⁴)

    δM_0 =  3/  4*n_0*J2  /p_0²*sqrt(1-e_0²)*(2 - 3sin_i_0²) +
            3/512*n_0*J2^2/p_0⁴/sqrt(1-e_0²)*((320e_0² - 280e_0^4) +
                                              ( 1600 - 1568e_0² + 328e_0^4)*sin_i_0² +
                                              (-2096 + 1072e_0² +  79e_0^4)*sin_i_0⁴) -
           45/128*n_0*J4  /p_0⁴*sqrt(1-e_0²)*e_0²*(-8 + 40sin_i_0 - 35sin_i_0²)

    # Create the output structure with the data.
    J4_Structure{T}(
        epoch, al_0, n_0, e_0, i_0, Ω_0, ω_0, f_0, M_0, 0, dn_o2, ddn_o6, al_0,
        e_0, i_0, Ω_0, ω_0, f_0, M_0, δa, δe, δΩ, δω, δM_0, j4_gc
       )
end

"""
    j4!(j4d::J4_Structure{T}, t::Number) where T

Propagate the orbit defined in `j4d` (see `J4_Structure`) until the time `t`
[s]. Notice that the values in `j4d` will be modified.

# Returns

* The position vector represented in the inertial frame at time `t` [m].
* The velocity vector represented in the inertial frame at time `t` [m/s]

###### Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. If the orbit parameters are obtained
from a TLE, then the inertial frame will be TEME. Notice, however, that the
perturbation theory requires an inertial frame with true equator.

"""
function j4!(j4d::J4_Structure{T}, t::Number) where T
    # Unpack the variables.
    @unpack_J4_Structure j4d
    @unpack_J4_GravCte   j4_gc

    # Time elapsed since epoch.
    Δt = t

    # Propagate the orbital elements.
    al_k = al_0 + δa*Δt
    e_k  = e_0  + δe*Δt
    i_k  = i_0
    Ω_k  = mod(Ω_0 + δΩ*Δt,         2π)
    ω_k  = mod(ω_0 + δω*Δt,         2π)
    M_k  = mod(M_0 + (δM_0+n_0)*Δt, 2π)
    f_k  = M_to_f(e_k, M_k)

    # Make sure that eccentricity is not lower than 0.
    (e_k < 0) && (e_k = T(0))

    # Compute the position and velocity vectors given the orbital elements.
    (r_i_k, v_i_k) = kepler_to_rv(al_k*R0, e_k, i_k, Ω_k, ω_k, f_k)

    # Update the J4 orbit propagator structure.
    j4d.Δt   = Δt
    j4d.al_k = al_k
    j4d.e_k  = e_k
    j4d.i_k  = i_k
    j4d.Ω_k  = Ω_k
    j4d.ω_k  = ω_k
    j4d.M_k  = M_k
    j4d.f_k  = f_k

    # Return the position and velocity vector represented in the inertial
    # reference frame.
    (r_i_k, v_i_k)
end
