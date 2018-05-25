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
#   Two Body orbit propagator.
#
#       This algorithm considers a perfect Keplerian orbit. In other words, no
#       perturbation is considered during the propagation and the Earth is
#       modeled as a perfect sphere.
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
# 2018-03-30: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export TwoBody_Structure
export twobody_init, twobody!

################################################################################
#                                  Overloads
################################################################################

# Copy for TwoBody_Structure.
function Base.copy(m::TwoBody_Structure)
    TwoBody_Structure(m.t_0,
                      m.n_0,
                      m.e_0,
                      m.i_0,
                      m.Ω_0,
                      m.ω_0,
                      m.M_0,
                      m.a,
                      m.M_k,
                      m.f_k,
                      m.μ)
end

################################################################################
#                                  Functions
################################################################################

"""
    function twobody_init(t_0::Number, n_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, M_0::Number, μ::T) where T

Initialize the data structure of Two Body orbit propagator algorithm.

# Args

* `t_0`: Epoch of the orbital elements [s].
* `n_0`: Mean motion at epoch [rad/s].
* `e_0`: Initial eccentricity.
* `i_0`: Initial inclination [rad].
* `Ω_0`: Initial right ascension of the ascending node [rad].
* `ω_0`: Initial argument of perigee [rad].
* `M_0`: Initial mean anomaly.
* `μ`: Standard gravitational parameter of the central body [m^3/s^2].

# Returns

The structure `TwoBody_Structure` with the initialized parameters.

"""
function twobody_init(t_0::Number,
                      n_0::Number,
                      e_0::Number,
                      i_0::Number,
                      Ω_0::Number,
                      ω_0::Number,
                      M_0::Number,
                      μ::T) where T
    # The propagator is only defined for 0 <= e < 1.
    if !(0 <= e_0 < 1)
        throw(ArgumentError("The Two Body propagator only supports eccentricities in the interval [0,1)"))
    end

    # Compute the semi-major axis using the angular velocity.
    a = (μ/n_0^2)^(1/3)

    # Compute the true anomaly.
    f = M_to_f(e_0, M_0)

    # Create and return the Two Body orbit propagator structure.
    TwoBody_Structure{T}(t_0, n_0, e_0, i_0, Ω_0, ω_0, M_0, a, M_0, f, μ)
end

"""
    function twobody!(tbd::TwoBody_Structure, t::Number)

Propagate the orbit defined in `tbd` until the time `t`. Notice that the values
in `tbd` will be modified.

# Args

* `tbd`: Two Body orbit propagator structure (see `TwoBody_Structure`).
* `t`: Time in which the elements will be propagated [s].

# Returns

* The position vector represented in the inertial frame at time `t` [m].
* The velocity vector represented in the inertial frame at time `t` [m/s]

###### Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. If the orbit parameters are obtained
from a TLE, then the inertial frame will be TEME.

"""
function twobody!(tbd::TwoBody_Structure, t::Number)
    # Time elapsed since epoch.
    Δt = t - tbd.t_0

    # Update the mean anomaly.
    tbd.M_k = tbd.M_0 + tbd.n_0*Δt

    # Convert the mean anomaly to true anomaly.
    tbd.f_k = M_to_f(tbd.e_0, tbd.M_k)

    # Compute the position and velocity vectors given the orbital elements.
    (r_i_k, v_i_k) =
        kepler_to_rv(tbd.a, tbd.e_0, tbd.i_0, tbd.Ω_0, tbd.ω_0, tbd.f_k)

    # Return the position and velocity vector represented in the inertial
    # reference frame.
    (r_i_k, v_i_k)
end
