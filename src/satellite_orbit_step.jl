#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    Propagate the satellite orbit.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References:
#
#    [1] Vallado, D. A., McClain, W. D (2013). Fundamentals of Astrodynamics
#        and Applications. Microcosm Press.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2014-12-08: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export satellite_orbit_compute_f

#==#
#
# @brief Compute the true anomaly given the mean anomaly.
#
# @param[in] a Semi-major axis [m].
# @param[in] e Eccentricity.
# @param[in] i Inclination [rad].
# @param[in] m Mean anomaly [rad].
# @param[in] tol (OPT) Tolerance for the Newton-Raphson method, default 1e-10.
#
# @retval SAT_LIGHTING_SUNLIGHT Satellite is under sunlight.
# @retval SAT_LIGHTING_PENUMBRA Satellite is at penumbra region.
# @retval SAT_LIGHTING_UMBRA Satellite is at umbra region.
#
# @note This function uses the Newton-Raphson method to compute the true
# anomaly.
#
#==#

function satellite_orbit_compute_f(a::Real,
                                   e::Real,
                                   i::Real,
                                   m::Real,
                                   tol::Float64 = 1e-10)

    # Update the eccentric anomaly using the Newton-Raphson method.
    # =============================================================

    # Initial guess.
    u = m

    sin_u = sin(u)
    cos_u = cos(u)

    # Newton-Raphson iterations.
    while ( abs(u - e*sin_u - m) > tol )
        u = u - (u - e*sin_u - m)/(1-e*cos_u)

        sin_u = sin(u)
        cos_u = cos(u)
    end

    # Transform u to the interval [0, 2*pi].
    u = mod(u, 2*pi)

    # Compute the true anomaly.
    half_f = atan( sqrt( (1+e)/(1-e) )*(sin_u/(1+cos_u))  )

    # f/2 and u/2 are in the same quadrant.
    if ( (u/2 > pi/2) && (u/2 <= 1.5*pi) )
        half_f += pi
    elseif ( u/2 > 1.5*pi )
        half_f = 2*pi - half_f
    end

    # Compute the true anomaly in the interval [0, 2*pi].
    f = mod(2*half_f, 2*pi)
end
