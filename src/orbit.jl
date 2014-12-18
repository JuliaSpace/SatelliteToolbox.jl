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
#    Many auxiliary functions for orbit computations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2014-12-18: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

#==#
# 
# @brief Return perturbation of the right accession of the ascending node.
#
# @param[in] a Semi-major axis [m].
# @param[in] e Eccentricity.
# @param[in] i Inclination [rad].
#
# @return The perturbation of the RAAN [rad/s].
#
#==#

function dRAAN_J2(a::Real, e::Real, i::Real)
    # Semi-lactum rectum.
    p = a*(1-e^2)

    # Unperturbed orbit period.
    n0 = n_J0(a)
    
    # Perturbation of the right ascension of the ascending node.
    -3.0/2.0*R0^2/(p^2)*n0*J2*cos(i)
end

#==#
# 
# @brief Return perturbation of the argument of perigee.
#
# @param[in] a Semi-major axis [m].
# @param[in] e Eccentricity.
# @param[in] i Inclination [rad].
#
# @return The perturbation of the argument of perigee [rad/s].
#
#==#

function dw_J2(a::Real, e::Real, i::Real)
    # Semi-lactum rectum.
    p = a*(1-e^2)

    # Unperturbed orbit period.
    n0 = n_J0(a)
    
    # Perturbation of the argument of perigee.
    3*R0^2*J2/(4*p^2)*n0*(5*cos(i)^2-1)
end

#==#
# 
# @brief Return the orbit angular velocity neglecting the perturbations
# (two-body).
#
# @param[in] a Semi-major axis [m].
#
# @return The unperturbed orbit angular velocity [rad/s].
#
#==#

n_J0(a::Real) = sqrt(m0/a^3)

#==#
# 
# @brief Return the orbit angular velocity considering the perturbations (up to
# J2 terms).
#
# @param[in] a Semi-major axis [m].
# @param[in] e Eccentricity.
# @param[in] i Inclination [rad].
#
# @return The perturbed orbit angular velocity [rad/s].
#
#==#

function n_J2(a::Real, e::Real, i::Real)
    # Semi-lactum rectum.
    p = a*(1-e^2)

    # Unperturbed orbit period.
    n0 = n_J0(a)

    # Orbit period considering the perturbations (up to J2).
    n0 + 3*R0^2*J2/(4*p^2)*n0*(sqrt(1-e^2)*(3*cos(i)^2-1) + (5*cos(i)^2-1))
end
