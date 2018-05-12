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
#   Compute the lighting condition of a satellite:
#       1) Sunlight,
#       2) penumbra, or
#       3) umbra.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Longo, C. R. O., Rickman, S. L (1995). Method for the Calculation of
#       Spacecraft Umbra and Penumbra Shadow Terminator Points. NASA Technical
#       Paper 3547.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2014-07-28: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export SAT_LIGHTING_SUNLIGHT, SAT_LIGHTING_UMBRA, SAT_LIGHTING_PENUMBRA
export satellite_lighting_condition

# Constants.
const SAT_LIGHTING_SUNLIGHT = 0
const SAT_LIGHTING_PENUMBRA = 1
const SAT_LIGHTING_UMBRA = 2

"""
### function satellite_lighting_condition(r_i::AbstractVector, s_i::AbstractVector)

Compute the satellite lighting condition given the Sun unitary vector `s_i` and
the satellite position vector `r_i`.

##### Args

* r_i: Satellite position vector represented in the Inertial coordinate frame
       [m].
* s_i: Sun position vector represented in the Inertial coordinate frame [m].

##### Returns

* `SAT_LIGHTING_SUNLIGHT`: Satellite is under sunlight.
* `SAT_LIGHTING_PENUMBRA`: Satellite is at penumbra region.
* `SAT_LIGHTING_UMBRA`: Satellite is at umbra region.

"""

function satellite_lighting_condition(r_i::AbstractVector, s_i::AbstractVector)
    # Norme of the sun position vector.
    norm_s_i = norm(s_i)

    # Umbra section [1].
    Xu = R0*norm_s_i/(Rs-R0)
    Au = asin(R0/Xu)

    # Penumbra section [1].
    Xp = R0*norm_s_i/(Rs+R0)
    Ap = asin(R0/Xp)

    # Projection of the satellite position vector on the sun direction [1].
    rs_i = dot(r_i, s_i)*s_i/(norm_s_i^2)

    # Distance of the umbral cone and the spacecraft [1].
    delta_i = r_i - rs_i

    # Location of the umbral cone terminator at the projected spacecraft
    # location [1].
    ep_i = (Xu - norm(rs_i))*tan(Au)

    # Location of the penumbral cone terminator at the projected spacecraft
    # location [1].
    Kp_i = (Xp + norm(rs_i))*tan(Ap)

    # Check if the satellite is under sunlight or umbra / penumbra.
    if dot(rs_i, s_i/norm_s_i) < 0
        if norm(delta_i) > Kp_i
            return SAT_LIGHTING_SUNLIGHT
        else
            if norm(delta_i) < ep_i
                return SAT_LIGHTING_UMBRA
            else
                return SAT_LIGHTING_PENUMBRA
            end
        end
    else
        return SAT_LIGHTING_SUNLIGHT
    end
end
