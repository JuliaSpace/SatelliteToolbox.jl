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
#   Functions created to compose rotations using Direction Cosine Matrices
#   (DCMs) and Quaternions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-05-09: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export compose_rotation

"""
### function compose_rotation(R1, R2, R3, ...) where T<:Number

Compose a set of rotations `R1`, `R2`, `R3`, ... such that:

    Result = R1 (First rotation) -> R2 (Second rotation) -> ...

The type of the rotations must be `Matrix` (DCM) or `Quaternion`.

##### Args

* R1, R2, R3, ...: Rotations (the first rotation is the first).

##### Returns

The full rotation. It will be represented on the same type as the input.

##### Remarks

The type of the elements in `R1`, `R2`, `R3`, ... **must** be the same.

"""

@inline function compose_rotation(D1::Matrix{T}, Ds::Matrix...) where T<:Number
    result = D1

    for Di in Ds
        result = Di*result
    end

    result
end

@inline function compose_rotation(q1::Quaternion{T}, qs::Quaternion...) where T<:Number
    result = q1

    for qi in qs
        result = result*qi
    end

    result
end
