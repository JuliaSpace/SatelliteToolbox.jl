#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
#   Rotations from an Earth-Centered Inertial (ECI) reference frame to another
#   ECI reference frame.
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
# 2018-05-25: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export rECItoECI

"""
    function rECEFtoECI([T,] [M,] ECIo, ECIf, JD_UTC::Number [, eop_data])
    function rECEFtoECI([T,] [M,] ECIo, JD_UTCo::Number, ECIf, JD_UTCf::Number [, eop_data])

Compute the rotation from an Earth-Centered Inertial (`ECI`) reference frame to
another ECI reference frame. If the origin and destination frame contain only
one *of date* frame, then the first signature is used and `JD_UTC` is the epoch
of this frame. On the other hand, if the origin and destination frame contain
two *of date* frame[^1], e.g. TOD => MOD, then the second signature must be used
in which `JD_UTCo` is the epoch of the origin frame and `JD_UTCf` is the epoch
of the destination frame.

The rotation description that will be used is given by `T`, which can be `DCM`
or `Quaternion`. The model used to compute the rotation is specified by `M`.
Currently, only IAU-76/FK5 is supported (`M = FK5()`). The origin ECI frame is
selected by the input `ECIo` and the destination ECI frame is selected by the
input `ECIf`.

[^1]: TEME is an *of date* frame.

# Rotation description

The rotations that aligns the origin ECI frame with the destination ECI frame
can be described by Direction Cosine Matrices or Quaternions. This is selected
by the parameter `T`.

The possible values are:

* `DCM`: The rotation will be described by a Direction Cosine Matrix.
* `Quaternion`: The rotation will be described by a Quaternion.

If no value is specified, then it falls back to `DCM`.

# Conversion model

The model that will be used to compute the rotation is given by `M`. The
possible values are:

* `FK5()`: Use the IAU-76/FK5 model.

If no value is specified, then it falls back to `FK5()`.

# ECI Frame

The supported ECI frames for both origin `ECIo` and destination `ECIf` are:

* `TEME()`: ECI will be selected as the True Equator Mean Equinox (TEME)
            reference frame.
* `TOD()`: ECI will be selected as the True of Date (TOD).
* `MOD()`: ECI will be selected as the Mean of Date (MOD).
* `J2000()`: ECI will be selected as the J2000 reference frame.
* `GCRF()`: ECI will be selected as the Geocentric Celestial Reference Frame
            (GCRF).

# EOP Data

The conversion between the frames depends on EOP Data (see `get_iers_eop` and
`read_iers_eop`). If IAU-76/FK5 model is used, then the type of `eop_data` must
be `EOPData_IAU1980`. The following table shows the requirements for EOP data
given the selected frames.

|   Model    |   ECIo  |   ECIf  |    EOP Data   | Function Signature |
|:-----------|:--------|:--------|:--------------|:-------------------|
| IAU-76/FK5 | `GCRF`  | `J2000` | EOP IAU1980   | First              |
| IAU-76/FK5 | `GCRF`  | `MOD`   | EOP IAU1980   | First              |
| IAU-76/FK5 | `GCRF`  | `TOD`   | EOP IAU1980   | First              |
| IAU-76/FK5 | `GCRF`  | `TEME`  | EOP IAU1980   | First              |
| IAU-76/FK5 | `J2000` | `GCRF`  | EOP IAU1980   | First              |
| IAU-76/FK5 | `J2000` | `MOD`   | EOP IAU1980   | First              |
| IAU-76/FK5 | `J2000` | `TOD`   | EOP IAU1980   | First              |
| IAU-76/FK5 | `J2000` | `TEME`  | Not required  | First              |
| IAU-76/FK5 | `MOD`   | `GCRF`  | EOP IAU1980   | First              |
| IAU-76/FK5 | `MOD`   | `J2000` | EOP IAU1980   | First              |
| IAU-76/FK5 | `MOD`   | `TOD`   | EOP IAU1980   | Second             |
| IAU-76/FK5 | `MOD`   | `TEME`  | EOP IAU1980   | Second             |
| IAU-76/FK5 | `TOD`   | `GCRF`  | EOP IAU1980   | First              |
| IAU-76/FK5 | `TOD`   | `J2000` | EOP IAU1980   | First              |
| IAU-76/FK5 | `TOD`   | `MOD`   | EOP IAU1980   | Second             |
| IAU-76/FK5 | `TOD`   | `TEME`  | EOP IAU1980   | Second             |
| IAU-76/FK5 | `TEME`  | `GCRF`  | EOP IAU1980   | First              |
| IAU-76/FK5 | `TEME`  | `J2000` | Not requrired | First              |
| IAU-76/FK5 | `TEME`  | `MOD`   | EOP IAU1980   | Second             |
| IAU-76/FK5 | `TEME`  | `TOD`   | EOP IAU1980   | Second             |

## MOD and TOD

In this function, MOD and TOD frames are always defined with IERS EOP
corrections. Hence, if one wants to obtain the MOD and TOD frames according to
the original IAU-76/FK5 theory, it is necessary to use the low-level functions
in file `./src/transformations/fk5/fk5.jl`.

# Args

* `T`: (OPTIONAL) Type of the rotation representation (**Default** = `DCM`).
* `M`: (OPTIONAL) Model used to compute the rotation (**Default** = `FK5()`).
* `ECIo`: Origin ECEF frame.
* `ECIf`: Destination ECI frame.
* `JD_UTC`: Julian day [UTC].
* `eop_data`: EOP Data.

* `JD_UTCo`: Julian day of the origin reference frame [UTC].
* `JD_UTCf`: Julian day of the destination reference frame [UTC].

# Returns

The rotation description represented by `T` that rotates the origin ECI
reference frame into alignment with the destination ECI reference frame.

# Examples

```julia-repl
julia> eop_IAU1980 = get_iers_eop(:IAU1980)
julia> rECItoECI(DCM, GCRF(), J2000(), DatetoJD(1986, 6, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
  1.0          -1.94091e-12   5.56251e-10
  1.94082e-12   1.0          -1.45796e-9
 -5.56251e-10   1.45796e-9    1.0

julia> rECItoECI(Quaternion, TEME(), GCRF(), DatetoJD(1986, 6, 19, 21, 35, 0), eop_IAU1980)
Quaternion{Float64}:
  + 0.9999986335698354 + 1.8300220800378406e-5.i + 0.0006653038777338896.j - 0.0015132396750687136.k

julia> rECItoECI(TOD(), DatetoJD(1986,6,19,21,35,0), TOD(), DatetoJD(1987,5,19,3,0,0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 1.0          -0.000224084  -9.7377e-5
 0.000224083   1.0          -5.80249e-6
 9.73783e-5    5.78067e-6    1.0

julia> rECItoECI(Quaternion, TOD(), JD_J2000, MOD(), JD_J2000, eop_IAU1980)
Quaternion{Float64}:
  + 0.9999999993282446 - 1.4002278413388636e-5.i + 1.3473863722202666e-5.j - 3.107896586417327e-5.k

julia> rECItoECI(J2000(), TEME(), DatetoJD(1986,6,19,21,35,0))
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
  0.999995    0.0030265    0.00133055
 -0.00302645  0.999995    -3.86125e-5
 -0.00133066  3.45854e-5   0.999999
```
"""
@inline rECItoECI(T_ECIo::T_ECIs,
                  T_ECIf::T_ECIs,
                  JD_UTC::Number,
                  eop_data::EOPData_IAU1980) =
    rECItoECI(DCM, Val{:FK5}, T_ECIo, T_ECIf, JD_UTC, eop_data)

@inline rECItoECI(T::Union{Type{DCM}, Type{Quaternion}},
                  T_ECIo::T_ECIs,
                  T_ECIf::T_ECIs,
                  JD_UTC::Number,
                  eop_data::EOPData_IAU1980) =
    rECItoECI(T, Val{:FK5}, T_ECIo, T_ECIf, JD_UTC, eop_data)

@inline rECItoECI(M::Type{Val{:FK5}},
                  T_ECIo::T_ECIs,
                  T_ECIf::T_ECIs,
                  JD_UTC::Number,
                  eop_data::EOPData_IAU1980) =
    rECItoECI(DCM, M, T_ECIo, T_ECIf, JD_UTC, eop_data)

# Specializations for those cases in which we have two *of dates* frames.
@inline rECItoECI(T_ECIo::T_ECIs_of_date,
                  JD_UTCo::Number,
                  T_ECIf::T_ECIs_of_date,
                  JD_UTCf::Number,
                  eop_data::EOPData_IAU1980) =
    rECItoECI(DCM, Val{:FK5}, T_ECIo, JD_UTCo, T_ECIf, JD_UTCf, eop_data)

@inline rECItoECI(T::Union{Type{DCM}, Type{Quaternion}},
                  T_ECIo::T_ECIs_of_date,
                  JD_UTCo::Number,
                  T_ECIf::T_ECIs_of_date,
                  JD_UTCf::Number,
                  eop_data::EOPData_IAU1980) =
    rECItoECI(T, Val{:FK5}, T_ECIo, JD_UTCo, T_ECIf, JD_UTCf, eop_data)

@inline rECItoECI(M::Type{Val{:FK5}},
                  T_ECIo::T_ECIs_of_date,
                  JD_UTCo::Number,
                  T_ECIf::T_ECIs_of_date,
                  JD_UTCf::Number,
                  eop_data::EOPData_IAU1980) =
    rECItoECI(DCM, M, T_ECIo, JD_UTCo, T_ECIf, JD_UTCf, eop_data)

# Specializations for those cases that EOP Data is not needed.
@inline rECItoECI(T_ECIo::Union{Type{Val{:J2000}},Type{Val{:TEME}}},
                  T_ECIf::Union{Type{Val{:J2000}},Type{Val{:TEME}}},
                  JD_UTC::Number) =
    rECItoECI(DCM, Val{:FK5}, T_ECIo, T_ECIf, JD_UTC)

@inline rECItoECI(T::Union{Type{DCM}, Type{Quaternion}},
                  T_ECIo::Union{Type{Val{:J2000}},Type{Val{:TEME}}},
                  T_ECIf::Union{Type{Val{:J2000}},Type{Val{:TEME}}},
                  JD_UTC::Number) =
    rECItoECI(T, Val{:FK5}, T_ECIo, T_ECIf, JD_UTC)

@inline rECItoECI(M::Type{Val{:FK5}},
                  T_ECIo::Union{Type{Val{:J2000}},Type{Val{:TEME}}},
                  T_ECIf::Union{Type{Val{:J2000}},Type{Val{:TEME}}},
                  JD_UTC::Number) =
    rECItoECI(DCM, M, T_ECIo, T_ECIf, JD_UTC)

################################################################################
#                                  IAU-76/FK5
################################################################################

#                                GCRF <=> J2000
# ==============================================================================

function rECItoECI(T::Type,
                   ::Type{Val{:FK5}},
                   ::Type{Val{:GCRF}},
                   ::Type{Val{:J2000}},
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980)

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    δΔϵ_1980 = eop_data.dEps[JD_UTC]*pi/648000
    δΔψ_1980 = eop_data.dPsi[JD_UTC]*pi/648000

    # In this case, we need to convert GCRF back to PEF and then convert to
    # J2000, which is the same conversion from PEF to GCRF **without** the EOP
    # data.
    #
    # TODO: Can I simplify the rotation from TOD with corrections to TOD without
    # corrections?
    r_MOD_GCRF  = rGCRFtoMOD_fk5(T, JD_TT)
    r_PEF_MOD   = rMODtoPEF_fk5(T, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)
    r_MOD_PEF   = rPEFtoMOD_fk5(T, JD_UT1, JD_TT, 0, 0)
    r_J2000_MOD = inv_rotation(r_MOD_GCRF)

    # Compose the full rotation.
    compose_rotation(r_MOD_GCRF, r_PEF_MOD, r_MOD_PEF, r_J2000_MOD)
end

@inline rECItoECI(T::Type,
                  M::Type{Val{:FK5}},
                  T_ECIo::Type{Val{:J2000}},
                  T_ECIf::Type{Val{:GCRF}},
                  JD_UTC::Number,
                  eop_data::EOPData_IAU1980) =
    inv_rotation(rECItoECI(T, M, T_ECIf, T_ECIo, JD_UTC, eop_data))

#                                 GCRF <=> MOD
# ==============================================================================

function rECItoECI(T::Type,
                   ::Type{Val{:FK5}},
                   ::Type{Val{:GCRF}},
                   ::Type{Val{:MOD}},
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980)

    # Get the time in TT.
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Return the rotation.
    rGCRFtoMOD_fk5(T, JD_TT)
end

@inline rECItoECI(T::Type,
                  M::Type{Val{:FK5}},
                  T_ECIo::Type{Val{:MOD}},
                  T_ECIf::Type{Val{:GCRF}},
                  JD_UTC::Number,
                  eop_data::EOPData_IAU1980) =
    inv_rotation(rECItoECI(T, M, T_ECIf, T_ECIo, JD_UTC, eop_data))

#                                 GCRF <=> TOD
# ==============================================================================

function rECItoECI(T::Type,
                   ::Type{Val{:FK5}},
                   ::Type{Val{:GCRF}},
                   ::Type{Val{:TOD}},
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980)

    # Get the time in TT.
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    δΔϵ_1980 = eop_data.dEps[JD_UTC]*pi/648000
    δΔψ_1980 = eop_data.dPsi[JD_UTC]*pi/648000

    # Return the rotation.
    r_MOD_GCRF = rGCRFtoMOD_fk5(T, JD_TT)
    r_TOD_MOD  = rMODtoTOD_fk5(T, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Compose the full rotation.
    compose_rotation(r_MOD_GCRF, r_TOD_MOD)
end

@inline rECItoECI(T::Type,
                  M::Type{Val{:FK5}},
                  T_ECIo::Type{Val{:TOD}},
                  T_ECIf::Type{Val{:GCRF}},
                  JD_UTC::Number,
                  eop_data::EOPData_IAU1980) =
    inv_rotation(rECItoECI(T, M, T_ECIf, T_ECIo, JD_UTC, eop_data))

#                                GCRF <=> TEME
# ==============================================================================

function rECItoECI(T::Type,
                   ::Type{Val{:FK5}},
                   ::Type{Val{:GCRF}},
                   ::Type{Val{:TEME}},
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980)

    # Get the time in TT.
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    δΔϵ_1980 = eop_data.dEps[JD_UTC]*pi/648000
    δΔψ_1980 = eop_data.dPsi[JD_UTC]*pi/648000

    # Return the rotation.
    r_MOD_GCRF = rGCRFtoMOD_fk5(T, JD_TT)
    r_TEME_MOD = rMODtoTEME(T, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Compose the full rotation.
    compose_rotation(r_MOD_GCRF, r_TEME_MOD)
end

@inline rECItoECI(T::Type,
                  M::Type{Val{:FK5}},
                  T_ECIo::Type{Val{:TEME}},
                  T_ECIf::Type{Val{:GCRF}},
                  JD_UTC::Number,
                  eop_data::EOPData_IAU1980) =
    inv_rotation(rECItoECI(T, M, T_ECIf, T_ECIo, JD_UTC, eop_data))

#                                J2000 <=> MOD
# ==============================================================================

function rECItoECI(T::Type,
                   ::Type{Val{:FK5}},
                   ::Type{Val{:J2000}},
                   ::Type{Val{:MOD}},
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980)

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    δΔϵ_1980 = eop_data.dEps[JD_UTC]*pi/648000
    δΔψ_1980 = eop_data.dPsi[JD_UTC]*pi/648000

    # In this case, we need to convert J2000 back to PEF and then convert to
    # MOD. This is necessary because we need to apply EOP corrections to convert
    # to MOD and just a `rGCRFtoMOD_fk5` would yield a frame according to the
    # original IAU-76/FK5 theory.
    #
    # TODO: Can I simplify this rotation?
    r_MOD_J2000 = rGCRFtoMOD_fk5(T, JD_TT)
    r_PEF_MOD   = rMODtoPEF_fk5(T, JD_UT1, JD_TT, 0, 0)
    r_MOD_PEF   = rPEFtoMOD_fk5(T, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Compose the full rotation.
    compose_rotation(r_MOD_J2000, r_PEF_MOD, r_MOD_PEF)
end

@inline rECItoECI(T::Type,
                  M::Type{Val{:FK5}},
                  T_ECIo::Type{Val{:MOD}},
                  T_ECIf::Type{Val{:J2000}},
                  JD_UTC::Number,
                  eop_data::EOPData_IAU1980) =
    inv_rotation(rECItoECI(T, M, T_ECIf, T_ECIo, JD_UTC, eop_data))

#                                J2000 <=> TOD
# ==============================================================================

function rECItoECI(T::Type,
                   ::Type{Val{:FK5}},
                   ::Type{Val{:J2000}},
                   ::Type{Val{:TOD}},
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980)

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    δΔϵ_1980 = eop_data.dEps[JD_UTC]*pi/648000
    δΔψ_1980 = eop_data.dPsi[JD_UTC]*pi/648000

    # In this case, we need to convert J2000 back to PEF and then convert to
    # TOD. This is necessary because we need to apply EOP corrections to convert
    # to TOD.
    #
    # TODO: Can I simplify this rotation?
    r_TOD_J2000 = rGCRFtoMOD_fk5(T, JD_TT)
    r_PEF_MOD   = rMODtoPEF_fk5(T, JD_UT1, JD_TT, 0, 0)
    r_MOD_PEF   = rPEFtoMOD_fk5(T, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)
    r_TOD_MOD   = rMODtoTOD_fk5(T, JD_UT1, δΔϵ_1980, δΔψ_1980)

    # Compose the full rotation.
    compose_rotation(r_TOD_J2000, r_PEF_MOD, r_MOD_PEF, r_TOD_MOD)
end

@inline rECItoECI(T::Type,
                  M::Type{Val{:FK5}},
                  T_ECIo::Type{Val{:TOD}},
                  T_ECIf::Type{Val{:J2000}},
                  JD_UTC::Number,
                  eop_data::EOPData_IAU1980) =
    inv_rotation(rECItoECI(T, M, T_ECIf, T_ECIo, JD_UTC, eop_data))

#                                J2000 <=> TEME
# ==============================================================================

@inline rECItoECI(T::Type,
                  M::Type{Val{:FK5}},
                  T_ECIo::Type{Val{:J2000}},
                  T_ECEFf::Type{Val{:TEME}},
                  JD_UTC::Number,
                  eop_data::EOPData_IAU1980) =
    rECItoECI(T, M, T_ECIo, T_ECEFf, JD_UTC)

@inline rECItoECI(T::Type,
                  M::Type{Val{:FK5}},
                  T_ECIo::Type{Val{:TEME}},
                  T_ECEFf::Type{Val{:J2000}},
                  JD_UTC::Number,
                  eop_data::EOPData_IAU1980) =
    rECItoECI(T, M, T_ECIo, T_ECEFf, JD_UTC)

function rECItoECI(T::Type,
                   ::Type{Val{:FK5}},
                   ::Type{Val{:J2000}},
                   ::Type{Val{:TEME}},
                   JD_UTC::Number)

    # Get the time in TT.
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Return the rotation.
    rGCRFtoTEME(T, JD_TT, 0, 0)
end

@inline rECItoECI(T::Type,
                  M::Type{Val{:FK5}},
                  T_ECIo::Type{Val{:TEME}},
                  T_ECIf::Type{Val{:J2000}},
                  JD_UTC::Number) =
    inv_rotation(rECItoECI(T, M, T_ECIf, T_ECIo, JD_UTC))

#                          Between MOD, TOD, and TEME
# ==============================================================================

function rECItoECI(T::Type,
                   M::Type{Val{:FK5}},
                   T_ECIo::T_ECIs_of_date,
                   JD_UTCo::Number,
                   T_ECIf::T_ECIs_of_date,
                   JD_UTCf::Number,
                   eop_data::EOPData_IAU1980)

    # In this case, we convert origin to GCRF and then convert back to the
    # destination. This is necessary because the user may want to change the
    # epoch.
    r_GCRF_ECIo = rECItoECI(T, T_ECIo,     Val{:GCRF}, JD_UTCo, eop_data)
    r_ECIf_GCRF = rECItoECI(T, Val{:GCRF}, T_ECIf,     JD_UTCf, eop_data)

    # Return the full rotation.
    compose_rotation(r_GCRF_ECIo, r_ECIf_GCRF)
end
