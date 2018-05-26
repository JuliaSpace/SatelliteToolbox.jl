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
two *of date* frame, e.g. TOD => MOD, then the second signature must be used in
which `JD_UTCo` is the epoch of the origin frame and `JD_UTCf` is the epoch of
the destination frame.

The rotation description that will be used is given by `T`, which can be `DCM`
or `Quaternion`. The model used to compute the rotation is specified by `M`.
Currently, only IAU-76/FK5 is supported (`M = FK5()`). The input ECI frame is
selected by the input `ECIo` and the destination ECI frame is selected by the
input `ECIf`.

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

|   Model    |   ECIo  |   ECIf  |    EOP Data   |
|:-----------|:--------|:--------|:--------------|
| IAU-76/FK5 | `GCRF`  | `J2000` | EOP IAU1980   |
| IAU-76/FK5 | `GCRF`  | `MOD`   | EOP IAU1980   |
| IAU-76/FK5 | `GCRF`  | `TOD`   | EOP IAU1980   |
| IAU-76/FK5 | `GCRF`  | `TEME`  | EOP IAU1980   |
| IAU-76/FK5 | `J2000` | `GCRF`  | EOP IAU1980   |
| IAU-76/FK5 | `J2000` | `MOD`   | EOP IAU1980   |
| IAU-76/FK5 | `J2000` | `TOD`   | EOP IAU1980   |
| IAU-76/FK5 | `J2000` | `TEME`  | EOP IAU1980   |
| IAU-76/FK5 | `TOD`   | `GCRF`  | EOP IAU1980   |
| IAU-76/FK5 | `TOD`   | `J2000` | EOP IAU1980   |
| IAU-76/FK5 | `TOD`   | `MOD`   | EOP IAU1980   |
| IAU-76/FK5 | `TOD`   | `TEME`  | EOP IAU1980   |
| IAU-76/FK5 | `MOD`   | `GCRF`  | EOP IAU1980   |
| IAU-76/FK5 | `MOD`   | `J2000` | EOP IAU1980   |
| IAU-76/FK5 | `MOD`   | `TOD`   | EOP IAU1980   |
| IAU-76/FK5 | `MOD`   | `TEME`  | EOP IAU1980   |
| IAU-76/FK5 | `TEME`  | `GCRF`  | EOP IAU1980   |
| IAU-76/FK5 | `TEME`  | `J2000` | EOP IAU1980   |
| IAU-76/FK5 | `TEME`  | `MOD`   | EOP IAU1980   |
| IAU-76/FK5 | `TEME`  | `TOD`   | EOP IAU1980   |

`*`: In this case, the Julian Time UTC will be assumed equal to Julian Time UT1
to compute the Greenwich Mean Sidereal Time. This is an approximation, but
should be sufficiently accurate for some applications. Notice that, if EOP Data
is provided, the Julian Day UT1 will be accurately computed.

## MOD and TOD

In this function, MOD and TOD frames are always defined with IERS EOP
corrections. Hence, if one wants to obtain the MOD and TOD frames according to
the original IAU-76/FK5 theory, it is necessary to use the low-level functions
in file `./src/transformations/fk5/fk5.jl`.

# Args

* `T`: (OPTIONAL) Type of the rotation representation (**Default** = `DCM`).
* `M`: (OPTIONAL) Model used to compute the rotation (**Default** = `FK5()`).
* `ECIo`: Input ECEF frame.
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
