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
#   Rotations from an Earth-Fixed Inertial (ECI) reference frame to an
#   Earth-Fixed, Earth-Centered (ECEF) reference frame.
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
#   Add support to TEME.
#
# 2018-05-09: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export rECItoECEF

"""
    function rECItoECEF([T,] [M,] ECI, ECEF, JD_UTC::Number [, eop_data])

Compute the rotation from an Earth-Centered Inertial (`ECI`) reference frame to
an Earth-Centered, Earth-Fixed (`ECEF`) reference frame at the Julian Day [UTC]
`JD_UTC`. The rotation description that will be used is given by `T`, which can
be `DCM` or `Quaternion`. The model used to compute the rotation is specified by
`M`. Currently, only IAU-76/FK5 is supported (`M = FK5()`). The ECI frame is
selected by the input `ECI` and the `ECEF` frame is selected by the input
`ECEF`. The possible values are listed below.

# Rotation description

The rotations that aligns the ECI with ECEF can be described by Direction Cosine
Matrices or Quaternion. This is selected by the parameter `T`. The possible
values are:

* `DCM`: The rotation will be described by a Direction Cosine Matrix.
* `Quaternion`: The rotation will be described by a Quaternion.

If no value is specified, then it falls back to `DCM`.

# Conversion model

The model that will be used to compute the rotation is given by `M`. The
possible values are:

* `FK5()`: Use the IAU-76/FK5 model.

If no value is specified, then it falls back to `FK5()`.

# ECI Frame

The ECI frame is selected by the parameter `ECI`. The possible values are:

* `TEME()`: ECI will be selected as the True Equator Mean Equinox (TEME)
            reference frame.
* `TOD()`: ECI will be selected as the True of Date (TOD).
* `MOD()`: ECI will be selected as the Mean of Date (MOD).
* `J2000()`: ECI will be selected as the J2000 reference frame.
* `GCRF()`: ECI will be selected as the Geocentric Celestial Reference Frame
            (GCRF).

# ECEF Frame

The ECEF frame is selected by the parameter `ECEF`. The possible values are:

* `ITRF()`: ECEF will be selected as the International Terrestrial Reference
            Frame (ITRF).
* `PEF()`: ECEF will be selected as the Pseudo-Earth Fixed (PEF) reference
           frame.

# EOP Data

The conversion between the frames depends on EOP Data (see `get_iers_eop` and
`read_iers_eop`). If IAU-76/FK5 model is used, then the type of `eop_data` must
be `EOPData_IAU1980`. The following table shows the requirements for EOP data
given the selected frames.

|   Model    |   ECI   |  ECEF  |    EOP Data   |
|:-----------|:--------|:-------|:--------------|
| IAU-76/FK5 | `GCRF`  | `ITRF` | EOP IAU1980   |
| IAU-76/FK5 | `J2000` | `ITRF` | EOP IAU1980   |
| IAU-76/FK5 | `MOD`   | `ITRF` | EOP IAU1980   |
| IAU-76/FK5 | `TOD`   | `ITRF` | EOP IAU1980   |
| IAU-76/FK5 | `GCRF`  | `PEF`  | EOP IAU1980   |
| IAU-76/FK5 | `J2000` | `PEF`  | Not required* |
| IAU-76/FK5 | `MOD`   | `PEF`  | EOP IAU1980   |
| IAU-76/FK5 | `TOD`   | `PEF`  | EOP IAU1980   |
| IAU-76/FK5 | `TEME`  | `PEF`  | Not required* |

`*`: In this case, the Julian Time UTC will be assumed equal to Julian Time UT1
to compute the Greenwich Mean Sidereal Time. This is an approximation, but
should be sufficiently accurate for some applications. Notice that, if EOP Data
is provided, the Julian Day UT1 will be accurately computed.

# Args

* `T`: (OPTIONAL) Type of the rotation representation (**Default** = `DCM`).
* `M`: (OPTIONAL) Model used to compute the rotation (**Default** = `FK5()`).
* `ECI`: ECI frame.
* `ECEF`: ECEF frame.
* `JD_UTC`: Julian day [UTC].
* `eop_data`: EOP Data.

# Returns

The rotation description represented by `T` that rotates the ECI reference frame
into alignment with the ECEF reference frame.

# Examples

```julia-repl
julia> eop_IAU1980 = get_iers_eop(:IAU1980)
julia> rECItoECEF(DCM, FK5(), GCRF(), ITRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267    -0.78518     -0.000797314
  0.78518     -0.619267     0.00106478
 -0.00132979   3.33483e-5   0.999999

julia> rECItoECEF(FK5(), GCRF(), ITRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267    -0.78518     -0.000797314
  0.78518     -0.619267     0.00106478
 -0.00132979   3.33483e-5   0.999999

julia> rECItoECEF(DCM, GCRF(), ITRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267    -0.78518     -0.000797314
  0.78518     -0.619267     0.00106478
 -0.00132979   3.33483e-5   0.999999

julia> rECItoECEF(GCRF(), ITRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267    -0.78518     -0.000797314
  0.78518     -0.619267     0.00106478
 -0.00132979   3.33483e-5   0.999999

julia> rECItoECEF(J2000(), PEF(), DatetoJD(1986, 06, 19, 21, 35, 0))
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619271    -0.785177    -0.000796885
  0.785176    -0.619272     0.00106622
 -0.00133066   3.45854e-5   0.999999

julia> rECItoECEF(J2000(), PEF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267    -0.78518     -0.000796879
  0.78518     -0.619267     0.00106623
 -0.00133066   3.45854e-5   0.999999

julia> rECItoECEF(Quaternion, ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
Quaternion{Float64}:
  + 0.43630989232629747 + 0.0005909971869613186.i - 0.00030510471843995434.j - 0.8997962188693898.k
```
"""
@inline rECItoECEF(T_ECI::T_ECIs,
                   T_ECEF::T_ECEFs,
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980) =
    rECItoECEF(DCM, Val{:FK5}, T_ECI, T_ECEF, JD_UTC, eop_data)

@inline rECItoECEF(T::Union{Type{DCM},Type{Quaternion}},
                   T_ECI::T_ECIs,
                   T_ECEF::T_ECEFs,
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980) =
    rECItoECEF(T, Val{:FK5}, T_ECI, T_ECEF, JD_UTC, eop_data)

@inline rECItoECEF(M::Type{Val{:FK5}},
                   T_ECI::T_ECIs,
                   T_ECEF::T_ECEFs,
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980) =
    rECItoECEF(DCM, M, T_ECI, T_ECEF, JD_UTC, eop_data)

@inline rECItoECEF(T::Union{Type{DCM},Type{Quaternion}},
                   M::Type{Val{:FK5}},
                   T_ECI::T_ECIs,
                   T_ECEF::T_ECEFs,
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980) =
    inv_rotation(rECEFtoECI(T, M, T_ECEF, T_ECI, JD_UTC, eop_data))

# Specializations for those cases that EOP Data is not needed.
@inline rECItoECEF(T_ECI::Union{Type{Val{:J2000}},Type{Val{:TEME}}},
                   T_ECEF::Type{Val{:PEF}},
                   JD_UTC::Number) =
    rECItoECEF(DCM, Val{:FK5}, T_ECI, T_ECEF, JD_UTC)

@inline rECItoECEF(M::Type{Val{:FK5}},
                   T_ECI::Union{Type{Val{:J2000}},Type{Val{:TEME}}},
                   T_ECEF::Type{Val{:PEF}},
                   JD_UTC::Number) =
    rECItoECEF(DCM, M, T_ECEF, T_ECI, JD_UTC)

@inline rECItoECEF(T::Union{Type{DCM},Type{Quaternion}},
                   T_ECI::Union{Type{Val{:J2000}},Type{Val{:TEME}}},
                   T_ECEF::Type{Val{:PEF}},
                   JD_UTC::Number) =
    rECItoECEF(T, Val{:FK5}, T_ECEF, T_ECI, JD_UTC)

@inline rECItoECEF(T::Union{Type{DCM},Type{Quaternion}},
                   M::Type{Val{:FK5}},
                   T_ECI::Union{Type{Val{:J2000}},Type{Val{:TEME}}},
                   T_ECEF::Type{Val{:PEF}},
                   JD_UTC::Number) =
    inv_rotation(rECEFtoECI(T, M, T_ECEF, T_ECI, JD_UTC))
