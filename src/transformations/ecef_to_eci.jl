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
#   Rotations from an Earth-Fixed, Earth-Centered (ECEF) reference frame to an
#   Earth-Fixed Inertial (ECI) reference frame.
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

export rECEFtoECI

"""
    function rECEFtoECI([T,] [M,] ECEF, ECI, JD_UTC::Number [, eop_data])

Compute the rotation from an Earth-Centered, Earth-Fixed (`ECEF`) reference
frame to an Earth-Centered Inertial (`ECI`) reference frame at the Julian Day
[UTC] `JD_UTC`. The rotation description that will be used is given by `T`,
which can be `DCM` or `Quaternion`. The model used to compute the rotation is
specified by `M`. Currently, only IAU-76/FK5 is supported (`M = FK5()`). The
ECEF frame is selected by the input `ECEF` and the `ECI` frame is selected by
the input `ECI`. The possible values are listed below.

# Rotation description

The rotations that aligns the ECEF with ECI can be described by Direction Cosine
Matrices or Quaternions. This is selected by the parameter `T`. The possible
values are:

* `DCM`: The rotation will be described by a Direction Cosine Matrix.
* `Quaternion`: The rotation will be described by a Quaternion.

If no value is specified, then it falls back to `DCM`.

# Conversion model

The model that will be used to compute the rotation is given by `M`. The
possible values are:

* `FK5()`: Use the IAU-76/FK5 model.

If no value is specified, then it falls back to `FK5()`.

# ECEF Frame

The ECEF frame is selected by the parameter `ECEF`. The possible values are:

* `ITRF()`: ECEF will be selected as the International Terrestrial Reference
            Frame (ITRF).
* `PEF()`: ECEF will be selected as the Pseudo-Earth Fixed (PEF) reference
           frame.

# ECI Frame

The ECI frame is selected by the parameter `ECI`. The possible values are:

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

|   Model    |  ECEF  |   ECI   |    EOP Data   |
|:-----------|:-------|:--------|:--------------|
| IAU-76/FK5 | `ITRF` | `GCRF`  | EOP IAU1980   |
| IAU-76/FK5 | `ITRF` | `J2000` | EOP IAU1980   |
| IAU-76/FK5 | `ITRF` | `MOD`   | EOP IAU1980   |
| IAU-76/FK5 | `ITRF` | `TOD`   | EOP IAU1980   |
| IAU-76/FK5 | `PEF`  | `GCRF`  | EOP IAU1980   |
| IAU-76/FK5 | `PEF`  | `J2000` | Not required* |
| IAU-76/FK5 | `PEF`  | `MOD`   | EOP IAU1980   |
| IAU-76/FK5 | `PEF`  | `TOD`   | EOP IAU1980   |
| IAU-76/FK5 | `PEF`  | `TEME`  | Not required* |

`*`: In this case, the Julian Time UTC will be assumed equal to Julian Time UT1
to compute the Greenwich Mean Sidereal Time. This is an approximation, but
should be sufficiently accurate for some applications. Notice that, if EOP Data
is provided, the Julian Day UT1 will be accurately computed.

# Args

* `T`: (OPTIONAL) Type of the rotation representation (**Default** = `DCM`).
* `M`: (OPTIONAL) Model used to compute the rotation (**Default** = `FK5()`).
* `ECEF`: ECEF frame.
* `ECI`: ECI frame.
* `JD_UTC`: Julian day [UTC].
* `eop_data`: EOP Data.

# Returns

The rotation description represented by `T` that rotates the ECEF reference
frame into alignment with the ECI reference frame.

# Examples

```julia-repl
julia> eop_IAU1980 = get_iers_eop(:IAU1980)
julia> rECEFtoECI(DCM, FK5(), ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267      0.78518     -0.00132979
 -0.78518      -0.619267     3.33483e-5
 -0.000797314   0.00106478   0.999999

julia> rECEFtoECI(FK5(), ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267      0.78518     -0.00132979
 -0.78518      -0.619267     3.33483e-5
 -0.000797314   0.00106478   0.999999

julia> rECEFtoECI(DCM, ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267      0.78518     -0.00132979
 -0.78518      -0.619267     3.33483e-5
 -0.000797314   0.00106478   0.999999

julia> rECEFtoECI(ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267      0.78518     -0.00132979
 -0.78518      -0.619267     3.33483e-5
 -0.000797314   0.00106478   0.999999

julia> rECEFtoECI(PEF(), J2000(), DatetoJD(1986, 06, 19, 21, 35, 0))
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619271      0.785176    -0.00133066
 -0.785177     -0.619272     3.45854e-5
 -0.000796885   0.00106622   0.999999

julia> rECEFtoECI(PEF(), J2000(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 -0.619267      0.78518     -0.00133066
 -0.78518      -0.619267     3.45854e-5
 -0.000796879   0.00106623   0.999999

julia> rECEFtoECI(Quaternion, ITRF(), GCRF(), DatetoJD(1986, 06, 19, 21, 35, 0), eop_IAU1980)
Quaternion{Float64}:
  + 0.43630989232629747 - 0.0005909971869613186.i + 0.00030510471843995434.j + 0.8997962188693898.k
```
"""
@inline rECEFtoECI(T_ECEF::T_ECEFs,
                   T_ECI::T_ECIs,
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980) =
    rECEFtoECI(DCM, Val{:FK5}, T_ECEF, T_ECI, JD_UTC, eop_data)

@inline rECEFtoECI(T::Type,
                   T_ECEF::T_ECEFs,
                   T_ECI::T_ECIs,
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980) =
    rECEFtoECI(T, Val{:FK5}, T_ECEF, T_ECI, JD_UTC, eop_data)

# Specializations for those cases that EOP Data is not needed.
@inline rECEFtoECI(T_ECEF::Type{Val{:PEF}},
                   T_ECI::Union{Type{Val{:J2000}}, Type{Val{:TEME}}},
                   JD_UTC::Number) =
    rECEFtoECI(DCM, Val{:FK5}, T_ECEF, T_ECI, JD_UTC)

@inline rECEFtoECI(T::Type,
                   T_ECEF::Type{Val{:PEF}},
                   T_ECI::Union{Type{Val{:J2000}}, Type{Val{:TEME}}},
                   JD_UTC::Number) =
    rECEFtoECI(T, Val{:FK5}, T_ECEF, T_ECI, JD_UTC)

################################################################################
#                                  IAU-76/FK5
################################################################################

# Specialization related to the default type of the rotation representation.
@inline rECEFtoECI(::Type{Val{:FK5}},
                   T_ECEF::T_ECEFs,
                   T_ECI::T_ECIs,
                   JD_UTC::Number,
                   eop_data::EOPData_IAU1980) =
    rECEFtoECI(DCM, Val{:FK5}, T_ECEF, T_ECI, JD_UTC, eop_data)

# Specializations for those cases that EOP Data is not needed.
@inline rECEFtoECI(::Type{Val{:FK5}},
                   T_ECEF::Type{Val{:PEF}},
                   T_ECI::Union{Type{Val{:J2000}}, Type{Val{:TEME}}},
                   JD_UTC::Number) =
    rECEFtoECI(DCM, Val{:FK5}, T_ECEF, T_ECI, JD_UTC)

#                                 ITRF => GCRF
# ==============================================================================

function rECEFtoECI(T::Type,
                    ::Type{Val{:FK5}},
                    ::Type{Val{:ITRF}},
                    ::Type{Val{:GCRF}},
                    JD_UTC::Number,
                    eop_data::EOPData_IAU1980)

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x[JD_UTC]*pi/648000
    y_p      = eop_data.y[JD_UTC]*pi/648000
    δΔϵ_1980 = eop_data.dEps[JD_UTC]*pi/648000
    δΔψ_1980 = eop_data.dPsi[JD_UTC]*pi/648000

    # Compute the rotation.
    rITRFtoGCRF_fk5(T, JD_UT1, JD_TT, x_p, y_p, δΔϵ_1980, δΔψ_1980)
end

#                                ITRF => J2000
# ==============================================================================

function rECEFtoECI(T::Type,
                    ::Type{Val{:FK5}},
                    ::Type{Val{:ITRF}},
                    ::Type{Val{:J2000}},
                    JD_UTC::Number,
                    eop_data::EOPData_IAU1980)

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x[JD_UTC]*pi/648000
    y_p      = eop_data.y[JD_UTC]*pi/648000

    # Compute the rotation.
    rITRFtoGCRF_fk5(T, JD_UT1, JD_TT, x_p, y_p, 0, 0)
end

#                                 ITRF => MOD
# ==============================================================================

function rECEFtoECI(T::Type,
                    ::Type{Val{:FK5}},
                    ::Type{Val{:ITRF}},
                    ::Type{Val{:MOD}},
                    JD_UTC::Number,
                    eop_data::EOPData_IAU1980)

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x[JD_UTC]*pi/648000
    y_p      = eop_data.y[JD_UTC]*pi/648000
    δΔϵ_1980 = eop_data.dEps[JD_UTC]*pi/648000
    δΔψ_1980 = eop_data.dPsi[JD_UTC]*pi/648000

    # Compute the rotation.
    r_PEF_ITRF = rITRFtoPEF_fk5(T, x_p, y_p)
    r_MOD_PEF  = rPEFtoMOD_fk5(T, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)

    compose_rotation(r_PEF_ITRF, r_MOD_PEF)
end

#                                 ITRF => TOD
# ==============================================================================

function rECEFtoECI(T::Type,
                    ::Type{Val{:FK5}},
                    ::Type{Val{:ITRF}},
                    ::Type{Val{:TOD}},
                    JD_UTC::Number,
                    eop_data::EOPData_IAU1980)

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x[JD_UTC]*pi/648000
    y_p      = eop_data.y[JD_UTC]*pi/648000
    δΔψ_1980 = eop_data.dPsi[JD_UTC]*pi/648000

    # Compute the rotation.
    r_PEF_ITRF = rITRFtoPEF_fk5(T, x_p, y_p)
    r_TOD_PEF  = rPEFtoTOD_fk5(T, JD_UT1, JD_TT, δΔψ_1980)

    compose_rotation(r_PEF_ITRF, r_TOD_PEF)
end


#                                 PEF => GCRF
# ==============================================================================

function rECEFtoECI(T::Type,
                    ::Type{Val{:FK5}},
                    ::Type{Val{:PEF}},
                    ::Type{Val{:GCRF}},
                    JD_UTC::Number,
                    eop_data::EOPData_IAU1980)

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x[JD_UTC]*pi/648000
    y_p      = eop_data.y[JD_UTC]*pi/648000
    δΔϵ_1980 = eop_data.dEps[JD_UTC]*pi/648000
    δΔψ_1980 = eop_data.dPsi[JD_UTC]*pi/648000

    # Compute the rotation.
    r_MOD_PEF  = rPEFtoMOD_fk5(T, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)
    r_GCRF_MOD = rMODtoGCRF_fk5(T, JD_TT)

    compose_rotation(r_MOD_PEF, r_GCRF_MOD)
end

#                                PEF => J2000
# ==============================================================================

function rECEFtoECI(T::Type,
                    ::Type{Val{:FK5}},
                    ::Type{Val{:PEF}},
                    ::Type{Val{:J2000}},
                    JD_UTC::Number,
                    eop_data::EOPData_IAU1980)

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Compute the rotation.
    r_MOD_PEF  = rPEFtoMOD_fk5(T, JD_UT1, JD_TT, 0, 0)
    r_GCRF_MOD = rMODtoGCRF_fk5(T, JD_TT)

    compose_rotation(r_MOD_PEF, r_GCRF_MOD)
end

function rECEFtoECI(T::Type,
                    ::Type{Val{:FK5}},
                    ::Type{Val{:PEF}},
                    ::Type{Val{:J2000}},
                    JD_UTC::Number)

    # Since we do not have EOP Data, assume that JD_UTC is equal to JD_UT1.
    JD_UT1 = JD_UTC
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Compute the rotation.
    r_MOD_PEF  = rPEFtoMOD_fk5(T, JD_UT1, JD_TT, 0, 0)
    r_GCRF_MOD = rMODtoGCRF_fk5(T, JD_TT)

    compose_rotation(r_MOD_PEF, r_GCRF_MOD)
end

#                                 PEF => MOD
# ==============================================================================

function rECEFtoECI(T::Type,
                    ::Type{Val{:FK5}},
                    ::Type{Val{:PEF}},
                    ::Type{Val{:MOD}},
                    JD_UTC::Number,
                    eop_data::EOPData_IAU1980)
    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x[JD_UTC]*pi/648000
    y_p      = eop_data.y[JD_UTC]*pi/648000
    δΔϵ_1980 = eop_data.dEps[JD_UTC]*pi/648000
    δΔψ_1980 = eop_data.dPsi[JD_UTC]*pi/648000

    # Compute the rotation.
    rPEFtoMOD_fk5(T, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)
end

#                                 PEF => TOD
# ==============================================================================

function rECEFtoECI(T::Type,
                    ::Type{Val{:FK5}},
                    ::Type{Val{:PEF}},
                    ::Type{Val{:TOD}},
                    JD_UTC::Number,
                    eop_data::EOPData_IAU1980)

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    x_p      = eop_data.x[JD_UTC]*pi/648000
    y_p      = eop_data.y[JD_UTC]*pi/648000
    δΔψ_1980 = eop_data.dPsi[JD_UTC]*pi/648000

    # Compute the rotation.
    rPEFtoTOD_fk5(T, JD_UT1, JD_TT, δΔψ_1980)
end


#                                 PEF => TEME
# ==============================================================================

function rECEFtoECI(T::Type,
                    ::Type{Val{:FK5}},
                    ::Type{Val{:PEF}},
                    ::Type{Val{:TEME}},
                    JD_UTC::Number,
                    eop_data::EOPData_IAU1980)
    # Get the time in UT1.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)

    # Compute the rotation.
    rPEFtoTEME(T, JD_UT1)
end

function rECEFtoECI(T::Type,
                    ::Type{Val{:FK5}},
                    ::Type{Val{:PEF}},
                    ::Type{Val{:TEME}},
                    JD_UTC::Number)
    # Since we do not have EOP Data, assume that JD_UTC is equal to JD_UT1.
    JD_UT1 = JD_UTC

    # Compute the rotation.
    rPEFtoTEME(T, JD_UT1)
end
