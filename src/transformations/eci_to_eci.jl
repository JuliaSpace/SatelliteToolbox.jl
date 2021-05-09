# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Rotations from an Earth-Centered Inertial (ECI) reference frame to another
#   ECI reference frame.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export r_eci_to_eci

"""
    r_ecef_to_eci([T,] ECIo, ECIf, JD_UTC::Number [, eop_data])
    r_ecef_to_eci([T,] ECIo, JD_UTCo::Number, ECIf, JD_UTCf::Number [, eop_data])

Compute the rotation from an Earth-Centered Inertial (`ECI`) reference frame to
another ECI reference frame. If the origin and destination frame contain only
one *of date* frame, then the first signature is used and `JD_UTC` is the epoch
of this frame. On the other hand, if the origin and destination frame contain
two *of date* frame[^1], e.g. TOD => MOD, then the second signature must be used
in which `JD_UTCo` is the epoch of the origin frame and `JD_UTCf` is the epoch
of the destination frame.

The rotation description that will be used is given by `T`, which can be `DCM`
or `Quaternion`. The origin ECI frame is selected by the input `ECIo` and the
destination ECI frame is selected by the input `ECIf`. The model used to compute
the rotation is specified by the selection of the origin and destination frames.
Currently, there are two models supported: IAU-76/FK5 and IAU-2006 with 2010
conventions (CIO and equinox approaches).

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

The model that will be used to compute the rotation is automatically inferred
given the selection of the origin and destination frames. **Notice that mixing
IAU-76/FK5 and IAU-2006/2010 frames is not supported yet.**

# ECI Frame

The supported ECI frames for both origin `ECIo` and destination `ECIf` are:

* `TEME()`: ECI will be selected as the True Equator Mean Equinox (TEME)
            reference frame.
* `TOD()`: ECI will be selected as the True of Date (TOD).
* `MOD()`: ECI will be selected as the Mean of Date (MOD).
* `J2000()`: ECI will be selected as the J2000 reference frame.
* `GCRF()`: ECI will be selected as the Geocentric Celestial Reference Frame
            (GCRF).
* `CIRS()`: ECEF will be selected as the Celestial Intermediate Reference System
            (CIRS).
* `ERS()`: ECI will be selected as the Earth Reference System (ERS).
* `MOD06()`: ECI will be selected as the Mean of Date (MOD) according to the
             definition in IAU-2006/2010 theory.
* `MJ2000()`: ECI will be selected as the J2000 mean equatorial frame (MJ2000).

!!! note

    The frames `MOD()` and `MOD06()` are virtually the same. However, we
    selected different names to make clear which theory are being used since
    mixing transformation between frames from IAU-76/FK5 and IAU-2006/2010 must
    be performed with caution.

# EOP Data

The conversion between the frames depends on EOP Data (see `get_iers_eop` and
`read_iers_eop`). If IAU-76/FK5 model is used, then the type of `eop_data` must
be `EOPData_IAU1980`. Otherwise, if IAU-2006/2010 model is used, then the type
of `eop_data` must be `EOPData_IAU2000A`. The following table shows the
requirements for EOP data given the selected frames.

|   Model                     |   ECIo   |   ECIf   |    EOP Data   | Function Signature |
|:----------------------------|:---------|:---------|:--------------|:-------------------|
| IAU-76/FK5                  | `GCRF`   | `J2000`  | EOP IAU1980   | First              |
| IAU-76/FK5                  | `GCRF`   | `MOD`    | EOP IAU1980   | First              |
| IAU-76/FK5                  | `GCRF`   | `TOD`    | EOP IAU1980   | First              |
| IAU-76/FK5                  | `GCRF`   | `TEME`   | EOP IAU1980   | First              |
| IAU-76/FK5                  | `J2000`  | `GCRF`   | EOP IAU1980   | First              |
| IAU-76/FK5                  | `J2000`  | `MOD`    | Not required  | First              |
| IAU-76/FK5                  | `J2000`  | `TOD`    | Not required  | First              |
| IAU-76/FK5                  | `J2000`  | `TEME`   | Not required  | First              |
| IAU-76/FK5                  | `MOD`    | `GCRF`   | EOP IAU1980   | First              |
| IAU-76/FK5                  | `MOD`    | `J2000`  | Not required  | First              |
| IAU-76/FK5                  | `MOD`    | `TOD`    | Not required  | Second             |
| IAU-76/FK5                  | `MOD`    | `TEME`   | Not required  | Second             |
| IAU-76/FK5                  | `TOD`    | `GCRF`   | EOP IAU1980   | First              |
| IAU-76/FK5                  | `TOD`    | `J2000`  | Not required  | First              |
| IAU-76/FK5                  | `TOD`    | `MOD`    | Not required  | Second             |
| IAU-76/FK5                  | `TOD`    | `TEME`   | Not required  | Second             |
| IAU-76/FK5                  | `TEME`   | `GCRF`   | EOP IAU1980   | First              |
| IAU-76/FK5                  | `TEME`   | `J2000`  | Not required  | First              |
| IAU-76/FK5                  | `TEME`   | `MOD`    | Not required  | Second             |
| IAU-76/FK5                  | `TEME`   | `TOD`    | Not required  | Second             |
| IAU-2006/2010 CIO-based     | `GCRF`   | `CIRS`   | Not required¹ | First              |
| IAU-2006/2010 CIO-based     | `CIRS`   | `CIRS`   | Not required¹ | Second             |
| IAU-2006/2010 Equinox-based | `GCRF`   | `MJ2000` | Not required  | First²             |
| IAU-2006/2010 Equinox-based | `GCRF`   | `MOD06`  | Not required  | First              |
| IAU-2006/2010 Equinox-based | `GCRF`   | `ERS`    | Not required³ | First              |
| IAU-2006/2010 Equinox-based | `MJ2000` | `GCRF`   | Not required  | First²             |
| IAU-2006/2010 Equinox-based | `MJ2000` | `MOD06`  | Not required  | First              |
| IAU-2006/2010 Equinox-based | `MJ2000` | `ERS`    | Not required³ | First              |
| IAU-2006/2010 Equinox-based | `MOD06`  | `GCRF`   | Not required  | First              |
| IAU-2006/2010 Equinox-based | `MOD06`  | `MJ2000` | Not required  | First              |
| IAU-2006/2010 Equinox-based | `MOD06`  | `ERS`    | Not required³ | First              |
| IAU-2006/2010 Equinox-based | `ERS`    | `GCRF`   | Not required³ | First              |
| IAU-2006/2010 Equinox-based | `ERS`    | `MJ2000` | Not required³ | First              |
| IAU-2006/2010 Equinox-based | `ERS`    | `MOD06`  | Not required³ | First              |

`¹`: In this case, the terms that account for the free-core nutation and time
dependent effects of the Celestial Intermediate Pole (CIP) position with respect
to the GCRF will not be available, reducing the precision.

`²`: The transformation between GCRF and MJ2000 is a constant rotation matrix
called bias. Hence, the date does not modify it. However, this signature was
kept to avoid complications in the API.

`³`: In this case, the terms that corrects the nutation in obliquity and in
longitude due to the free core nutation will not be available, reducing the
precision.

## MOD and TOD

In this function, if EOP corrections are not provided, then MOD and TOD frames
will be computed considering the original IAU-76/FK5 theory. Otherwise, the
corrected frame will be used.

# Returns

The rotation description represented by `T` that rotates the origin ECI
reference frame into alignment with the destination ECI reference frame.

# Examples

```julia-repl
julia> eop_IAU1980 = get_iers_eop(:IAU1980);

julia> r_eci_to_eci(DCM, GCRF(), J2000(), DatetoJD(1986, 6, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
  1.0          -2.45469e-12   4.56602e-10
  2.45466e-12   1.0          -1.84455e-9
 -4.56602e-10   1.84455e-9    1.0

julia> r_eci_to_eci(Quaternion, TEME(), GCRF(), DatetoJD(1986, 6, 19, 21, 35, 0), eop_IAU1980)
Quaternion{Float64}:
  + 0.9999986335698654 + 1.8300414020900853e-5.i + 0.0006653038276169474.j - 0.0015132396749411375.k

julia> r_eci_to_eci(TOD(), DatetoJD(1986,6,19,21,35,0), TOD(), DatetoJD(1987,5,19,3,0,0), eop_IAU1980)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 1.0          -0.000224087  -9.73784e-5
 0.000224086   1.0          -5.79859e-6
 9.73797e-5    5.77677e-6    1.0

julia> r_eci_to_eci(Quaternion, TOD(), JD_J2000, MOD(), JD_J2000, eop_IAU1980)
Quaternion{Float64}:
  + 0.9999999993282687 - 1.400220690336851e-5.i + 1.3473593746216003e-5.j - 3.107834312843103e-5.k

julia> r_eci_to_eci(J2000(), TEME(), DatetoJD(1986,6,19,21,35,0))
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
  0.999995    0.0030265    0.00133055
 -0.00302645  0.999995    -3.86125e-5
 -0.00133066  3.45854e-5   0.999999

julia> eop_IAU2000A = get_iers_eop(:IAU2000A);

julia> r_eci_to_eci(CIRS(), GCRF(), DatetoJD(1986,6,19,21,35,0), eop_IAU2000A)
3×3 StaticArrays.SArray{Tuple{3,3},Float64,2,9}:
 0.999999     3.88379e-8  -0.00133066
 7.18735e-9   1.0          3.45882e-5
 0.00133066  -3.45882e-5   0.999999

julia> r_eci_to_eci(Quaternion, CIRS(), GCRF(), DatetoJD(1986,6,19,21,35,0), eop_IAU2000A)
Quaternion{Float64}:
  + 0.9999997785177528 + 1.7294102099105917e-5.i + 0.0006653310148723835.j + 7.912627369563795e-9.k
```
"""
@inline function r_eci_to_eci(
    T_ECIo::T_ECIs,
    T_ECIf::T_ECIs,
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    return r_eci_to_eci(DCM, T_ECIo, T_ECIf, JD_UTC, eop_data)
end

@inline function r_eci_to_eci(
    T_ECIo::T_ECIs_IAU_2006,
    T_ECIf::T_ECIs_IAU_2006,
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(DCM, T_ECIo, T_ECIf, JD_UTC, eop_data)
end

# Specializations for those cases in which we have two *of dates* frames.
@inline function r_eci_to_eci(
    T_ECIo::T_ECIs_of_date,
    JD_UTCo::Number,
    T_ECIf::T_ECIs_of_date,
    JD_UTCf::Number,
    eop_data::EOPData_IAU1980
)
    return r_eci_to_eci(DCM, T_ECIo, JD_UTCo, T_ECIf, JD_UTCf, eop_data)
end

@inline function r_eci_to_eci(
    T_ECIo::Val{:CIRS},
    JD_UTCo::Number,
    T_ECIf::Val{:CIRS},
    JD_UTCf::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(DCM, T_ECIo, JD_UTCo, T_ECIf, JD_UTCf, eop_data)
end

@inline function r_eci_to_eci(
    T_ECIo::T_ECIs_IAU_2006_Equinox_of_date,
    JD_UTCo::Number,
    T_ECIf::T_ECIs_IAU_2006_Equinox_of_date,
    JD_UTCf::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(DCM, T_ECIo, JD_UTCo, T_ECIf, JD_UTCf, eop_data)
end

# Specializations for those cases that EOP Data is not needed.
@inline function r_eci_to_eci(
    T_ECIo::Val{:J2000},
    T_ECIf::Union{Val{:MOD}, Val{:TOD}, Val{:TEME}},
    JD_UTC::Number
)
    return r_eci_to_eci(DCM, T_ECIo, T_ECIf, JD_UTC)
end

@inline function r_eci_to_eci(
    T_ECIo::Union{Val{:MOD}, Val{:TOD}, Val{:TEME}},
    T_ECIf::Val{:J2000},
    JD_UTC::Number
)
    return r_eci_to_eci(DCM, T_ECIo, T_ECIf, JD_UTC)
end

@inline function r_eci_to_eci(
    T_ECIo::T_ECIs_of_date,
    JD_UTCo::Number,
    T_ECIf::T_ECIs_of_date,
    JD_UTCf::Number
)
    return r_eci_to_eci(DCM, T_ECIo, JD_UTCo, T_ECIf, JD_UTCf)
end

@inline function r_eci_to_eci(
    T_ECIo::T_ECIs_IAU_2006,
    T_ECIf::T_ECIs_IAU_2006,
    JD_UTC::Number
)
    return r_eci_to_eci(DCM, T_ECIo, T_ECIf, JD_UTC)
end

 @inline function r_eci_to_eci(
    T_ECIo::Val{:CIRS},
    JD_UTCo::Number,
    T_ECIf::Val{:CIRS},
    JD_UTCf::Number
 )
    return r_eci_to_eci(DCM, T_ECIo, JD_UTCo, T_ECIf, JD_UTCf)
end

@inline function r_eci_to_eci(
    T_ECIo::T_ECIs_IAU_2006_Equinox_of_date,
    T_ECIf::T_ECIs_IAU_2006_Equinox_of_date,
    JD_UTC::Number
)
    return r_eci_to_eci(DCM, T_ECIo, T_ECIf, JD_UTC)
end

################################################################################
#                                  IAU-76/FK5
################################################################################

#                                GCRF <=> J2000
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:J2000},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    δΔϵ_1980 = eop_data.dEps(JD_UTC)*arcsec2rad
    δΔψ_1980 = eop_data.dPsi(JD_UTC)*arcsec2rad

    # In this case, we need to convert GCRF back to PEF and then convert to
    # J2000, which is the same conversion from PEF to GCRF **without** the EOP
    # data.
    #
    # TODO: Can I simplify the rotation from TOD with corrections to TOD without
    # corrections?
    r_MOD_GCRF  = r_gcrf_to_mod_fk5(T, JD_TT)
    r_PEF_MOD   = r_mod_to_pef_fk5(T, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)
    r_MOD_PEF   = r_pef_to_mod_fk5(T, JD_UT1, JD_TT, 0, 0)
    r_J2000_MOD = inv_rotation(r_MOD_GCRF)

    # Compose the full rotation.
    return compose_rotation(r_MOD_GCRF, r_PEF_MOD, r_MOD_PEF, r_J2000_MOD)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:J2000},
    T_ECIf::Val{:GCRF},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, JD_UTC, eop_data))
end

#                                 GCRF <=> MOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:MOD},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    # Get the time in TT.
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Return the rotation.
    return r_gcrf_to_mod_fk5(T, JD_TT)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:MOD},
    T_ECIf::Val{:GCRF},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, JD_UTC, eop_data))
end

#                                 GCRF <=> TOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:TOD},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in TT.
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    δΔϵ_1980 = eop_data.dEps(JD_UTC)*arcsec2rad
    δΔψ_1980 = eop_data.dPsi(JD_UTC)*arcsec2rad

    # Return the rotation.
    r_MOD_GCRF = r_gcrf_to_mod_fk5(T, JD_TT)
    r_TOD_MOD  = r_mod_to_tod_fk5(T, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Compose the full rotation.
    return compose_rotation(r_MOD_GCRF, r_TOD_MOD)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:TOD},
    T_ECIf::Val{:GCRF},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, JD_UTC, eop_data))
end

#                                GCRF <=> TEME
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:TEME},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in TT.
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    δΔϵ_1980 = eop_data.dEps(JD_UTC)*arcsec2rad
    δΔψ_1980 = eop_data.dPsi(JD_UTC)*arcsec2rad

    # Return the rotation.
    r_MOD_GCRF = r_gcrf_to_mod_fk5(T, JD_TT)
    r_TEME_MOD = rMODtoTEME(T, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Compose the full rotation.
    return compose_rotation(r_MOD_GCRF, r_TEME_MOD)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:TEME},
    T_ECIf::Val{:GCRF},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, JD_UTC, eop_data))
end

#                                J2000 <=> MOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:J2000},
    ::Val{:MOD},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    δΔϵ_1980 = eop_data.dEps(JD_UTC)*arcsec2rad
    δΔψ_1980 = eop_data.dPsi(JD_UTC)*arcsec2rad

    # In this case, we need to convert J2000 back to PEF and then convert to
    # MOD. This is necessary because we need to apply EOP corrections to convert
    # to MOD and just a `r_gcrf_to_mod_fk5` would yield a frame according to the
    # original IAU-76/FK5 theory.
    #
    # TODO: Can I simplify this rotation?
    r_MOD_J2000 = r_gcrf_to_mod_fk5(T, JD_TT)
    r_PEF_MOD   = r_mod_to_pef_fk5(T, JD_UT1, JD_TT, 0, 0)
    r_MOD_PEF   = r_pef_to_mod_fk5(T, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Compose the full rotation.
    return compose_rotation(r_MOD_J2000, r_PEF_MOD, r_MOD_PEF)
end

function r_eci_to_eci(T::T_ROT, ::Val{:J2000}, ::Val{:MOD}, JD_UTC::Number)
    # Get the time in TT.
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Compute and return the rotation.
    r_MOD_J2000 = r_gcrf_to_mod_fk5(T, JD_TT)

    return r_MOD_J2000
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:MOD},
    T_ECIf::Val{:J2000},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, JD_UTC, eop_data))
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:MOD},
    T_ECIf::Val{:J2000},
    JD_UTC::Number
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, JD_UTC))
end

#                                J2000 <=> TOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:J2000},
    ::Val{:TOD},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    arcsec2rad = π/648000

    # Get the time in UT1 and TT.
    JD_UT1 = JD_UTCtoUT1(JD_UTC, eop_data)
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    #
    # TODO: The difference is small, but should it be `JD_TT` or `JD_UTC`?
    δΔϵ_1980 = eop_data.dEps(JD_UTC)*arcsec2rad
    δΔψ_1980 = eop_data.dPsi(JD_UTC)*arcsec2rad

    # In this case, we need to convert J2000 back to PEF and then convert to
    # TOD. This is necessary because we need to apply EOP corrections to convert
    # to TOD.
    #
    # TODO: Can I simplify this rotation?
    r_MOD_J2000 = r_gcrf_to_mod_fk5(T, JD_TT)
    r_PEF_MOD   = r_mod_to_pef_fk5(T, JD_UT1, JD_TT, 0, 0)
    r_MOD_PEF   = r_pef_to_mod_fk5(T, JD_UT1, JD_TT, δΔϵ_1980, δΔψ_1980)
    r_TOD_MOD   = r_mod_to_tod_fk5(T, JD_TT, δΔϵ_1980, δΔψ_1980)

    # Compose the full rotation.
    return compose_rotation(r_MOD_J2000, r_PEF_MOD, r_MOD_PEF, r_TOD_MOD)
end

function r_eci_to_eci(T::T_ROT, ::Val{:J2000}, ::Val{:TOD}, JD_UTC::Number)
    # Get the time in TT.
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Compute and return the composed rotation.
    r_MOD_J2000 = r_gcrf_to_mod_fk5(T, JD_TT)
    r_TOD_MOD   = r_mod_to_tod_fk5(T, JD_TT, 0, 0)

    return compose_rotation(r_MOD_J2000, r_TOD_MOD)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:TOD},
    T_ECIf::Val{:J2000},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, JD_UTC, eop_data))
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:TOD},
    T_ECIf::Val{:J2000},
    JD_UTC::Number
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, JD_UTC))
end

#                                J2000 <=> TEME
# ==============================================================================

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:J2000},
    T_ECEFf::Val{:TEME},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    return r_eci_to_eci(T, T_ECIo, T_ECEFf, JD_UTC)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:TEME},
    T_ECEFf::Val{:J2000},
    JD_UTC::Number,
    eop_data::EOPData_IAU1980
)
    return r_eci_to_eci(T, T_ECIo, T_ECEFf, JD_UTC)
end

function r_eci_to_eci(T::T_ROT, ::Val{:J2000}, ::Val{:TEME}, JD_UTC::Number)
    # Get the time in TT.
    JD_TT  = JD_UTCtoTT(JD_UTC)

    # Return the rotation.
    return rGCRFtoTEME(T, JD_TT, 0, 0)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:TEME},
    T_ECIf::Val{:J2000},
    JD_UTC::Number
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, JD_UTC))
end

#                          Between MOD, TOD, and TEME
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::T_ECIs_of_date,
    JD_UTCo::Number,
    T_ECIf::T_ECIs_of_date,
    JD_UTCf::Number,
    eop_data::EOPData_IAU1980
)
    # In this case, we convert origin to GCRF and then convert back to the
    # destination. This is necessary because the user may want to change the
    # epoch.
    r_GCRF_ECIo = r_eci_to_eci(T, T_ECIo,     Val(:GCRF), JD_UTCo, eop_data)
    r_ECIf_GCRF = r_eci_to_eci(T, Val(:GCRF), T_ECIf,     JD_UTCf, eop_data)

    # Return the full rotation.
    return compose_rotation(r_GCRF_ECIo, r_ECIf_GCRF)
end

function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::T_ECIs_of_date,
    JD_UTCo::Number,
    T_ECIf::T_ECIs_of_date,
    JD_UTCf::Number
)
    # In this case, in which we do not have EOP data, we convert origin to J2000
    # and then convert back to the destination. This is necessary because the
    # user may want to change the epoch.
    r_GCRF_ECIo = r_eci_to_eci(T, T_ECIo,      Val(:J2000), JD_UTCo)
    r_ECIf_GCRF = r_eci_to_eci(T, Val(:J2000), T_ECIf,      JD_UTCf)

    # Return the full rotation.
    return compose_rotation(r_GCRF_ECIo, r_ECIf_GCRF)
end

################################################################################
#                           IAU-2006/2010 CIO-based
################################################################################

#                                GCRF <=> CIRS
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:CIRS},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec2rad = π/648000

    # Get the time in TT.
    JD_TT = JD_UTCtoTT(JD_UTC)

    # Get the EOP data related to the desired epoch.
    dX = eop_data.dX(JD_UTC)*arcsec2rad
    dY = eop_data.dY(JD_UTC)*arcsec2rad

    # Compute and return the rotation.
    return r_gcrf_to_cirs_iau2006(T, JD_TT, dX, dY)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:CIRS},
    T_ECIf::Val{:GCRF},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, JD_UTC, eop_data))
end

function r_eci_to_eci(T::T_ROT, ::Val{:GCRF}, ::Val{:CIRS}, JD_UTC::Number)
    # Get the time in TT.
    JD_TT = JD_UTCtoTT(JD_UTC)

    # Compute and return the rotation.
    return r_gcrf_to_cirs_iau2006(T, JD_TT)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:CIRS},
    T_ECIf::Val{:GCRF},
    JD_UTC::Number
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, JD_UTC))
end

#                                 Between CIRS
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:CIRS},
    JD_UTCo::Number,
    T_ECIf::Val{:CIRS},
    JD_UTCf::Number,
    eop_data::EOPData_IAU2000A
)
    # In this case, we convert origin to GCRF and then convert back to the
    # destination. This is necessary because the user may want to change the
    # epoch.
    r_GCRF_ECIo = r_eci_to_eci(T, T_ECIo,     Val(:GCRF), JD_UTCo, eop_data)
    r_ECIf_GCRF = r_eci_to_eci(T, Val(:GCRF), T_ECIf,     JD_UTCf, eop_data)

    # Return the full rotation.
    return compose_rotation(r_GCRF_ECIo, r_ECIf_GCRF)
end

function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:CIRS},
    JD_UTCo::Number,
    T_ECIf::Val{:CIRS},
    JD_UTCf::Number
)
    # In this case, we convert origin to GCRF and then convert back to the
    # destination. This is necessary because the user may want to change the
    # epoch.
    r_GCRF_ECIo = r_eci_to_eci(T, T_ECIo,     Val(:GCRF), JD_UTCo)
    r_ECIf_GCRF = r_eci_to_eci(T, Val(:GCRF), T_ECIf,     JD_UTCf)

    # Return the full rotation.
    return compose_rotation(r_GCRF_ECIo, r_ECIf_GCRF)
end

################################################################################
#                         IAU-2006/2010 equinox-based
################################################################################

#                               GCRF <=> MJ2000
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MJ2000},
    ::Val{:GCRF},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    return r_mj2000_to_gcrf_iau2006(T)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MJ2000},
    ::Val{:GCRF},
    JD_UTC::Number
)
    return r_mj2000_to_gcrf_iau2006(T)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:MJ2000},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    return r_gcrf_to_mj2000_iau2006(T)
end

function r_eci_to_eci(T::T_ROT, ::Val{:GCRF}, ::Val{:MJ2000}, JD_UTC::Number)
    return r_gcrf_to_mj2000_iau2006(T)
end

#                                 GCRF <=> MOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MOD06},
    ::Val{:GCRF},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(T, Val(:MOD06), Val(:GCRF), JD_UTC)
end

function r_eci_to_eci(T::T_ROT, ::Val{:MOD06}, ::Val{:GCRF}, JD_UTC::Number)
    # Get the time in TT.
    JD_TT = JD_UTCtoTT(JD_UTC)

    # Compute and return the composed rotation.
    r_MJ2000_MOD  = r_mod_to_mj2000_iau2006(T, JD_TT)
    r_GCRF_MJ2000 = r_mj2000_to_gcrf_iau2006(T)

    return compose_rotation(r_MJ2000_MOD, r_GCRF_MJ2000)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:MOD06},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(T, Val(:GCRF), Val(:MOD06), JD_UTC)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:MOD06},
    JD_UTC::Number
)
    return inv_rotation(r_eci_to_eci(T, Val(:MOD06), Val(:GCRF), JD_UTC))
end

#                                 GCRF <=> ERS
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:ERS},
    ::Val{:GCRF},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec2rad = π/648000

    # Get the time in TT.
    JD_TT = JD_UTCtoTT(JD_UTC)

    # Obtain the correction of the nutation in obliquity and longitude.
    δΔϵ_2000, δΔΨ_2000 = dEps_dPsi(eop_data, JD_UTC)
    δΔϵ_2000 *= arcsec2rad
    δΔΨ_2000 *= arcsec2rad

    # Compute and return the composed rotation.
    r_MOD_ERS     = r_ers_to_mod_iau2006(T, JD_TT, δΔϵ_2000, δΔΨ_2000)
    r_MJ2000_MOD  = r_mod_to_mj2000_iau2006(T, JD_TT)
    r_GCRF_MJ2000 = r_mj2000_to_gcrf_iau2006(T)

    return compose_rotation(r_MOD_ERS, r_MJ2000_MOD, r_GCRF_MJ2000)
end

function r_eci_to_eci(T::T_ROT, ::Val{:ERS}, ::Val{:GCRF}, JD_UTC::Number)
    # Get the time in TT.
    JD_TT = JD_UTCtoTT(JD_UTC)

    # Compute and return the composed rotation.
    r_MOD_ERS     = r_ers_to_mod_iau2006(T, JD_TT, 0, 0)
    r_MJ2000_MOD  = r_mod_to_mj2000_iau2006(T, JD_TT)
    r_GCRF_MJ2000 = r_mj2000_to_gcrf_iau2006(T)

    return compose_rotation(r_MOD_ERS, r_MJ2000_MOD, r_GCRF_MJ2000)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:ERS},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    return inv_rotation(r_eci_to_eci(T, Val(:ERS), Val(:GCRF), JD_UTC, eop_data))
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:ERS},
    JD_UTC::Number
)
    return inv_rotation(r_eci_to_eci(T, Val(:ERS), Val(:GCRF), JD_UTC))
end

#                                MJ2000 <=> MOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MOD06},
    ::Val{:MJ2000},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(T, Val(:MOD06), Val(:MJ2000), JD_UTC)
end

function r_eci_to_eci(T::T_ROT, ::Val{:MOD06}, ::Val{:MJ2000}, JD_UTC::Number)
    # Get the time in TT.
    JD_TT = JD_UTCtoTT(JD_UTC)

    # Compute the rotation.
    return r_mod_to_mj2000_iau2006(T, JD_TT)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MJ2000},
    ::Val{:MOD06},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(T, Val(:MJ2000), Val(:MOD06), JD_UTC)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MJ2000},
    ::Val{:MOD06},
    JD_UTC::Number
)
    return inv_rotation(r_eci_to_eci(T, Val(:MOD06), Val(:MJ2000), JD_UTC))
end

#                                 MJ2000 <=> ERS
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:ERS},
    ::Val{:MJ2000},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    arcsec2rad = π/648000

    # Get the time in TT.
    JD_TT = JD_UTCtoTT(JD_UTC)

    # Obtain the correction of the nutation in obliquity and longitude.
    δΔϵ_2000, δΔΨ_2000 = dEps_dPsi(eop_data, JD_UTC)
    δΔϵ_2000 *= arcsec2rad
    δΔΨ_2000 *= arcsec2rad

    # Compute and return the composed rotation.
    r_MOD_ERS     = r_ers_to_mod_iau2006(T, JD_TT, δΔϵ_2000, δΔΨ_2000)
    r_MJ2000_MOD  = r_mod_to_mj2000_iau2006(T, JD_TT)

    return compose_rotation(r_MOD_ERS, r_MJ2000_MOD)
end

function r_eci_to_eci(T::T_ROT, ::Val{:ERS}, ::Val{:MJ2000}, JD_UTC::Number)
    # Get the time in TT.
    JD_TT = JD_UTCtoTT(JD_UTC)

    # Compute and return the composed rotation.
    r_MOD_ERS     = r_ers_to_mod_iau2006(T, JD_TT, 0, 0)
    r_MJ2000_MOD  = r_mod_to_mj2000_iau2006(T, JD_TT)

    return compose_rotation(r_MOD_ERS, r_MJ2000_MOD)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MJ2000},
    ::Val{:ERS},
    JD_UTC::Number,
    eop_data::EOPData_IAU2000A
)
    return inv_rotation(r_eci_to_eci(T, Val(:ERS), Val(:MJ2000), JD_UTC, eop_data))
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MJ2000},
    ::Val{:ERS},
    JD_UTC::Number
)
    return inv_rotation(r_eci_to_eci(T, Val(:ERS), Val(:MJ2000), JD_UTC))
end

#                             Between ERS and MOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::T_ECIs_IAU_2006_Equinox_of_date,
    JD_UTCo::Number,
    T_ECIf::T_ECIs_IAU_2006_Equinox_of_date,
    JD_UTCf::Number,
    eop_data::EOPData_IAU2000A
)
    # In this case, we convert origin to GCRF and then convert back to the
    # destination. This is necessary because the user may want to change the
    # epoch.
    r_GCRF_ECIo = r_eci_to_eci(T, T_ECIo,     Val(:GCRF), JD_UTCo, eop_data)
    r_ECIf_GCRF = r_eci_to_eci(T, Val(:GCRF), T_ECIf,     JD_UTCf, eop_data)

    # Return the full rotation.
    return compose_rotation(r_GCRF_ECIo, r_ECIf_GCRF)
end

function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::T_ECIs_IAU_2006_Equinox_of_date,
    JD_UTCo::Number,
    T_ECIf::T_ECIs_IAU_2006_Equinox_of_date,
    JD_UTCf::Number
)
    # In this case, we convert origin to GCRF and then convert back to the
    # destination. This is necessary because the user may want to change the
    # epoch. Notice that, differently from IAU-76/FK5, we can convert to GCRF
    # without using EOP data with minor degradation in precsion.
    r_GCRF_ECIo = r_eci_to_eci(T, T_ECIo,      Val(:GCRF), JD_UTCo)
    r_ECIf_GCRF = r_eci_to_eci(T, Val(:GCRF), T_ECIf,      JD_UTCf)

    # Return the full rotation.
    return compose_rotation(r_GCRF_ECIo, r_ECIf_GCRF)
end
