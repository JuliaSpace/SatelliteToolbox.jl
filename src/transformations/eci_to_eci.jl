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
    r_ecef_to_eci([T,] ECIo, ECIf, jd_utc::Number [, eop_data])
    r_ecef_to_eci([T,] ECIo, jd_utco::Number, ECIf, jd_utcf::Number [, eop_data])

Compute the rotation from an Earth-Centered Inertial (`ECI`) reference frame to
another ECI reference frame. If the origin and destination frame contain only
one *of date* frame, then the first signature is used and `jd_utc` is the epoch
of this frame. On the other hand, if the origin and destination frame contain
two *of date* frame`¹`, e.g. TOD => MOD, then the second signature must be used
in which `jd_utco` is the epoch of the origin frame and `jd_utcf` is the epoch
of the destination frame.

The rotation description that will be used is given by `T`, which can be `DCM`
or `Quaternion`. The origin ECI frame is selected by the input `ECIo` and the
destination ECI frame is selected by the input `ECIf`. The model used to compute
the rotation is specified by the selection of the origin and destination frames.
Currently, there are two models supported: IAU-76/FK5 and IAU-2006 with 2010
conventions (CIO and equinox approaches).

`¹`: TEME is an *of date* frame.

# Rotation description

The rotations that aligns the origin ECI frame with the destination ECI frame
can be described by Direction Cosine Matrices or Quaternions. This is selected
by the parameter `T`.

The possible values are:

- `DCM`: The rotation will be described by a Direction Cosine Matrix.
- `Quaternion`: The rotation will be described by a Quaternion.

If no value is specified, then it falls back to `DCM`.

# Conversion model

The model that will be used to compute the rotation is automatically inferred
given the selection of the origin and destination frames. **Notice that mixing
IAU-76/FK5 and IAU-2006/2010 frames is not supported.**

# ECI Frame

The supported ECI frames for both origin `ECIo` and destination `ECIf` are:

- `TEME()`: ECI will be selected as the True Equator Mean Equinox (TEME)
    reference frame.
- `TOD()`: ECI will be selected as the True of Date (TOD).
- `MOD()`: ECI will be selected as the Mean of Date (MOD).
- `J2000()`: ECI will be selected as the J2000 reference frame.
- `GCRF()`: ECI will be selected as the Geocentric Celestial Reference Frame
    (GCRF).
- `CIRS()`: ECEF will be selected as the Celestial Intermediate Reference System
    (CIRS).
- `ERS()`: ECI will be selected as the Earth Reference System (ERS).
- `MOD06()`: ECI will be selected as the Mean of Date (MOD) according to the
    definition in IAU-2006/2010 theory.
- `MJ2000()`: ECI will be selected as the J2000 mean equatorial frame (MJ2000).

!!! note
    The frames `MOD()` and `MOD06()` are virtually the same. However, we
    selected different names to make clear which theory are being used since
    mixing transformation between frames from IAU-76/FK5 and IAU-2006/2010 must
    be performed with caution.

# EOP Data

The conversion between the frames depends on EOP Data (see
[`get_iers_eop`](@ref) and [`read_iers_eop`](@ref)). If IAU-76/FK5 model is
used, then the type of `eop_data` must be [`EOPData_IAU1980`](@ref). Otherwise,
if IAU-2006/2010 model is used, then the type of `eop_data` must be
[`EOPData_IAU2000A`](@ref). The following table shows the requirements for EOP
data given the selected frames.

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
julia> eop_IAU1980 = get_iers_eop(Val(:IAU1980));

julia> r_eci_to_eci(DCM, GCRF(), J2000(), date_to_jd(1986, 6, 19, 21, 35, 0), eop_IAU1980)
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
  1.0          -4.71326e-12   1.53474e-9
  4.71332e-12   1.0          -3.53979e-9
 -1.53474e-9    3.53979e-9    1.0

julia> r_eci_to_eci(Quaternion, TEME(), GCRF(), date_to_jd(1986, 6, 19, 21, 35, 0), eop_IAU1980)
Quaternion{Float64}:
  + 0.999999 + 1.83013e-5⋅i + 0.000665304⋅j + 0.000665304⋅k

julia> r_eci_to_eci(TOD(), date_to_jd(1986,6,19,21,35,0), TOD(), date_to_jd(1987,5,19,3,0,0), eop_IAU1980)
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 1.0          -0.000224088  -9.73787e-5
 0.000224087   1.0          -5.80065e-6
 9.738e-5      5.77883e-6    1.0

julia> r_eci_to_eci(Quaternion, TOD(), JD_J2000, MOD(), JD_J2000, eop_IAU1980)
Quaternion{Float64}:
  + 1.0 - 1.40025e-5⋅i + 1.34736e-5⋅j + 1.34736e-5⋅k

julia> r_eci_to_eci(J2000(), TEME(), date_to_jd(1986,6,19,21,35,0))
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
  0.999995    0.0030265    0.00133055
 -0.00302645  0.999995    -3.86125e-5
 -0.00133066  3.45854e-5   0.999999

julia> eop_IAU2000A = get_iers_eop(Val(:IAU2000A));

julia> r_eci_to_eci(CIRS(), GCRF(), date_to_jd(1986,6,19,21,35,0), eop_IAU2000A)
3×3 StaticArrays.SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
 0.999999     3.88389e-8  -0.00133066
 7.18837e-9   1.0          3.45897e-5
 0.00133066  -3.45897e-5   0.999999

julia> r_eci_to_eci(Quaternion, CIRS(), GCRF(), date_to_jd(1986,6,19,21,35,0), eop_IAU2000A)
Quaternion{Float64}:
  + 1.0 + 1.72949e-5⋅i + 0.000665332⋅j + 0.000665332⋅k
```
"""
@inline function r_eci_to_eci(
    T_ECIo::T_ECIs,
    T_ECIf::T_ECIs,
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    return r_eci_to_eci(DCM, T_ECIo, T_ECIf, jd_utc, eop_data)
end

@inline function r_eci_to_eci(
    T_ECIo::T_ECIs_IAU_2006,
    T_ECIf::T_ECIs_IAU_2006,
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(DCM, T_ECIo, T_ECIf, jd_utc, eop_data)
end

# Specializations for those cases in which we have two *of dates* frames.
@inline function r_eci_to_eci(
    T_ECIo::T_ECIs_of_date,
    jd_utco::Number,
    T_ECIf::T_ECIs_of_date,
    jd_utcf::Number,
    eop_data::EOPData_IAU1980
)
    return r_eci_to_eci(DCM, T_ECIo, jd_utco, T_ECIf, jd_utcf, eop_data)
end

@inline function r_eci_to_eci(
    T_ECIo::Val{:CIRS},
    jd_utco::Number,
    T_ECIf::Val{:CIRS},
    jd_utcf::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(DCM, T_ECIo, jd_utco, T_ECIf, jd_utcf, eop_data)
end

@inline function r_eci_to_eci(
    T_ECIo::T_ECIs_IAU_2006_Equinox_of_date,
    jd_utco::Number,
    T_ECIf::T_ECIs_IAU_2006_Equinox_of_date,
    jd_utcf::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(DCM, T_ECIo, jd_utco, T_ECIf, jd_utcf, eop_data)
end

# Specializations for those cases that EOP Data is not needed.
@inline function r_eci_to_eci(
    T_ECIo::Val{:J2000},
    T_ECIf::Union{Val{:MOD}, Val{:TOD}, Val{:TEME}},
    jd_utc::Number
)
    return r_eci_to_eci(DCM, T_ECIo, T_ECIf, jd_utc)
end

@inline function r_eci_to_eci(
    T_ECIo::Union{Val{:MOD}, Val{:TOD}, Val{:TEME}},
    T_ECIf::Val{:J2000},
    jd_utc::Number
)
    return r_eci_to_eci(DCM, T_ECIo, T_ECIf, jd_utc)
end

@inline function r_eci_to_eci(
    T_ECIo::T_ECIs_of_date,
    jd_utco::Number,
    T_ECIf::T_ECIs_of_date,
    jd_utcf::Number
)
    return r_eci_to_eci(DCM, T_ECIo, jd_utco, T_ECIf, jd_utcf)
end

@inline function r_eci_to_eci(
    T_ECIo::T_ECIs_IAU_2006,
    T_ECIf::T_ECIs_IAU_2006,
    jd_utc::Number
)
    return r_eci_to_eci(DCM, T_ECIo, T_ECIf, jd_utc)
end

 @inline function r_eci_to_eci(
    T_ECIo::Val{:CIRS},
    jd_utco::Number,
    T_ECIf::Val{:CIRS},
    jd_utcf::Number
 )
    return r_eci_to_eci(DCM, T_ECIo, jd_utco, T_ECIf, jd_utcf)
end

@inline function r_eci_to_eci(
    T_ECIo::T_ECIs_IAU_2006_Equinox_of_date,
    T_ECIf::T_ECIs_IAU_2006_Equinox_of_date,
    jd_utc::Number
)
    return r_eci_to_eci(DCM, T_ECIo, T_ECIf, jd_utc)
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
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    milliarcsec_to_rad = π / 648000000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    δΔϵ_1980 = eop_data.dEps(jd_utc) * milliarcsec_to_rad
    δΔψ_1980 = eop_data.dPsi(jd_utc) * milliarcsec_to_rad

    # In this case, we need to convert GCRF back to PEF and then convert to
    # J2000, which is the same conversion from PEF to GCRF **without** the EOP
    # data.
    #
    # TODO: Can I simplify the rotation from TOD with corrections to TOD without
    # corrections?
    r_mod_gcrf  = r_gcrf_to_mod_fk5(T, jd_tt)
    r_pef_mod   = r_mod_to_pef_fk5(T, jd_ut1, jd_tt, δΔϵ_1980, δΔψ_1980)
    r_mod_pef   = r_pef_to_mod_fk5(T, jd_ut1, jd_tt, 0, 0)
    r_j2000_mod = inv_rotation(r_mod_gcrf)

    # Compose the full rotation.
    return compose_rotation(r_mod_gcrf, r_pef_mod, r_mod_pef, r_j2000_mod)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:J2000},
    T_ECIf::Val{:GCRF},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, jd_utc, eop_data))
end

#                                 GCRF <=> MOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:MOD},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    # Get the time in TT.
    jd_tt = jd_utc_to_tt(jd_utc)

    # Return the rotation.
    return r_gcrf_to_mod_fk5(T, jd_tt)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:MOD},
    T_ECIf::Val{:GCRF},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, jd_utc, eop_data))
end

#                                 GCRF <=> TOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:TOD},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    milliarcsec_to_rad = π / 648000000

    # Get the time in TT.
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    δΔϵ_1980 = eop_data.dEps(jd_utc) * milliarcsec_to_rad
    δΔψ_1980 = eop_data.dPsi(jd_utc) * milliarcsec_to_rad

    # Return the rotation.
    r_mod_gcrf = r_gcrf_to_mod_fk5(T, jd_tt)
    r_tod_mod  = r_mod_to_tod_fk5(T, jd_tt, δΔϵ_1980, δΔψ_1980)

    # Compose the full rotation.
    return compose_rotation(r_mod_gcrf, r_tod_mod)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:TOD},
    T_ECIf::Val{:GCRF},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, jd_utc, eop_data))
end

#                                GCRF <=> TEME
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:TEME},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    milliarcsec_to_rad = π / 648000000

    # Get the time in TT.
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    δΔϵ_1980 = eop_data.dEps(jd_utc) * milliarcsec_to_rad
    δΔψ_1980 = eop_data.dPsi(jd_utc) * milliarcsec_to_rad

    # Return the rotation.
    r_mod_gcrf = r_gcrf_to_mod_fk5(T, jd_tt)
    r_teme_mod = r_mod_to_teme(T, jd_tt, δΔϵ_1980, δΔψ_1980)

    # Compose the full rotation.
    return compose_rotation(r_mod_gcrf, r_teme_mod)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:TEME},
    T_ECIf::Val{:GCRF},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, jd_utc, eop_data))
end

#                                J2000 <=> MOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:J2000},
    ::Val{:MOD},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    milliarcsec_to_rad = π / 648000000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    δΔϵ_1980 = eop_data.dEps(jd_utc) * milliarcsec_to_rad
    δΔψ_1980 = eop_data.dPsi(jd_utc) * milliarcsec_to_rad

    # In this case, we need to convert J2000 back to PEF and then convert to
    # MOD. This is necessary because we need to apply EOP corrections to convert
    # to MOD and just a `r_gcrf_to_mod_fk5` would yield a frame according to the
    # original IAU-76/FK5 theory.
    #
    # TODO: Can I simplify this rotation?
    r_mod_j2000 = r_gcrf_to_mod_fk5(T, jd_tt)
    r_pef_mod   = r_mod_to_pef_fk5(T, jd_ut1, jd_tt, 0, 0)
    r_mod_pef   = r_pef_to_mod_fk5(T, jd_ut1, jd_tt, δΔϵ_1980, δΔψ_1980)

    # Compose the full rotation.
    return compose_rotation(r_mod_j2000, r_pef_mod, r_mod_pef)
end

function r_eci_to_eci(T::T_ROT, ::Val{:J2000}, ::Val{:MOD}, jd_utc::Number)
    # Get the time in TT.
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Compute and return the rotation.
    r_mod_j2000 = r_gcrf_to_mod_fk5(T, jd_tt)

    return r_mod_j2000
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:MOD},
    T_ECIf::Val{:J2000},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, jd_utc, eop_data))
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:MOD},
    T_ECIf::Val{:J2000},
    jd_utc::Number
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, jd_utc))
end

#                                J2000 <=> TOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:J2000},
    ::Val{:TOD},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    milliarcsec_to_rad = π / 648000000

    # Get the time in UT1 and TT.
    jd_ut1 = jd_utc_to_ut1(jd_utc, eop_data)
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    δΔϵ_1980 = eop_data.dEps(jd_utc) * milliarcsec_to_rad
    δΔψ_1980 = eop_data.dPsi(jd_utc) * milliarcsec_to_rad

    # In this case, we need to convert J2000 back to PEF and then convert to
    # TOD. This is necessary because we need to apply EOP corrections to convert
    # to TOD.
    #
    # TODO: Can I simplify this rotation?
    r_mod_j2000 = r_gcrf_to_mod_fk5(T, jd_tt)
    r_pef_mod   = r_mod_to_pef_fk5(T, jd_ut1, jd_tt, 0, 0)
    r_mod_pef   = r_pef_to_mod_fk5(T, jd_ut1, jd_tt, δΔϵ_1980, δΔψ_1980)
    r_tod_mod   = r_mod_to_tod_fk5(T, jd_tt, δΔϵ_1980, δΔψ_1980)

    # Compose the full rotation.
    return compose_rotation(r_mod_j2000, r_pef_mod, r_mod_pef, r_tod_mod)
end

function r_eci_to_eci(T::T_ROT, ::Val{:J2000}, ::Val{:TOD}, jd_utc::Number)
    # Get the time in TT.
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Compute and return the composed rotation.
    r_mod_j2000 = r_gcrf_to_mod_fk5(T, jd_tt)
    r_tod_mod   = r_mod_to_tod_fk5(T, jd_tt, 0, 0)

    return compose_rotation(r_mod_j2000, r_tod_mod)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:TOD},
    T_ECIf::Val{:J2000},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, jd_utc, eop_data))
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:TOD},
    T_ECIf::Val{:J2000},
    jd_utc::Number
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, jd_utc))
end

#                                J2000 <=> TEME
# ==============================================================================

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:J2000},
    T_ECEFf::Val{:TEME},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    return r_eci_to_eci(T, T_ECIo, T_ECEFf, jd_utc)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:TEME},
    T_ECEFf::Val{:J2000},
    jd_utc::Number,
    eop_data::EOPData_IAU1980
)
    return r_eci_to_eci(T, T_ECIo, T_ECEFf, jd_utc)
end

function r_eci_to_eci(T::T_ROT, ::Val{:J2000}, ::Val{:TEME}, jd_utc::Number)
    # Get the time in TT.
    jd_tt  = jd_utc_to_tt(jd_utc)

    # Return the rotation.
    return r_gcrf_to_teme(T, jd_tt, 0, 0)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:TEME},
    T_ECIf::Val{:J2000},
    jd_utc::Number
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, jd_utc))
end

#                          Between MOD, TOD, and TEME
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::T_ECIs_of_date,
    jd_utco::Number,
    T_ECIf::T_ECIs_of_date,
    jd_utcf::Number,
    eop_data::EOPData_IAU1980
)
    # In this case, we convert origin to GCRF and then convert back to the
    # destination. This is necessary because the user may want to change the
    # epoch.
    r_gcrf_ecio = r_eci_to_eci(T, T_ECIo,     Val(:GCRF), jd_utco, eop_data)
    r_ecif_gcrf = r_eci_to_eci(T, Val(:GCRF), T_ECIf,     jd_utcf, eop_data)

    # Return the full rotation.
    return compose_rotation(r_gcrf_ecio, r_ecif_gcrf)
end

function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::T_ECIs_of_date,
    jd_utco::Number,
    T_ECIf::T_ECIs_of_date,
    jd_utcf::Number
)
    # In this case, in which we do not have EOP data, we convert origin to J2000
    # and then convert back to the destination. This is necessary because the
    # user may want to change the epoch.
    r_gcrf_ecio = r_eci_to_eci(T, T_ECIo,      Val(:J2000), jd_utco)
    r_ecif_gcrf = r_eci_to_eci(T, Val(:J2000), T_ECIf,      jd_utcf)

    # Return the full rotation.
    return compose_rotation(r_gcrf_ecio, r_ecif_gcrf)
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
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    milliarcsec_to_rad = π / 648000000

    # Get the time in TT.
    jd_tt = jd_utc_to_tt(jd_utc)

    # Get the EOP data related to the desired epoch.
    dx = eop_data.dX(jd_utc) * milliarcsec_to_rad
    dy = eop_data.dY(jd_utc) * milliarcsec_to_rad

    # Compute and return the rotation.
    return r_gcrf_to_cirs_iau2006(T, jd_tt, dx, dy)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:CIRS},
    T_ECIf::Val{:GCRF},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, jd_utc, eop_data))
end

function r_eci_to_eci(T::T_ROT, ::Val{:GCRF}, ::Val{:CIRS}, jd_utc::Number)
    # Get the time in TT.
    jd_tt = jd_utc_to_tt(jd_utc)

    # Compute and return the rotation.
    return r_gcrf_to_cirs_iau2006(T, jd_tt)
end

@inline function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:CIRS},
    T_ECIf::Val{:GCRF},
    jd_utc::Number
)
    return inv_rotation(r_eci_to_eci(T, T_ECIf, T_ECIo, jd_utc))
end

#                                 Between CIRS
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:CIRS},
    jd_utco::Number,
    T_ECIf::Val{:CIRS},
    jd_utcf::Number,
    eop_data::EOPData_IAU2000A
)
    # In this case, we convert origin to GCRF and then convert back to the
    # destination. This is necessary because the user may want to change the
    # epoch.
    r_gcrf_ecio = r_eci_to_eci(T, T_ECIo,     Val(:GCRF), jd_utco, eop_data)
    r_ecif_gcrf = r_eci_to_eci(T, Val(:GCRF), T_ECIf,     jd_utcf, eop_data)

    # Return the full rotation.
    return compose_rotation(r_gcrf_ecio, r_ecif_gcrf)
end

function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::Val{:CIRS},
    jd_utco::Number,
    T_ECIf::Val{:CIRS},
    jd_utcf::Number
)
    # In this case, we convert origin to GCRF and then convert back to the
    # destination. This is necessary because the user may want to change the
    # epoch.
    r_gcrf_ecio = r_eci_to_eci(T, T_ECIo,     Val(:GCRF), jd_utco)
    r_ecif_gcrf = r_eci_to_eci(T, Val(:GCRF), T_ECIf,     jd_utcf)

    # Return the full rotation.
    return compose_rotation(r_gcrf_ecio, r_ecif_gcrf)
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
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    return r_mj2000_to_gcrf_iau2006(T)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MJ2000},
    ::Val{:GCRF},
    jd_utc::Number
)
    return r_mj2000_to_gcrf_iau2006(T)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:MJ2000},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    return r_gcrf_to_mj2000_iau2006(T)
end

function r_eci_to_eci(T::T_ROT, ::Val{:GCRF}, ::Val{:MJ2000}, jd_utc::Number)
    return r_gcrf_to_mj2000_iau2006(T)
end

#                                 GCRF <=> MOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MOD06},
    ::Val{:GCRF},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(T, Val(:MOD06), Val(:GCRF), jd_utc)
end

function r_eci_to_eci(T::T_ROT, ::Val{:MOD06}, ::Val{:GCRF}, jd_utc::Number)
    # Get the time in TT.
    jd_tt = jd_utc_to_tt(jd_utc)

    # Compute and return the composed rotation.
    r_mj2000_mod  = r_mod_to_mj2000_iau2006(T, jd_tt)
    r_gcrf_mj2000 = r_mj2000_to_gcrf_iau2006(T)

    return compose_rotation(r_mj2000_mod, r_gcrf_mj2000)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:MOD06},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(T, Val(:GCRF), Val(:MOD06), jd_utc)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:MOD06},
    jd_utc::Number
)
    return inv_rotation(r_eci_to_eci(T, Val(:MOD06), Val(:GCRF), jd_utc))
end

#                                 GCRF <=> ERS
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:ERS},
    ::Val{:GCRF},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    milliarcsec_to_rad = π / 648000000

    # Get the time in TT.
    jd_tt = jd_utc_to_tt(jd_utc)

    # Obtain the correction of the nutation in obliquity and longitude.
    δΔϵ_2000, δΔΨ_2000 = deps_dpsi(eop_data, jd_utc)
    δΔϵ_2000 *= milliarcsec_to_rad
    δΔΨ_2000 *= milliarcsec_to_rad

    # Compute and return the composed rotation.
    r_mod_ers     = r_ers_to_mod_iau2006(T, jd_tt, δΔϵ_2000, δΔΨ_2000)
    r_mj2000_mod  = r_mod_to_mj2000_iau2006(T, jd_tt)
    r_gcrf_mj2000 = r_mj2000_to_gcrf_iau2006(T)

    return compose_rotation(r_mod_ers, r_mj2000_mod, r_gcrf_mj2000)
end

function r_eci_to_eci(T::T_ROT, ::Val{:ERS}, ::Val{:GCRF}, jd_utc::Number)
    # Get the time in TT.
    jd_tt = jd_utc_to_tt(jd_utc)

    # Compute and return the composed rotation.
    r_mod_ers     = r_ers_to_mod_iau2006(T, jd_tt, 0, 0)
    r_mj2000_mod  = r_mod_to_mj2000_iau2006(T, jd_tt)
    r_gcrf_mj2000 = r_mj2000_to_gcrf_iau2006(T)

    return compose_rotation(r_mod_ers, r_mj2000_mod, r_gcrf_mj2000)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:ERS},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    return inv_rotation(r_eci_to_eci(T, Val(:ERS), Val(:GCRF), jd_utc, eop_data))
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:GCRF},
    ::Val{:ERS},
    jd_utc::Number
)
    return inv_rotation(r_eci_to_eci(T, Val(:ERS), Val(:GCRF), jd_utc))
end

#                                MJ2000 <=> MOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MOD06},
    ::Val{:MJ2000},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(T, Val(:MOD06), Val(:MJ2000), jd_utc)
end

function r_eci_to_eci(T::T_ROT, ::Val{:MOD06}, ::Val{:MJ2000}, jd_utc::Number)
    # Get the time in TT.
    jd_tt = jd_utc_to_tt(jd_utc)

    # Compute the rotation.
    return r_mod_to_mj2000_iau2006(T, jd_tt)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MJ2000},
    ::Val{:MOD06},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    return r_eci_to_eci(T, Val(:MJ2000), Val(:MOD06), jd_utc)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MJ2000},
    ::Val{:MOD06},
    jd_utc::Number
)
    return inv_rotation(r_eci_to_eci(T, Val(:MOD06), Val(:MJ2000), jd_utc))
end

#                                 MJ2000 <=> ERS
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:ERS},
    ::Val{:MJ2000},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    milliarcsec_to_rad = π / 648000000

    # Get the time in TT.
    jd_tt = jd_utc_to_tt(jd_utc)

    # Obtain the correction of the nutation in obliquity and longitude.
    δΔϵ_2000, δΔΨ_2000 = deps_dpsi(eop_data, jd_utc)
    δΔϵ_2000 *= milliarcsec_to_rad
    δΔΨ_2000 *= milliarcsec_to_rad

    # Compute and return the composed rotation.
    r_mod_ers    = r_ers_to_mod_iau2006(T, jd_tt, δΔϵ_2000, δΔΨ_2000)
    r_mj2000_mod = r_mod_to_mj2000_iau2006(T, jd_tt)

    return compose_rotation(r_mod_ers, r_mj2000_mod)
end

function r_eci_to_eci(T::T_ROT, ::Val{:ERS}, ::Val{:MJ2000}, jd_utc::Number)
    # Get the time in TT.
    jd_tt = jd_utc_to_tt(jd_utc)

    # Compute and return the composed rotation.
    r_mod_ers    = r_ers_to_mod_iau2006(T, jd_tt, 0, 0)
    r_mj2000_mod = r_mod_to_mj2000_iau2006(T, jd_tt)

    return compose_rotation(r_mod_ers, r_mj2000_mod)
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MJ2000},
    ::Val{:ERS},
    jd_utc::Number,
    eop_data::EOPData_IAU2000A
)
    return inv_rotation(r_eci_to_eci(T, Val(:ERS), Val(:MJ2000), jd_utc, eop_data))
end

function r_eci_to_eci(
    T::T_ROT,
    ::Val{:MJ2000},
    ::Val{:ERS},
    jd_utc::Number
)
    return inv_rotation(r_eci_to_eci(T, Val(:ERS), Val(:MJ2000), jd_utc))
end

#                             Between ERS and MOD
# ==============================================================================

function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::T_ECIs_IAU_2006_Equinox_of_date,
    jd_utco::Number,
    T_ECIf::T_ECIs_IAU_2006_Equinox_of_date,
    jd_utcf::Number,
    eop_data::EOPData_IAU2000A
)
    # In this case, we convert origin to GCRF and then convert back to the
    # destination. This is necessary because the user may want to change the
    # epoch.
    r_gcrf_ecio = r_eci_to_eci(T, T_ECIo,     Val(:GCRF), jd_utco, eop_data)
    r_ecif_gcrf = r_eci_to_eci(T, Val(:GCRF), T_ECIf,     jd_utcf, eop_data)

    # Return the full rotation.
    return compose_rotation(r_gcrf_ecio, r_ecif_gcrf)
end

function r_eci_to_eci(
    T::T_ROT,
    T_ECIo::T_ECIs_IAU_2006_Equinox_of_date,
    jd_utco::Number,
    T_ECIf::T_ECIs_IAU_2006_Equinox_of_date,
    jd_utcf::Number
)
    # In this case, we convert origin to GCRF and then convert back to the
    # destination. This is necessary because the user may want to change the
    # epoch. Notice that, differently from IAU-76/FK5, we can convert to GCRF
    # without using EOP data with minor degradation in precsion.
    r_gcrf_ecio = r_eci_to_eci(T, T_ECIo,      Val(:GCRF), jd_utco)
    r_ecif_gcrf = r_eci_to_eci(T, Val(:GCRF), T_ECIf,      jd_utcf)

    # Return the full rotation.
    return compose_rotation(r_gcrf_ecio, r_ecif_gcrf)
end
