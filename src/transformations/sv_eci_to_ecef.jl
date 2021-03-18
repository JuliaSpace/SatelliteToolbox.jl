# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Convert an orbit state vector from an Earth-Centered Inertial (ECI)
#   reference frame to an Earth-Centered, Earth-Fixed (ECEF) reference frame.
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

export svECItoECEF

"""
    svECItoECEF(sv::OrbitStateVector, ECI, ECEF, JD_UTC [, eop_data])

Convert the orbit state vector `sv` from the Earth-Centered Inertial (ECI)
reference frame `ECI` to the Earth-Centered, Earth-Fixed (ECEF) reference frame
at the Julian day `JD_UTC` [UTC]. The `eop_data` may be required depending on
the selection of the input and output reference system. For more information,
see the documentation of the function `rECItoECEF`.

!!! info

    It is assumed that the input velocity and acceleration in `sv` are obtained
    by an observer on the ECI frame. Thus, the output will contain the velocity
    and acceleration as measured by an observer on the ECEF frame.

"""
function svECItoECEF(sv::OrbitStateVector, T_ECI::T_ECIs,
                     T_ECEF::Val{:ITRF}, JD_UTC::Number,
                     eop_data::EOPData_IAU1980)

    # First, convert from ECI to PEF.
    sv_pef = svECItoECEF(sv, T_ECI, Val(:PEF), JD_UTC, eop_data)

    # Now, convert from PEF to ITRF.
    return svECEFtoECEF(sv_pef, Val(:PEF), Val(:ITRF), JD_UTC, eop_data)
end

function svECItoECEF(sv::OrbitStateVector, T_ECI::T_ECIs_IAU_2006,
                     T_ECEF::Val{:ITRF}, JD_UTC::Number,
                     eop_data::EOPData_IAU2000A)

    # First, convert from ECI to TIRS.
    sv_tirs = svECItoECEF(sv, T_ECI, Val(:TIRS), JD_UTC, eop_data)

    # Now, convert from TIRS to ITRF.
    return svECEFtoECEF(sv_tirs, Val(:TIRS), Val(:ITRF), JD_UTC, eop_data)
end

function svECItoECEF(sv::OrbitStateVector,
                     T_ECI::Union{T_ECIs, T_ECIs_IAU_2006},
                     T_ECEF::Union{Val{:PEF},Val{:TIRS}},
                     JD_UTC::Number,
                     eop_data::Union{Nothing,EOPData_IAU1980, EOPData_IAU2000A} = nothing)

    # Get the matrix that converts the ECI to the ECEF.
    if eop_data == nothing
        D = rECItoECEF(DCM, T_ECI, T_ECEF, JD_UTC)
    else
        D = rECItoECEF(DCM, T_ECI, T_ECEF, JD_UTC, eop_data)
    end

    # Since the ECI and ECEF frames have a relative velocity between them, then
    # we must account from it when converting the velocity and acceleration. The
    # angular velocity between those frames is computed using `we` and corrected
    # by the length of day (LOD) parameter of the EOP data, if available.
    ω  = we * ( 1 - (eop_data != nothing ? eop_data.LOD(JD_UTC)/86400 : 0 ) )
    vω = SVector{3}(0,0,ω)

    # Compute the position in the ECEF frame.
    r_ecef = D*sv.r

    # Compute the velocity in the ECEF frame.
    vω_x_r = vω × r_ecef
    v_ecef = D*sv.v - vω_x_r

    # Compute the acceleration in the ECI frame.
    a_ecef = D*sv.a - vω × vω_x_r - 2vω × v_ecef

    return satsv(sv.t, r_ecef, v_ecef, a_ecef)
end
