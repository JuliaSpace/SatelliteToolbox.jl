# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Convert an orbit state vector from an Earth-Centered, Earth-Fixed (ECEF)
#   reference frame to a Earth-Centered Inertial (ECI) reference frame.
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

export svECEFtoECI

"""
    svECEFtoECI(sv::OrbitStateVector, ECEF, ECI, JD_UTC [, eop_data])

Convert the orbit state vector `sv` from the Earth-Centered, Earth-Fixed (ECEF)
reference frame `ECEF` to the Earth-Centered Inertial (ECI) reference frame at
the Julian day `JD_UTC` [UTC]. The `eop_data` may be required depending on the
selection of the input and output reference system. For more information,
see the documentation of the function `r_ecef_to_eci`.

!!! info

    It is assumed that the input velocity and acceleration in `sv` are obtained
    by an observer on the ECEF frame. Thus, the output will contain the velocity
    and acceleration as measured by an observer on the ECI frame.

"""
function svECEFtoECI(sv::OrbitStateVector, T_ECEF::Val{:ITRF},
                     T_ECI::T_ECIs, JD_UTC::Number, eop_data::EOPData_IAU1980)

    # First, convert from the ITRF to PEF.
    sv_pef = sv_ecef_to_ecef(sv, Val(:ITRF), Val(:PEF), JD_UTC, eop_data)

    # Now, convert from PEF to ECI.
    return svECEFtoECI(sv_pef, Val(:PEF), T_ECI, JD_UTC, eop_data)
end

function svECEFtoECI(sv::OrbitStateVector, T_ECEF::Val{:ITRF},
                     T_ECI::T_ECIs_IAU_2006, JD_UTC::Number,
                     eop_data::EOPData_IAU2000A)

    # First, convert from the ITRF to TIRS.
    sv_tirs = sv_ecef_to_ecef(sv, Val(:ITRF), Val(:TIRS), JD_UTC, eop_data)

    # Now, convert from TIRS to ECI.
    return svECEFtoECI(sv_tirs, Val(:TIRS), T_ECI, JD_UTC, eop_data)
end

function svECEFtoECI(sv::OrbitStateVector,
                     T_ECEF::Union{Val{:PEF},Val{:TIRS}},
                     T_ECI::Union{T_ECIs, T_ECIs_IAU_2006}, JD_UTC::Number,
                     eop_data::Union{Nothing, EOPData_IAU1980, EOPData_IAU2000A} = nothing)

    # Get the matrix that converts the ECEF to the ECI.
    if eop_data == nothing
        D = r_ecef_to_eci(DCM, T_ECEF, T_ECI, JD_UTC)
    else
        D = r_ecef_to_eci(DCM, T_ECEF, T_ECI, JD_UTC, eop_data)
    end

    # Since the ECI and ECEF frames have a relative velocity between them, then
    # we must account from it when converting the velocity and acceleration. The
    # angular velocity between those frames is computed using `we` and corrected
    # by the length of day (LOD) parameter of the EOP data, if available.
    ω  = we * ( 1 - (eop_data != nothing ? eop_data.LOD(JD_UTC)/86400 : 0 ) )
    vω = SVector{3}(0,0,ω)

    # Compute the position in the ECI frame.
    r_eci = D*sv.r

    # Compute the velocity in the ECI frame.
    vω_x_r = vω × sv.r
    v_eci = D*(sv.v + vω_x_r )

    # Compute the acceleration in the ECI frame.
    a_eci = D*(sv.a + vω × vω_x_r + 2vω × sv.v)

    return orbsv(sv.t, r_eci, v_eci, a_eci)
end
