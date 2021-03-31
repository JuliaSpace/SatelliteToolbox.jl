# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Structures related to reference frames.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export T_ECEFs, T_ECIs, T_ECIs_of_date, T_ECEFs_IAU_2006, T_ECIs_IAU_2006,
       T_ECIs_IAU_2006_Equinox_of_date, T_ROT

"""
    T_ECEFs

Union of all Earth-Centered Earth-Fixed (ECEF) frames supported by the
IAU-76/FK5 theory.

"""
T_ECEFs = Union{Val{:ITRF}, Val{:PEF}}

"""
    T_ECIs

Union of all Earth-Centered Inertial (ECI) frames supported by the IAU-76/FK5
theory.

"""
T_ECIs = Union{Val{:GCRF}, Val{:J2000}, Val{:TOD}, Val{:MOD}, Val{:TEME}}

"""
    T_ECIs_of_date

Union of all *of date* Earth-Centered Inertial (ECI) frames supported by the
IAU-76/FK5 theory.

"""
T_ECIs_of_date = Union{Val{:TOD}, Val{:MOD}, Val{:TEME}}

"""
    T_ECEFs_IAU_2006

Union of all Earth-Centered Earth-Fixed (ECEF) frames supported by IAU-2006/2010
theory.

"""
T_ECEFs_IAU_2006 = Union{Val{:ITRF}, Val{:TIRS}}

"""
    T_ECIs_IAU_2006_CIO

Union of all Earth-Centered Inertial (ECI) frames supported by CIO-based
IAU-2006/2010 theory.

"""
T_ECIs_IAU_2006_CIO = Union{Val{:GCRF}, Val{:CIRS}}

"""
    T_ECIs_IAU_2006_Equinox

Union of all Earth-Centered Inertial (ECI) frames supported by Equinox-based
IAU-2006/2010 theory.

"""
T_ECIs_IAU_2006_Equinox = Union{Val{:GCRF}, Val{:MJ2000}, Val{:MOD06}, Val{:ERS}}

"""
    T_ECIs_IAU_2006

Union of all Earth-Centered Inertial (ECI) frames supported by IAU-2006/2010
theory.

"""
T_ECIs_IAU_2006 = Union{T_ECIs_IAU_2006_CIO, T_ECIs_IAU_2006_Equinox}

"""
    T_ECIs_IAU_2006_Equinox_of_date

Union of all *of date* Earth-Centered Inertial (ECI) frames supported by the
equinox-based IAU-2006/2010 theory.

"""
T_ECIs_IAU_2006_Equinox_of_date = Union{Val{:MOD06}, Val{:ERS}}

"""
    T_ROT

Union of all supported rotation descriptions.

"""
T_ROT = Union{Type{DCM}, Type{Quaternion}}
