# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Structures related to the IERS EOP support.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export EOPData_IAU1980, EOPData_IAU2000A

"""
    EOPData_IAU1980{T}

EOP Data for IAU 1980.

# Fields

* `x, y`: Polar motion with respect to the crust [arcsec].
* `UT1_UTC`: Irregularities of the rotation angle [s].
* `LOD`: Length of day offset [ms].
* `dPsi, dEps`: Celestial pole offsets referred to the model IAU1980 [milliarcsec].
* `*_err`: Errors in the components [same unit as the component].

# Remarks

Each field will be an `AbstractInterpolation` indexed by the Julian Day. Hence,
if one want to obtain, for example, the X component of the polar motion with
respect to the crust at 19 June 2018, the following can be used:

    x[DatestoJD(2018,19,06,0,0,0)]

"""
struct EOPData_IAU1980{T}
    x::T
    y::T
    UT1_UTC::T
    LOD::T
    dPsi::T
    dEps::T

    # Errors
    x_err::T
    y_err::T
    UT1_UTC_err::T
    LOD_err::T
    dPsi_err::T
    dEps_err::T
end

"""
    EOPData_IAU2000A{T}

EOP Data for IAU 2000A.

# Fields

* `x, y`: Polar motion with respect to the crust [arcsec].
* `UT1_UTC`: Irregularities of the rotation angle [s].
* `LOD`: Length of day offset [s].
* `dX, dY`: Celestial pole offsets referred to the model IAU2000A [arcsec].
* `*_err`: Errors in the components [same unit as the component].

# Remarks

Each field will be an `AbstractInterpolation` indexed by the Julian Day. Hence,
if one want to obtain, for example, the X component of the polar motion with
respect to the crust at 19 June 2018, the following can be used:

    x[DatestoJD(2018,19,06,0,0,0)]

"""
struct EOPData_IAU2000A{T}
    x::T
    y::T
    UT1_UTC::T
    LOD::T
    dX::T
    dY::T

    # Errors
    x_err::T
    y_err::T
    UT1_UTC_err::T
    LOD_err::T
    dX_err::T
    dY_err::T
end
