#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Types and structures of SatelliteToolbox.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-05-14: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

################################################################################
#                             Earth Gravity Model
################################################################################

export EGM_Coefs

"""
Structure to store the EGM coefficients.

"""
struct EGM_Coefs{T1,T2,T3}
    C::Matrix{T1}
    S::Matrix{T1}
    μ::T2
    R0::T3
end

################################################################################
#                                   IERS EOP
################################################################################

export EOPData_IAU1980, EOPData_IAU2000A

"""
EOP Data for IAU 1980. The fields are described as follows:

* `x, y`: Polar motion with respect to the crust [arcsec].
* `UT1_UTC`: Irregularities of the rotation angle [s].
* `LOD`: Length of day offset [s].
* `dPsi, dEps`: Celestial pole offsets referred to the model IAU1980 [arcsec].
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
EOP Data for IAU 2000A. The files are described as follows:

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


################################################################################
#                                    Orbit
################################################################################

export Orbit, TLE

"""
This structure defines the orbit in terms of the Keplerian elements.

* `t`: Orbit epoch.
* `a`: Semi-major axis [m].
* `e`: Eccentricity.
* `i`: Inclination [rad].
* `Ω`: Right ascension of the ascending node [rad].
* `ω`: Argument of perigee [rad].
* `f`: True anomaly [rad].

"""
mutable struct Orbit{T1,T2,T3,T4,T5,T6,T7}
    t::T1
    a::T2
    e::T3
    i::T4
    Ω::T5
    ω::T6
    f::T7
end

"""
This structure contains the same elements of the TLE with the same units.

* `name`: Name of the satellite.

# First line

* `sat_num`: Satellite number.
* `classification`: Classification ('U', 'C', or 'S').
* `int_designator`: International designator.
* `epoch_year`: Epoch year (two digits).
* `epoch_day`: Epoch day (day + fraction of the day).
* `epoch`: The epoch represented in Julian Day.
* `dn_o2`: 1st time derivative of mean motion / 2 [rev/day²].
* `ddn_o6`: 2nd time derivative of mean motion / 6 [rev/day³].
* `bstar`: B* drag term.
* `elem_set_number`: Element set number.
* `checksum_l1`: Checksum of the line 1 (modulo 10).

# Second line

* `i`: Inclination [deg].
* `Ω`: Right ascension of the ascending node [deg].
* `e`: Eccentricity.
* `ω`: Argument of perigee [deg].
* `M`: Mean anomaly [deg].
* `n`: Mean motion [rev/day].
* `rev_num`: Revolution number at epoch [rev].
* `checksum_l2`: Checksum of the line 2 (modulo 10).

"""
@with_kw struct TLE
    name::String

    # First line
    # ==========
    sat_num::Int
    classification::Char
    int_designator::String
    epoch_year::Int
    epoch_day::Float64
    epoch::Float64
    dn_o2::Float64
    ddn_o6::Float64
    bstar::Float64
    elem_set_number::Int
    checksum_l1::Int

    # Second line
    # ===========

    i::Float64
    Ω::Float64
    e::Float64
    ω::Float64
    M::Float64
    n::Float64
    rev_num::Int
    checksum_l2
end

################################################################################
#                              Orbit Propagators
################################################################################

#                          Two Body Orbit Propagator
# ==============================================================================

export TwoBody_Structure

"""
Low level Two Body orbit propagator structure.

"""
mutable struct TwoBody_Structure{T}
    # Initial parameters.
    epoch::T
    n_0::T
    e_0::T
    i_0::T
    Ω_0::T
    ω_0::T
    M_0::T
    # Propagation time from epoch.
    Δt::T
    # Auxiliary parameters.
    a::T
    # Current parameters.
    M_k::T
    f_k::T
    # Standard gravitational parameter of the central body [m^3/s^2].
    μ::T
end

"""
Structure that holds the information related to the Two Body orbit propagator.

* `orb`: Current orbit (see `Orbit`).
* `tbd`: Structure that stores the Two Body orbit propagator data (see
        `TwoBody_Structure`).

"""
mutable struct OrbitPropagatorTwoBody{T}
    orb::Orbit{T,T,T,T,T,T,T}

    # Two Body orbit propagator related fields.
    tbd::TwoBody_Structure{T}
end


#                             J2 Orbit Propagator
# ==============================================================================

export J2_GravCte, J2_Structure, OrbitPropagatorJ2

"""
Gravitational constants for J2 orbit propagator.

* `R0`: Earth equatorial radius [m].
* `μm`: √GM [er/s]^(3/2).
* `J2`: The second gravitational zonal harmonic of the Earth.

"""
@with_kw struct J2_GravCte{T}
    R0::T
    μm::T
    J2::T
end

"""
Low level J2 orbit propagator structure.

"""
@with_kw mutable struct J2_Structure{T}
    # Orbit parameters.
    epoch::T
    a_0::T
    n_0::T
    e_0::T
    i_0::T
    Ω_0::T
    ω_0::T
    M_0::T
    # Propagation time from epoch.
    Δt::T
    # Drag parameters.
    dn_o2::T   # First time derivative of mean motion [rad/s²].
    ddn_o6::T  # Second time derivative of mean motion [rad/s³].
    # Current parameters.
    a_k::T
    e_k::T
    i_k::T
    Ω_k::T
    ω_k::T
    M_k::T
    n_k::T
    f_k::T
    # Useful constants to decrease the computational burden.
    C1::T
    C2::T
    C3::T
    C4::T
    # J2 orbit propagator gravitational constants.
    j2_gc::J2_GravCte{T}
end

"""
Structure that holds the information related to the J2 orbit propagator.

* `orb`: Current orbit (see `Orbit`).
* `j2d`: Structure that stores the J2 orbit propagator data (see
         `J2_Structure`).

"""
mutable struct OrbitPropagatorJ2{T}
    orb::Orbit{T,T,T,T,T,T,T}

    # J2 orbit propagator related fields.
    j2d::J2_Structure{T}
end

#                                     SGP4
# ==============================================================================

export SGP4_GravCte, SGP4_Structure, OrbitPropagatorSGP4

"""
Gravitational constants for SGP4.

* `R0`: Earth equatorial radius [km].
* `XKE`: √GM [er/s]^(3/2).
* `J2`: The second gravitational zonal harmonic of the Earth.
* `J3`: The thrid gravitational zonal harmonic of the Earth.
* `J4`: The fourth gravitational zonal harmonic of the Earth.

"""
@with_kw struct SGP4_GravCte{T<:Real}
    R0::T
    XKE::T
    J2::T
    J3::T
    J4::T
end

"""
Low level SGP4 structure.

"""
@with_kw mutable struct SGP4_Structure{T<:Real}
    # TLE parameters.
    epoch::T
    n_0::T
    e_0::T
    i_0::T
    Ω_0::T
    ω_0::T
    M_0::T
    bstar::T
    # Propagation time from epoch.
    Δt::T
    # Current parameters.
    a_k::T
    e_k::T
    i_k::T
    Ω_k::T
    ω_k::T
    M_k::T
    n_k::T
    # Parameters related with the orbit.
    all_0::T
    nll_0::T
    # Useful constants to decrease the computational burden.
    AE::T
    QOMS2T::T
    β_0::T
    ξ::T
    η::T
    sin_i_0::T
    θ::T
    θ2::T
    θ3::T
    θ4::T
    A_30::T
    k_2::T
    k_4::T
    C1::T
    C3::T
    C4::T
    C5::T
    D2::T
    D3::T
    D4::T
    dotM::T
    dotω::T
    dotΩ1::T
    dotΩ::T
    # Others.
    isimp::Bool
    # SGP4 gravitational constants.
    sgp4_gc::SGP4_GravCte{T}
end

"""
Structure that holds the information related to the SGP4 propagator.

* `orb`: Current orbit (see `Orbit`).
* `sgp4_gc`: Gravitational contents of the SGP4 algorithm (see `SGP4_GravCte`).
* `sgp4d`: Structure that stores the SGP4 data (see `SGP4_Structure`).

"""
mutable struct OrbitPropagatorSGP4{T}
    orb::Orbit{T,T,T,T,T,T,T}

    # SGP4 related fields.
    sgp4_gc::SGP4_GravCte{T}
    sgp4d::SGP4_Structure{T}
end

################################################################################
#                               Reference Frames
################################################################################

export T_ECEFs, T_ECIs, T_ECIs_of_date

"""
Union of all Earth-Centered Earth-Fixed (ECEF) frames supported.

"""
T_ECEFs = Union{Type{Val{:ITRF}}, Type{Val{:PEF}}}

"""
Union of all Earth-Centered Inertial (ECI) frames supported.

"""
T_ECIs = Union{Type{Val{:GCRF}},
               Type{Val{:J2000}},
               Type{Val{:TOD}},
               Type{Val{:MOD}},
               Type{Val{:TEME}}}

"""
Union of all *of date* Earth-Centered Inertial (ECI) frames supported.

"""
T_ECIs_of_date = Union{Type{Val{:TOD}},
                       Type{Val{:MOD}},
                       Type{Val{:TEME}}}
