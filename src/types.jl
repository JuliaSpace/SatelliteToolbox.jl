#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Types and structures of SatelliteToolbox.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

################################################################################
#                              Atmospheric Models
################################################################################

#                             Jacchia-Bowman 2008
# ==============================================================================

export JB2008_Output

"""
    JB2008_Output

Output structure of the Jacchia-Bowman 2008.

# Fields

* `nN2`: Number density of N₂ [1/m³].
* `nO2`: Number density of O₂ [1/m³].
* `nO`: Number density of O [1/m³].
* `nAr`: Number density of Ar [1/m³].
* `nHe`: Number density of He [1/m³].
* `nH`: Number density of H [1/m³].
* `rho`: Total density [kg/m³].
* `T_exo`: Exospheric temperature [K].
* `Tz`: Temperature at the selected altitude [K].

"""
@with_kw struct JB2008_Output{T}
    nN2::T
    nO2::T
    nO::T
    nAr::T
    nHe::T
    nH::T
    rho::T
    T_exo::T
    Tz::T
end

#                             Jacchia-Roberts 1971
# ==============================================================================

export JR1971_Output

"""
    JR1971_Output

Output structure of the Jacchia-Roberts 1971 model.

# Fields

* `nN2`: Number density of N₂ [1/m³].
* `nO2`: Number density of O₂ [1/m³].
* `nO`: Number density of O [1/m³].
* `nAr`: Number density of Ar [1/m³].
* `nHe`: Number density of He [1/m³].
* `nH`: Number density of H [1/m³].
* `rho`: Total density [kg/m³].
* `T_exo`: Exospheric temperature [K].
* `Tz`: Temperature at the selected altitude [K].

"""
@with_kw struct JR1971_Output{T}
    nN2::T
    nO2::T
    nO::T
    nAr::T
    nHe::T
    nH::T
    rho::T
    T_exo::T
    Tz::T
end

#                                 NRLMSISE-00
# ==============================================================================

export NRLMSISE00_Structure, NRLMSISE00_Flags, NRLMSISE00_Output

"""
    NRLMSISE00_Flags

Flags to configure NRLMSISE-00.

# Fields

* `output_m_kg`
* `F107_Mean`
* `time_independent`
* `sym_annual`
* `sym_semiannual`
* `asym_annual`
* `asyn_semiannual`
* `diurnal`
* `semidiurnal`
* `daily_ap`
* `all_ut_long_effects`
* `longitudinal`
* `ut_mixed_ut_long`
* `mixed_ap_ut_long`
* `terdiurnal`
* `departures_from_eq`
* `all_tinf_var`
* `all_tlb_var`
* `all_tn1_var`
* `all_s_var`
* `all_tn2_var`
* `all_nlb_var`
* `all_tn3_var`
* `turbo_scale_height`
* `use_ap_array`

"""
@with_kw mutable struct NRLMSISE00_Flags
    output_m_kg::Bool         = false
    F107_Mean::Bool           = true
    time_independent::Bool    = true
    sym_annual::Bool          = true
    sym_semiannual::Bool      = true
    asym_annual::Bool         = true
    asym_semiannual::Bool     = true
    diurnal::Bool             = true
    semidiurnal::Bool         = true
    daily_ap::Bool            = true
    all_ut_long_effects::Bool = true
    longitudinal::Bool        = true
    ut_mixed_ut_long::Bool    = true
    mixed_ap_ut_long::Bool    = true
    terdiurnal::Bool          = true
    departures_from_eq::Bool  = true
    all_tinf_var::Bool        = true
    all_tlb_var::Bool         = true
    all_tn1_var::Bool         = true
    all_s_var::Bool           = true
    all_tn2_var::Bool         = true
    all_nlb_var::Bool         = true
    all_tn3_var::Bool         = true
    turbo_scale_height::Bool  = true
    use_ap_array::Bool        = false
end

"""
    NRLMSISE00_Structure{T}

Structure with the configuration parameters for NRLMSISE-00 model. It can be
created using the function `conf_nrlmsise00`.

"""
@with_kw mutable struct NRLMSISE00_Structure{T}
    # Inputs
    # ======
    year::Int = 0
    doy::Int  = 0
    sec::T    = 0
    alt::T    = 0
    g_lat::T  = 0
    g_long::T = 0
    lst::T    = 0
    f107A::T  = 0
    f107::T   = 0
    ap::T     = 0
    ap_array::Vector{T} = zeros(T, 7)
    flags::NRLMSISE00_Flags = NRLMSISE00_Flags()

    # Auxiliary variables to improve code performance
    # ===============================================
    re::T          = T(6367.088132098377173)
    gsurf::T       = T(0)
    df::T          = T(0)
    dfa::T         = T(0)
    plg::Matrix{T} = zeros(T, 8, 8)
    ctloc::T       = T(0)
    stloc::T       = T(0)
    c2tloc::T      = T(0)
    s2tloc::T      = T(0)
    s3tloc::T      = T(0)
    c3tloc::T      = T(0)

    # In the original source code, it has 4 components, but only 1 is used.
    apt::T         = T(0)
    apdf::T        = T(0)
    dm28::T        = T(0)
    # The original code declared all the `meso_*` vectors as global variables.
    # However, only two values really need to be shared between the functions
    # `gts7` and `gtd7`.
    meso_tn1_5::T  = T(0)
    meso_tgn1_2::T = T(0)
end

"""
    NRLMSISE00_Output

Output structure for NRLMSISE00 model.

# Fields

* `den_N`: Nitrogen number density [U].
* `den_N2`: N₂ number density [U].
* `den_O`: Oxygen number density [U].
* `den_aO`: Anomalous Oxygen number density [U].
* `den_O2`: O₂ number density [U].
* `den_H`: Hydrogen number density [U].
* `den_He`: Helium number density [U].
* `den_Ar`: Argon number density [U].
* `den_Total`: Total mass density \\[T/U] (this value has different meanings for
               routines `gtd7` and `gtd7d`).
* `T_exo`: Exospheric temperature [K].
* `T_alt`: Temperature at the selected altitude [K].
* `flags`: Flags used to compute NRLMSISE-00 model.

Notice that:

* If `flags.output_m_kg` is `false`, then [U] is [cm⁻³] and [T] is [g/cm⁻³].
* If `flags.output_m_kg` is `true`, then [U] is [m⁻³] and [T] is [kg/m⁻³].

# Remarks

Anomalous oxygen is defined as hot atomic oxygen or ionized oxygen that can
become appreciable at high altitudes (`> 500 km`) for some ranges of inputs,
thereby affection drag on satellites and debris. We group these species under
the term **Anomalous Oxygen**, since their individual variations are not
presently separable with the drag data used to define this model component.

"""
@with_kw struct NRLMSISE00_Output{T}
    # Densities.
    den_N::T
    den_N2::T
    den_O::T
    den_aO::T
    den_O2::T
    den_H::T
    den_He::T
    den_Ar::T
    den_Total::T

    # Temperatures.
    T_exo::T
    T_alt::T

    # Flags.
    flags::NRLMSISE00_Flags
end

################################################################################
#                                Gravity Model
################################################################################

export ICGEM, GravityModel_Coefs

"""
    ICGEM

Structure to store the information contained in ICGEM files.

"""
@with_kw struct ICGEM
    # Mandatory keywords.
    product_type::String
    modelname::String
    gravity_constant::Float64
    radius::Float64
    max_degree::Int
    errors::Symbol

    # Optional keywords.
    tide_system::Symbol
    norm::Symbol

    # Coefficients.
    Clm::Matrix{Float64}           = Matrix{Float64}(undef,0,0)
    Slm::Matrix{Float64}           = Matrix{Float64}(undef,0,0)
    sigmaC::Matrix{Float64}        = Matrix{Float64}(undef,0,0)
    sigmaC_formal::Matrix{Float64} = Matrix{Float64}(undef,0,0)
    sigmaS::Matrix{Float64}        = Matrix{Float64}(undef,0,0)
    sigmaS_formal::Matrix{Float64} = Matrix{Float64}(undef,0,0)
end

"""
    GravityModel_Coefs{T}

Structure to store the information about a gravity model.

"""
@with_kw struct GravityModel_Coefs{T}
    name::String
    μ::T
    R0::T
    legendre_norm::Symbol
    n_max::Int
    C::Matrix{T}
    S::Matrix{T}
end

################################################################################
#                                   IERS EOP
################################################################################

export EOPData_IAU1980, EOPData_IAU2000A

"""
    EOPData_IAU1980{T}

EOP Data for IAU 1980.

# Fields

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


################################################################################
#                                    Orbit
################################################################################

export Orbit

"""
    Orbit{T1,T2}

This structure defines the orbit in terms of the Keplerian elements.

# Fields

* `t`: Orbit epoch.
* `a`: Semi-major axis [m].
* `e`: Eccentricity.
* `i`: Inclination [rad].
* `Ω`: Right ascension of the ascending node [rad].
* `ω`: Argument of perigee [rad].
* `f`: True anomaly [rad].

"""
mutable struct Orbit{T1,T2}
    t::T1
    a::T2
    e::T2
    i::T2
    Ω::T2
    ω::T2
    f::T2
end

function Orbit(t::T1, a::T2, e::T3, i::T4, Ω::T5, ω::T6, f::T7) where
    {T1<:Number, T2<:Number, T3<:Number, T4<:Number, T5<:Number, T6<:Number, T7<:Number}

    T = promote_type(T2,T3,T4,T5,T6,T7)

    return Orbit{T1,T}(t,a,e,i,Ω,ω,f)
end

"""
    SatelliteStateVector{T}

Store the state vector of the satellite.

# Fields

* `t`: Epoch [Julian Day].
* `r`: Position vector [m].
* `v`: Velocity vector [m/s].
* `a`: Acceleration vector [m/s²].

"""
@with_kw_noshow mutable struct SatelliteStateVector{T}
    t::T            = 0
    r::SVector{3,T} = SVector{3,T}(0,0,0)
    v::SVector{3,T} = SVector{3,T}(0,0,0)
    a::SVector{3,T} = SVector{3,T}(0,0,0)
end

################################################################################
#                              Orbit Propagators
################################################################################

export OrbitPropagator

"""
    OrbitPropagator{T}

Abstract type of the orbit propagator. Every propagator structure must be a
subtype of this type and must implement the following API functions:

    propagate!(orbp, t::Number)
    propagate!(orbp, t::AbstractVector)
    propagate_to_epoch!(orbp, JD::Number)
    propagate_to_epoch!(orbp, JD::AbstractVector)
    step!(orbp, Δt::Number)

"""
abstract type OrbitPropagator{T} end

#                          Two Body Orbit Propagator
# ==============================================================================

export TwoBody_Structure, OrbitPropagatorTwoBody

"""
    TwoBody_Structure{T}

Low level Two Body orbit propagator structure.

"""
mutable struct TwoBody_Structure{T}
    # Initial parameters.
    epoch::T
    a_0::T
    e_0::T
    i_0::T
    Ω_0::T
    ω_0::T
    M_0::T
    f_0::T
    # Propagation time from epoch.
    Δt::T
    # Mean motion.
    n_0::T
    # Current parameters.
    M_k::T
    f_k::T
    # Standard gravitational parameter of the central body [m^3/s^2].
    μ::T
end

"""
    OrbitPropagatorTwoBody{T} <: OrbitPropagator{T}

Structure that holds the information related to the Two Body orbit propagator.

# Fields

* `orb`: Mean orbital elements (see `Orbit`).
* `tbd`: Structure that stores the Two Body orbit propagator data (see
        `TwoBody_Structure`).

"""
mutable struct OrbitPropagatorTwoBody{T} <: OrbitPropagator{T}
    orb::Orbit{T,T}

    # Two Body orbit propagator related fields.
    tbd::TwoBody_Structure{T}
end


#                             J2 Orbit Propagator
# ==============================================================================

export J2_GravCte, J2_Structure, OrbitPropagatorJ2

"""
    J2_GravCte{T}

Gravitational constants for J2 orbit propagator.

# Fields

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
    J2_Structure{T}

Low level J2 orbit propagator structure.

"""
@with_kw mutable struct J2_Structure{T}
    # Initial mean orbital elements.
    epoch::T
    al_0::T    # Normalized semi-major axis [er].
    n_0::T
    e_0::T
    i_0::T
    Ω_0::T
    ω_0::T
    f_0::T
    M_0::T
    # Propagation time from epoch.
    Δt::T
    # Drag parameters.
    dn_o2::T   # First time derivative of mean motion [rad/s²].
    ddn_o6::T  # Second time derivative of mean motion [rad/s³].
    # Current mean orbital elements.
    al_k::T    # Normalized semi-major axis [er].
    e_k::T
    i_k::T
    Ω_k::T
    ω_k::T
    f_k::T
    M_k::T
    # First-order time-derivative of the orbital elements.
    δa::T
    δe::T
    δΩ::T
    δω::T
    δM_0::T
    # J2 orbit propagator gravitational constants.
    j2_gc::J2_GravCte{T}
end

"""
    OrbitPropagatorJ2{T} <: OrbitPropagator{T}

Structure that holds the information related to the J2 orbit propagator.

# Fields

* `orb`: Mean orbital elements (see `Orbit`).
* `j2d`: Structure that stores the J2 orbit propagator data (see
         `J2_Structure`).

"""
mutable struct OrbitPropagatorJ2{T} <: OrbitPropagator{T}
    orb::Orbit{T,T}

    # J2 orbit propagator related fields.
    j2d::J2_Structure{T}
end

#                        J2 osculating orbit propagator
# ==============================================================================

export J2osc_Strutcture

"""
    J2osc_Structure{T}

Low level J2 osculating orbit propagator structure.

"""
@with_kw mutable struct J2osc_Structure{T}
    # J2 orbit propagator to propagate the mean elements.
    j2d::J2_Structure{T}
    # Propagation time from epoch.
    Δt::T
    # Current osculating Keplerian elements.
    a_k::T
    e_k::T
    i_k::T
    Ω_k::T
    ω_k::T
    f_k::T
    M_k::T
end

#                             J4 orbit propagator
# ==============================================================================

export J4_GravCte, J4_Structure, OrbitPropagatorJ4

"""
    J4_GravCte{T}

Gravitational constants for J4 orbit propagator.

# Fields

* `R0`: Earth equatorial radius [m].
* `μm`: √GM [er/s]^(3/2).
* `J2`: The second gravitational zonal harmonic of the Earth.
* `J4`: The fourth gravitational zonal harmonic of the Earth.

"""
@with_kw struct J4_GravCte{T}
    R0::T
    μm::T
    J2::T
    J4::T
end

"""
    J4_Structure{T}

Low level J4 orbit propagator structure.

"""
@with_kw mutable struct J4_Structure{T}
    # Initial mean orbital elements.
    epoch::T
    al_0::T
    n_0::T
    e_0::T
    i_0::T
    Ω_0::T
    ω_0::T
    f_0::T
    M_0::T
    # Propagation time from epoch.
    Δt::T
    # Drag parameters.
    dn_o2::T   # First time derivative of mean motion [rad/s²].
    ddn_o6::T  # Second time derivative of mean motion [rad/s³].
    # Current mean orbital elements.
    al_k::T
    e_k::T
    i_k::T
    Ω_k::T
    ω_k::T
    f_k::T
    M_k::T
    # First-order time-derivative of the orbital elements.
    δa::T
    δe::T
    δΩ::T
    δω::T
    δM_0::T
    # J4 orbit propagator gravitational constants.
    j4_gc::J4_GravCte{T}
end

"""
    OrbitPropagatorJ4{T} <: OrbitPropagator{T}

Structure that holds the information related to the J4 orbit propagator.

# Fields

* `orb`: Mean orbital elements (see `Orbit`).
* `j4d`: Structure that stores the J4 orbit propagator data (see
         `J4_Structure`).

"""
mutable struct OrbitPropagatorJ4{T} <: OrbitPropagator{T}
    orb::Orbit{T,T}

    # J4 orbit propagator related fields.
    j4d::J4_Structure{T}
end

#                                     SGP4
# ==============================================================================

export OrbitPropagatorSGP4

"""
    OrbitPropagatorSGP4{T} <: OrbitPropagator{T}

Structure that holds the information related to the SGP4 propagator.

# Fields

* `orb`: Mean orbital elements (see `Orbit`).
* `sgp4_gc`: Gravitational contents of the SGP4 algorithm (see `SGP4_GravCte`).
* `sgp4d`: Structure that stores the SGP4 data (see `SGP4_Structure`).
"""
mutable struct OrbitPropagatorSGP4{T} <: OrbitPropagator{T}
    orb::Orbit{T,T}

    # SGP4 related fields.
    sgp4_gc::SGP4_GravCte{T}
    sgp4d::SGP4_Structure{T}
end

################################################################################
#                               Reference Frames
################################################################################

export T_ECEFs, T_ECIs, T_ECIs_of_date, T_ECEFs_IAU_2006, T_ECIs_IAU_2006, T_ROT

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
    T_ECIs_IAU_2006

Union of all Earth-Centered Inertial (ECI) frames supported by IAU-2006/2010
theory.

"""
T_ECIs_IAU_2006 = Union{Val{:GCRF}, Val{:CIRS}}

"""
    T_ROT

Union of all supported rotation descriptions.

"""
T_ROT = Union{Type{DCM}, Type{Quaternion}}

