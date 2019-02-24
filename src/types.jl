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

export Orbit, TLE

"""
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

# Fields

* `name`: Name of the satellite.

## First line

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

## Second line

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

export OrbitPropagator

"""
Abstract type of the orbit propagator. Every propagator structure must be a
subtype of this type and must implement the following API functions:

    function propagate!(orbp, t::Number)
    function propagate!(orbp, t::AbstractVector)
    function propagate_to_epoch!(orbp, JD::Number)
    function propagate_to_epoch!(orbp, JD::AbstractVector)
    function step!(orbp, Δt::Number)

"""
abstract type OrbitPropagator{T} end

#                          Two Body Orbit Propagator
# ==============================================================================

export TwoBody_Structure, OrbitPropagatorTwoBody

"""
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
Structure that holds the information related to the Two Body orbit propagator.

# Fields

* `orb`: Current orbit (see `Orbit`).
* `tbd`: Structure that stores the Two Body orbit propagator data (see
        `TwoBody_Structure`).

"""
mutable struct OrbitPropagatorTwoBody{T} <: OrbitPropagator{T}
    orb::Orbit{T,T,T,T,T,T,T}

    # Two Body orbit propagator related fields.
    tbd::TwoBody_Structure{T}
end


#                             J2 Orbit Propagator
# ==============================================================================

export J2_GravCte, J2_Structure, OrbitPropagatorJ2

"""
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
Low level J2 orbit propagator structure.

"""
@with_kw mutable struct J2_Structure{T}
    # Initial orbit parameters.
    epoch::T
    al_0::T    # Normalized semi-major axis [er].
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
    al_k::T    # Normalized semi-major axis [er].
    e_k::T
    i_k::T
    Ω_k::T
    ω_k::T
    M_k::T
    n_k::T
    f_k::T
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
Structure that holds the information related to the J2 orbit propagator.

# Fields

* `orb`: Current orbit (see `Orbit`).
* `j2d`: Structure that stores the J2 orbit propagator data (see
         `J2_Structure`).

"""
mutable struct OrbitPropagatorJ2{T} <: OrbitPropagator{T}
    orb::Orbit{T,T,T,T,T,T,T}

    # J2 orbit propagator related fields.
    j2d::J2_Structure{T}
end

#                                     SGP4
# ==============================================================================

export SGP4_GravCte, SGP4_Structure, OrbitPropagatorSGP4

"""
Gravitational constants for SGP4.

# Fields

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

@with_kw mutable struct SGP4_DeepSpace{T<:Real}
    atime::T  = 0.0
    xli::T    = 0.0
    xni::T    = 0.0
    xnq::T    = 0.0
    xfact::T  = 0.0
    ssl::T    = 0.0
    ssg::T    = 0.0
    ssh::T    = 0.0
    sse::T    = 0.0
    ssi::T    = 0.0
    xlamo::T  = 0.0
    omegaq::T = 0.0
    omgdt::T  = 0.0
    gmst::T   = 0.0
    del1::T   = 0.0
    del2::T   = 0.0
    del3::T   = 0.0
    fasx2::T  = 0.0
    fasx4::T  = 0.0
    fasx6::T  = 0.0
    d2201::T  = 0.0
    d2211::T  = 0.0
    d3210::T  = 0.0
    d3222::T  = 0.0
    d4410::T  = 0.0
    d4422::T  = 0.0
    d5220::T  = 0.0
    d5232::T  = 0.0
    d5421::T  = 0.0
    d5433::T  = 0.0
    xnddt::T  = 0.0
    xndot::T  = 0.0
    xldot::T  = 0.0
    zmos::T   = 0.0
    se2::T    = 0.0
    se3::T    = 0.0
    si2::T    = 0.0
    si3::T    = 0.0
    sl2::T    = 0.0
    sl3::T    = 0.0
    sl4::T    = 0.0
    sgh2::T   = 0.0
    sgh3::T   = 0.0
    sgh4::T   = 0.0
    sh2::T    = 0.0
    sh3::T    = 0.0
    zmol::T   = 0.0
    ee2::T    = 0.0
    e3::T     = 0.0
    xi2::T    = 0.0
    xi3::T    = 0.0
    xl2::T    = 0.0
    xl3::T    = 0.0
    xl4::T    = 0.0
    xgh2::T   = 0.0
    xgh3::T   = 0.0
    xgh4::T   = 0.0
    xh2::T    = 0.0
    xh3::T    = 0.0
    pe::T     = 0.0
    pinc::T   = 0.0
    pgh::T    = 0.0
    ph::T     = 0.0
    pl::T     = 0.0
    pgh0::T   = 0.0
    ph0::T    = 0.0
    pe0::T    = 0.0
    pinc0::T  = 0.0
    pl0::T    = 0.0

    isynfl::Bool = false
    iresfl::Bool = false
    ilsz::Bool   = false
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
    θ²::T
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
    # Selected algorithm.
    algorithm::Symbol
    # SGP4 gravitational constants.
    sgp4_gc::SGP4_GravCte{T}
    # SGP4 deep space structure.
    sgp4_ds::SGP4_DeepSpace{T}
end

"""
Structure that holds the information related to the SGP4 propagator.

# Fields

* `orb`: Current orbit (see `Orbit`).
* `sgp4_gc`: Gravitational contents of the SGP4 algorithm (see `SGP4_GravCte`).
* `sgp4d`: Structure that stores the SGP4 data (see `SGP4_Structure`).

"""
mutable struct OrbitPropagatorSGP4{T} <: OrbitPropagator{T}
    orb::Orbit{T,T,T,T,T,T,T}

    # SGP4 related fields.
    sgp4_gc::SGP4_GravCte{T}
    sgp4d::SGP4_Structure{T}
end

################################################################################
#                               Reference Frames
################################################################################

export T_ECEFs, T_ECIs, T_ECIs_of_date, T_ECEFs_IAU_2006, T_ECIs_IAU_2006, T_ROT

"""
Union of all Earth-Centered Earth-Fixed (ECEF) frames supported by the
IAU-76/FK5 theory.

"""
T_ECEFs = Union{Type{Val{:ITRF}}, Type{Val{:PEF}}}

"""
Union of all Earth-Centered Inertial (ECI) frames supported by the IAU-76/FK5
theory.

"""
T_ECIs = Union{Type{Val{:GCRF}},
               Type{Val{:J2000}},
               Type{Val{:TOD}},
               Type{Val{:MOD}},
               Type{Val{:TEME}}}

"""
Union of all *of date* Earth-Centered Inertial (ECI) frames supported by the
IAU-76/FK5 theory.

"""
T_ECIs_of_date = Union{Type{Val{:TOD}},
                       Type{Val{:MOD}},
                       Type{Val{:TEME}}}
"""
Union of all Earth-Centered Earth-Fixed (ECEF) frames supported by IAU-2006/2010
theory.

"""
T_ECEFs_IAU_2006 = Union{Type{Val{:ITRF}}, Type{Val{:TIRS}}}

"""
Union of all Earth-Centered Inertial (ECI) frames supported by IAU-2006/2010
theory.

"""
T_ECIs_IAU_2006 = Union{Type{Val{:GCRF}}, Type{Val{:CIRS}}}

"""
Union of all supported rotation descriptions.

"""
T_ROT = Union{Type{DCM}, Type{Quaternion}}

