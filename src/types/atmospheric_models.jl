# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Structures related to the atmospheric models.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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
