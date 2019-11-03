#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Types and structures of SGP4 model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export SGP4_GravCte, SGP4_Structure, TLE

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
    # Current mean orbital elements.
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
@with_kw_noshow struct TLE
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
