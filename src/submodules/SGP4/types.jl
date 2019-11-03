#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Types and structures of SGP4 model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export SGP4_GravCte, SGP4_Structure

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

