# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Types and structures of SGP4 model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export Sgp4Constants, Sgp4Propagator

"""
    Sgp4Constants{T<:Real}

Gravitational constants for SGP4.

# Fields

- `R0::T`: Earth equatorial radius [km].
- `XKE::T`: 60 ⋅ √(GM / R0^3) [er/min]^(3/2).
- `J2::T`: The second gravitational zonal harmonic of the Earth.
- `J3::T`: The thrid gravitational zonal harmonic of the Earth.
- `J4::T`: The fourth gravitational zonal harmonic of the Earth.
"""
@with_kw struct Sgp4Constants{T}
    R0::T
    XKE::T
    J2::T
    J3::T
    J4::T
end

@with_kw mutable struct Sgp4DeepSpace{T}
    atime::T  = T(0)
    xli::T    = T(0)
    xni::T    = T(0)
    xnq::T    = T(0)
    xfact::T  = T(0)
    ssl::T    = T(0)
    ssg::T    = T(0)
    ssh::T    = T(0)
    sse::T    = T(0)
    ssi::T    = T(0)
    xlamo::T  = T(0)
    omegaq::T = T(0)
    omgdt::T  = T(0)
    gmst::T   = T(0)
    del1::T   = T(0)
    del2::T   = T(0)
    del3::T   = T(0)
    fasx2::T  = T(0)
    fasx4::T  = T(0)
    fasx6::T  = T(0)
    d2201::T  = T(0)
    d2211::T  = T(0)
    d3210::T  = T(0)
    d3222::T  = T(0)
    d4410::T  = T(0)
    d4422::T  = T(0)
    d5220::T  = T(0)
    d5232::T  = T(0)
    d5421::T  = T(0)
    d5433::T  = T(0)
    xnddt::T  = T(0)
    xndot::T  = T(0)
    xldot::T  = T(0)
    zmos::T   = T(0)
    se2::T    = T(0)
    se3::T    = T(0)
    si2::T    = T(0)
    si3::T    = T(0)
    sl2::T    = T(0)
    sl3::T    = T(0)
    sl4::T    = T(0)
    sgh2::T   = T(0)
    sgh3::T   = T(0)
    sgh4::T   = T(0)
    sh2::T    = T(0)
    sh3::T    = T(0)
    zmol::T   = T(0)
    ee2::T    = T(0)
    e3::T     = T(0)
    xi2::T    = T(0)
    xi3::T    = T(0)
    xl2::T    = T(0)
    xl3::T    = T(0)
    xl4::T    = T(0)
    xgh2::T   = T(0)
    xgh3::T   = T(0)
    xgh4::T   = T(0)
    xh2::T    = T(0)
    xh3::T    = T(0)
    pe::T     = T(0)
    pinc::T   = T(0)
    pgh::T    = T(0)
    ph::T     = T(0)
    pl::T     = T(0)
    pgh0::T   = T(0)
    ph0::T    = T(0)
    pe0::T    = T(0)
    pinc0::T  = T(0)
    pl0::T    = T(0)

    isynfl::Bool = false
    iresfl::Bool = false
    ilsz::Bool   = false
end

"""
    Sgp4Propagator{Tepoch, T}

Low level SGP4 proapgator structure.
"""
@with_kw mutable struct Sgp4Propagator{Tepoch, T}
    # TLE parameters.
    epoch::Tepoch
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
    sgp4c::Sgp4Constants{T}
    # SGP4 deep space structure.
    sgp4ds::Sgp4DeepSpace{T}
end
