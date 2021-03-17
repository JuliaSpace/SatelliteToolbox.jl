# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   General structures related to orbits.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

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
