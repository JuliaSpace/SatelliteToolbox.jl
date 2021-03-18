# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   General structures related to orbits.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export Orbit, KeplerianElements, OrbitStateVector

"""
    Orbit

Abstract type of an orbit representation.

"""
abstract type Orbit{T} end

"""
    KeplerianElements{T1,T2}

This structure defines the orbit in terms of the Keplerian elements.

# Fields

* `t`: Epoch.
* `a`: Semi-major axis [m].
* `e`: Eccentricity [ ].
* `i`: Inclination [rad].
* `Ω`: Right ascension of the ascending node [rad].
* `ω`: Argument of perigee [rad].
* `f`: True anomaly [rad].

"""
@with_kw_noshow struct KeplerianElements{T} <: Orbit{T}
    t::T
    a::T
    e::T
    i::T
    Ω::T
    ω::T
    f::T
end

function KeplerianElements(t::T1, a::T2, e::T3, i::T4, Ω::T5, ω::T6, f::T7) where
    {T1<:Number, T2<:Number, T3<:Number, T4<:Number, T5<:Number, T6<:Number, T7<:Number}
    T = promote_type(T1, T2, T3, T4, T5, T6, T7)
    return KeplerianElements{T}(t,a,e,i,Ω,ω,f)
end

"""
    OrbitStateVector{T}

Store the state vector representation of an orbit.

# Fields

* `t`: Epoch [Julian Day].
* `r`: Position vector [m].
* `v`: Velocity vector [m/s].
* `a`: Acceleration vector [m/s²].

"""
@with_kw_noshow struct OrbitStateVector{T} <: Orbit{T}
    t::T            = 0
    r::SVector{3,T} = SVector{3,T}(0,0,0)
    v::SVector{3,T} = SVector{3,T}(0,0,0)
    a::SVector{3,T} = SVector{3,T}(0,0,0)
end
