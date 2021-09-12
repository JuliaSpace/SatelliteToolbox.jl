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
abstract type Orbit{Tepoch, T} end

"""
    KeplerianElements{Tepoch, T}

This structure defines the orbit in terms of the Keplerian elements.

# Fields

- `t::Tepoch`: Epoch.
- `a::T`: Semi-major axis [m].
- `e::T`: Eccentricity [ ].
- `i::T`: Inclination [rad].
- `Ω::T`: Right ascension of the ascending node [rad].
- `ω::T`: Argument of perigee [rad].
- `f::T`: True anomaly [rad].
"""
@with_kw_noshow struct KeplerianElements{Tepoch, T} <: Orbit{Tepoch, T}
    t::Tepoch
    a::T
    e::T
    i::T
    Ω::T
    ω::T
    f::T
end

function KeplerianElements(
    t::Number,
    a::T1,
    e::T2,
    i::T3,
    Ω::T4,
    ω::T5,
    f::T6
) where {T1<:Number, T2<:Number, T3<:Number, T4<:Number, T5<:Number, T6<:Number}
    T = promote_type(T1, T2, T3, T4, T5, T6)
    return KeplerianElements{typeof(t), T}(t, a, e, i, Ω, ω, f)
end

"""
    OrbitStateVector{Tepoch, T}

Store the state vector representation of an orbit.

# Fields

- `t::Tepoch`: Epoch [Julian Day].
- `r::SVector{3, T}`: Position vector [m].
- `v::SVector{3, T}`: Velocity vector [m/s].
- `a::SVector{3, T}`: Acceleration vector [m/s²].
"""
@with_kw_noshow struct OrbitStateVector{Tepoch, T} <: Orbit{Tepoch, T}
    t::Tepoch        = 0
    r::SVector{3, T} = SVector{3, T}(0, 0, 0)
    v::SVector{3, T} = SVector{3, T}(0, 0, 0)
    a::SVector{3, T} = SVector{3, T}(0, 0, 0)
end
