# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Conversions between representations using Julia built-in system.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

function Base.convert(
    ::Type{KeplerianElements{Tepoch, T}},
    k::KeplerianElements
) where {Tepoch, T}
    return KeplerianElements{Tepoch, T}(k.t, k.a, k.e, k.i, k.Ω, k.ω, k.f)
end

function Base.convert(
    ::Type{OrbitStateVector{Tepoch, T}},
    sv::OrbitStateVector
) where {Tepoch, T}
    return OrbitStateVector{Tepoch, T}(sv.t, sv.r, sv.v, sv.a)
end

function Base.convert(
    ::Type{KeplerianElements},
    sv::OrbitStateVector{Tepoch, T}
) where {Tepoch, T}
    return Base.convert(KeplerianElements{Tepoch, T}, sv)
end

function Base.convert(
    KT::Type{KeplerianElements{Tepoch, T}},
    sv::OrbitStateVector
) where {Tepoch, T}
    k = sv_to_kepler(sv)
    return convert(KT, k)
end

function Base.convert(
    ::Type{OrbitStateVector},
    k::KeplerianElements{Tepoch, T}
) where {Tepoch, T}
    return Base.convert(OrbitStateVector{Tepoch, T}, k)
end

function Base.convert(
    ST::Type{OrbitStateVector{Tepoch, T}},
    k::KeplerianElements
) where {Tepoch, T}
    sv = kepler_to_sv(k)
    return convert(ST, sv)
end
