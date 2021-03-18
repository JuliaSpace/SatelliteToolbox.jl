# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Conversions between representations using Julia built-in system.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

Base.convert(::Type{KeplerianElements{T}}, k::KeplerianElements) where T =
    KeplerianElements{T}(k.t, k.a, k.e, k.i, k.Ω, k.ω, k.f)

Base.convert(::Type{OrbitStateVector{T}}, sv::OrbitStateVector) where T =
    OrbitStateVector{T}(sv.t, sv.r, sv.v, sv.a)

function Base.convert(KT::Type{KeplerianElements{T}}, sv::OrbitStateVector) where T
    k = sv_to_kepler(sv)
    return convert(KT, k)
end

function Base.convert(ST::Type{OrbitStateVector{T}}, k::KeplerianElements) where T
    sv = kepler_to_sv(k)
    return convert(ST, sv)
end
