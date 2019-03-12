#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions related to the state vector.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export satsv

"""
    function satsv(t::T1, r::AbstractVector{T2}, v::AbstractVector{T3} = [0,0,0], a::AbstractVector{T4} = [0,0,0]) where {T1<:Number, T2<:Number, T3<:Number, T4<:Number}
    function satsv(t::T1, vec::AbstractVector{T2}) where {T1<:Number, T2<:Number}

Create a new satellite state vector (see `SatelliteStateVector`) using the
position `r`, velocity `v`, and acceleration `a`. It is also possible to pass a
vector `vec` with the information concatenated.

!!! info

    The vectors `r`, `v`, and `a` must have at least 3 elements. In the case
    more elements are available, they will be neglected. On the other hand, the
    vector `v` must have 6 or 9 dimensions, indicating `[r;v]`, or `[r;v;a]`.

"""
function satsv(t::T1, r::AbstractVector{T2}, v::AbstractVector{T3},
               a::AbstractVector{T4} = [0,0,0]) where {T1<:Number, T2<:Number,
                                                       T3<:Number, T4<:Number}

    len_r = length(r)
    len_v = length(v)
    len_a = length(a)

    len_r  < 3 && error("The vector `r` must have at least 3 elements.")
    len_v  < 3 && error("The vector `v` must have at least 3 elements.")
    len_a  < 3 && error("The vector `a` must have at least 3 elements.")
    len_r != 3 && @warn("Only the first 3 elements of the vector `r` will be used!")
    len_v != 3 && @warn("Only the first 3 elements of the vector `v` will be used!")
    len_a != 3 && @warn("Only the first 3 elements of the vector `a` will be used!")

    T = promote_type(T1,T2,T3,T4)
    return SatelliteStateVector{T}(t = t, r = r[1:3], v = v[1:3], a = a[1:3])
end

function satsv(t::T1, vec::AbstractVector{T2}) where {T1<:Number, T2<:Number}
    len_vec = length(vec)

    len_vec ∉ [6,9] && error("The length of input vector must be 6, or 9.")

    if len_vec == 6
        return SatelliteStateVector(t = t, r = vec[1:3], v = vec[4:6])
    else
        return SatelliteStateVector(t = t, r = vec[1:3], v = vec[4:6], a = vec[7:9])
    end
end

################################################################################
#                                  Overloads
################################################################################

getindex(sv::SatelliteStateVector{T}, ::Colon) where T<:Number =
    SVector{9,T}(sv.r..., sv.v..., sv.a...)
