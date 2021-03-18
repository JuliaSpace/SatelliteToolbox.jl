# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions related to the orbit state vector.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export orbsv, sv_to_kepler

################################################################################
#                                     API
################################################################################

# Direct getters.
@inline get_epoch(sv::OrbitStateVector) = sv.t
@inline get_r(sv::OrbitStateVector) = sv.r
@inline get_v(sv::OrbitStateVector) = sv.v
@inline get_rv(sv::OrbitStateVector) = sv.r, sv.r

# Getters that require conversions.
@inline get_a(sv::OrbitStateVector) = get_a(sv_to_kepler(sv))
@inline get_e(sv::OrbitStateVector) = get_e(sv_to_kepler(sv))
@inline get_i(sv::OrbitStateVector) = get_i(sv_to_kepler(sv))
@inline get_Ω(sv::OrbitStateVector) = get_Ω(sv_to_kepler(sv))
@inline get_ω(sv::OrbitStateVector) = get_ω(sv_to_kepler(sv))
@inline get_f(sv::OrbitStateVector) = get_f(sv_to_kepler(sv))
@inline get_M(sv::OrbitStateVector) = get_M(sv_to_kepler(sv))

# Conversions.
"""
    sv_to_kepler(sv::OrbitStateVector)

Convert the state vector `sv` to Keplerian elements represented by an instance
of the structure `KeplerianElements`.

"""
@inline sv_to_kepler(sv::OrbitStateVector) = rv_to_kepler(sv.r, sv.v, sv.t)

################################################################################
#                                Initialization
################################################################################

"""
    orbsv(t::T1, r::AbstractVector{T2}, v::AbstractVector{T3} = [0,0,0], a::AbstractVector{T4} = [0,0,0]) where {T1<:Number, T2<:Number, T3<:Number, T4<:Number}
    orbsv(t::T1, vec::AbstractVector{T2}) where {T1<:Number, T2<:Number}

Create a new satellite state vector (see `OrbitStateVector`) using the position
`r`, velocity `v`, and acceleration `a`. It is also possible to pass a vector
`vec` with the information concatenated.

!!! info

    The vectors `r`, `v`, and `a` must have at least 3 elements. In the case
    more elements are available, they will be neglected. On the other hand, the
    vector `v` must have 6 or 9 dimensions, indicating `[r;v]`, or `[r;v;a]`.

"""
function orbsv(t::T1, r::AbstractVector{T2}, v::AbstractVector{T3},
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

    T = promote_type(T1, T2, T3, T4)
    return OrbitStateVector{T}(t = t, r = r[1:3], v = v[1:3], a = a[1:3])
end

function orbsv(t::T1, vec::AbstractVector{T2}) where {T1<:Number, T2<:Number}
    len_vec = length(vec)

    len_vec ∉ [6,9] && error("The length of input vector must be 6, or 9.")

    if len_vec == 6
        return OrbitStateVector(t = t, r = vec[1:3], v = vec[4:6])
    else
        return OrbitStateVector(t = t, r = vec[1:3], v = vec[4:6], a = vec[7:9])
    end
end

################################################################################
#                                  Overloads
################################################################################

getindex(sv::OrbitStateVector{T}, ::Colon) where T<:Number =
    SVector{9,T}(sv.r..., sv.v..., sv.a...)

################################################################################
#                                      IO
################################################################################

function show(io::IO, sv::OrbitStateVector{T}) where T
    compact = get(io, :compact, true)
    epoch_str = sprint(print, sv.t, context = :compact => compact)
    jd_str = sprint(print, JDtoDate(DateTime, sv.t))
    print(io, "OrbitStateVector{", T, "}: Epoch = $epoch_str ($jd_str)")
end

function show(io::IO, mime::MIME"text/plain", sv::OrbitStateVector{T}) where T
    # Check if the `io` supports colors.
    color = get(io, :color, false)

    # Compact printing.
    compact = get(io, :compact, true)

    b = (color) ? _b : ""
    d = (color) ? _d : ""

    t_str  = sprint(print, sv.t, context = :compact => compact)
    JD_str = sprint(print, JDtoDate(DateTime, sv.t), context = :compact => compact)
    r_str  = sprint(print, sv.r./1000, context = :compact => compact)
    v_str  = sprint(print, sv.v./1000, context = :compact => compact)

    # Add units.
    max_length = max(length(r_str), length(v_str))
    r_str *= " "^(max_length - length(r_str) + 1) * " km"
    v_str *= " "^(max_length - length(v_str) + 1) * " km/s"

    println(io, "OrbitStateVector{", string(T), "}: ")
    println(io, "$b  epoch :$d ", t_str, " (", JD_str, ")")
    println(io, "$b      r :$d ", r_str)
    print(io,   "$b      v :$d ", v_str)

    return nothing
end
