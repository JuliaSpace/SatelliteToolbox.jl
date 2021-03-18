# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Orbit representation as Keplerian elements.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export kepler_to_sv

################################################################################
#                                     API
################################################################################

# Direct getters.
@inline get_epoch(k::KeplerianElements) = k.t
@inline get_a(k::KeplerianElements) = k.a
@inline get_e(k::KeplerianElements) = k.e
@inline get_i(k::KeplerianElements) = k.i
@inline get_Ω(k::KeplerianElements) = k.Ω
@inline get_ω(k::KeplerianElements) = k.ω
@inline get_f(k::KeplerianElements) = k.f

# Getters that require conversions.
@inline get_M(k::KeplerianElements) = f_to_M(k.e, k.f)
@inline get_r(k::KeplerianElements) = kepler_to_rv(k)[1]
@inline get_v(k::KeplerianElements) = kepler_to_rv(k)[2]
@inline get_rv(k::KeplerianElements) = kepler_to_rv(k)

# Conversions.
"""
    kepler_to_sv(k::KeplerianElements)

Convert the Keplerian elements `k` to a state vector.

"""
function kepler_to_sv(k::KeplerianElements{T}) where T
    r_i, v_i = kepler_to_rv(k)
    return OrbitStateVector{T}(t = k.t, r = r_i, v = v_i)
end

################################################################################
#                                  Overloads
################################################################################

getindex(k::KeplerianElements{T}, ::Colon) where T<:Number =
    SVector{7,T}(k.t, k.a, k.e, k.i, k.Ω, k.ω, k.f)

################################################################################
#                                      IO
################################################################################

function show(io::IO, k::KeplerianElements{T}) where T
    compact = get(io, :compact, true)
    epoch_str = sprint(print, k.t, context = :compact => compact)
    jd_str = sprint(print, JDtoDate(DateTime, k.t))
    print(io, "KeplerianElements{", T, "}: Epoch = $epoch_str ($jd_str)")
end

function show(io::IO, mime::MIME"text/plain", k::KeplerianElements{T}) where T
    d2r = 180/π

    # Check if the IO supports color.
    color = get(io, :color, false)

    # Compact printing.
    compact = get(io, :compact, true)

    # Definition of colors that will be used for printing.
    b = color ? _b : ""
    d = color ? _d : ""
    g = color ? _g : ""
    y = color ? _y : ""

    # Convert the data to string.
    date_str  = sprint(print, JDtoDate(DateTime, k.t))
    epoch_str = sprint(print, k.t, context = :compact => compact)
    a_str     = sprint(print, k.a/1000, context = :compact => compact)
    e_str     = sprint(print, k.e, context = :compact => compact)
    i_str     = sprint(print, rad2deg(k.i), context = :compact => compact)
    Ω_str     = sprint(print, rad2deg(k.Ω), context = :compact => compact)
    ω_str     = sprint(print, rad2deg(k.ω), context = :compact => compact)
    f_str     = sprint(print, rad2deg(k.f), context = :compact => compact)

    # Padding to align in the floating point.
    Δepoch = findfirst('.', epoch_str)
    Δepoch === nothing && (Δepoch = length(epoch_str) + 1)

    Δa = findfirst('.', a_str)
    Δa === nothing && (Δa = length(a_str) + 1)

    Δe = findfirst('.', e_str)
    Δe === nothing && (Δe = length(e_str) + 1)

    Δi = findfirst('.', i_str)
    Δi === nothing && (Δi = length(i_str) + 1)

    ΔΩ = findfirst('.', Ω_str)
    ΔΩ === nothing && (ΔΩ = length(Ω_str) + 1)

    Δω = findfirst('.', ω_str)
    Δω === nothing && (Δω = length(ω_str) + 1)

    Δf = findfirst('.', f_str)
    Δf === nothing && (Δf = length(f_str) + 1)

    dp_pos = max(Δepoch, Δa, Δe, Δi, ΔΩ, Δω, Δf)

    epoch_str = " "^(dp_pos - Δepoch) * epoch_str
    a_str     = " "^(dp_pos - Δa) * a_str
    e_str     = " "^(dp_pos - Δe) * e_str
    i_str     = " "^(dp_pos - Δi) * i_str
    Ω_str     = " "^(dp_pos - ΔΩ) * Ω_str
    ω_str     = " "^(dp_pos - Δω) * ω_str
    f_str     = " "^(dp_pos - Δf) * f_str

    max_length = max(length(a_str),
                     length(e_str),
                     length(i_str),
                     length(Ω_str),
                     length(ω_str),
                     length(f_str))

    # Add the units.
    a_str *= " "^(max_length - length(a_str)) * " km"
    i_str *= " "^(max_length - length(i_str)) * " °"
    Ω_str *= " "^(max_length - length(Ω_str)) * " °"
    ω_str *= " "^(max_length - length(ω_str)) * " °"
    f_str *= " "^(max_length - length(f_str)) * " °"

    # Print the Keplerian elements.
    println(io, "KeplerianElements{", T, "}:")
    println(io, "$(b)           Epoch : $(d)" * epoch_str * " (" * date_str * ")");
    println(io, "$(b) Semi-major axis : $(d)" * a_str)
    println(io, "$(b)    Eccentricity : $(d)" * e_str)
    println(io, "$(b)     Inclination : $(d)" * i_str)
    println(io, "$(b)            RAAN : $(d)" * Ω_str)
    println(io, "$(b) Arg. of Perigee : $(d)" * ω_str)
    print(io,   "$(b)    True Anomaly : $(d)" * f_str)
end
