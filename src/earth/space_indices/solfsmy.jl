#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   SOLFSMY.txt
#
#   This file contains the following indices:
#
#       F10, F81c, S10, S81c, M10, M81c, Y10, Y81c
#
#   in which 81c means the 81-day averaged centered value.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

################################################################################
#                       Private Structures and Variables
################################################################################

"""
Structure to store the interpolations of the data in `SOLFSMY.TXT` file.

# Fields

* `F10`: 10.7-cm solar flux [10⁻²² W/(m² Hz)].
* `F10ₐ`: 10.7-cm averaged solar flux, 81-day centered on input time.
* `S10`: EUV index.
* `S10ₐ`: EUV 81-day averaged centered index.
* `XM10`: MG2 index scaled to F10.
* `XM10ₐ`: MG2 81-day averaged centered index.
* `Y10ₐ`: Solar X-ray & Lya 81-day averaged centered index.
* `Y10ₐ`: Solar X-ray & Lya 81-day averaged centered index.

"""
struct _SOLFSMY_Structure{T}
    F10::T
    F10ₐ::T
    S10::T
    S10ₐ::T
    M10::T
    M10ₐ::T
    Y10::T
    Y10ₐ::T
end

# Remote file: SOLFSMY.TXT
_solfsmy = @RemoteFile(
    "http://sol.spacenvironment.net/jb2008/indices/SOLFSMY.TXT",
    file="SOLFSMY.TXT",
    updates=:daily
   )

# Optional variable that will store the `SOLFSMY.TXT` data.
@OptionalData _solfsmy_data _SOLFSMY_Structure "Run `init_space_indices()` to initialize the space indices structures."

################################################################################
#                               Public Functions
################################################################################

#                                   Getters
# ==============================================================================

# Create a function `get_S` for each `S` field in `_SOLFSMY_Structure`.
for field in fieldnames(_SOLFSMY_Structure)
    field_name = string(field)
    func_name  = "get_" * field_name
    func       = Symbol("get_" * string(field))
    qfield     = Meta.quot(field)

    @eval begin
        export $func
"""
    $($func_name)(JD::Number)

Get the value of the index `$($field_name)` at Julian Day `JD`.

This function requires the initialization of the variable `_solfsmy_data`.
Otherwise, an exception will be raised. To initialize it, run
`init_space_indices()`.

"""
        ($func)(JD::Number) = getfield(get(_solfsmy_data), $qfield)(JD)
    end
end

################################################################################
#                              Private Functions
################################################################################

"""
    function _parse_solfsmy(path::AbstractString)

Parse the `SOLFSMY.TXT` file in `path` and returning an instance of the
structure `_SOLFSMY_Structure` with the initialized interpolations.

The format of the file `SOLFSMY.TXT` must be:

    YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c  Ssrc

"""
function _parse_solfsmy(path::AbstractString)
    # Allocate the raw data.
    JD   = Float64[]
    F10  = Float64[]
    F10ₐ = Float64[]
    S10  = Float64[]
    S10ₐ = Float64[]
    M10  = Float64[]
    M10ₐ = Float64[]
    Y10  = Float64[]
    Y10ₐ = Float64[]

    # Read the raw data in the file.
    open(path) do file
        line_num = 0

        for ln in eachline(file)
            line_num += 1

            # Ignore comments.
            (ln[1] == '#') && continue

            tokens   = split(ln)
            num_toks = length(tokens)

            # Ignore blank lines.
            (num_toks == 0) && continue

            # Ignore and print warning about lines with bad format.
            if num_toks != 12
                @warn "Line $line_num of file SOLFSMY.TXT has invalid format!"
                continue
            end

            # Parse data.
            push!(JD  , parse(Float64, tokens[ 3]))
            push!(F10 , parse(Float64, tokens[ 4]))
            push!(F10ₐ, parse(Float64, tokens[ 5]))
            push!(S10 , parse(Float64, tokens[ 6]))
            push!(S10ₐ, parse(Float64, tokens[ 7]))
            push!(M10 , parse(Float64, tokens[ 8]))
            push!(M10ₐ, parse(Float64, tokens[ 9]))
            push!(Y10 , parse(Float64, tokens[10]))
            push!(Y10ₐ, parse(Float64, tokens[11]))
        end
    end

    # Create the interpolations for each parameter.
    #
    # The parameters in `SOLFSMY.TXT` are computed at 12 UT. Hence, if we
    # interpolate by the nearest-neighbor, then it will always return the data
    # related to that day.
    knots    = (JD,)
    itp_F10  = interpolate(knots, F10 , Gridded(Constant()))
    itp_F10ₐ = interpolate(knots, F10ₐ, Gridded(Constant()))
    itp_S10  = interpolate(knots, S10 , Gridded(Constant()))
    itp_S10ₐ = interpolate(knots, S10ₐ, Gridded(Constant()))
    itp_M10  = interpolate(knots, M10 , Gridded(Constant()))
    itp_M10ₐ = interpolate(knots, M10ₐ, Gridded(Constant()))
    itp_Y10  = interpolate(knots, Y10 , Gridded(Constant()))
    itp_Y10ₐ = interpolate(knots, Y10ₐ, Gridded(Constant()))

    _SOLFSMY_Structure(itp_F10, itp_F10ₐ,
                       itp_S10, itp_S10ₐ,
                       itp_M10, itp_M10ₐ,
                       itp_Y10, itp_Y10ₐ)
end
