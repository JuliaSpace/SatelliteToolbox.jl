# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   SOLFSMY.txt
#
#   This file contains the following indices:
#
#       F10, F81c, S10, S81c, M10, M81c, Y10, Y81c
#
#   in which 81c means the 81-day averaged centered value.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################################################################
#                       Private Structures and Variables
################################################################################

"""
    _SOLFSMY_Structure

Structure to store the interpolations of the data in `SOLFSMY.TXT` file.

# Fields

* `F10`: 10.7-cm solar flux [10⁻²² W/(m² Hz)].
* `F81a`: 10.7-cm averaged solar flux, 81-day centered on input time.
* `S10`: EUV index.
* `S81a`: EUV 81-day averaged centered index.
* `M10`: MG2 index scaled to F10.
* `M81a`: MG2 81-day averaged centered index.
* `Y81a`: Solar X-ray & Lya 81-day averaged centered index.
* `Y81a`: Solar X-ray & Lya 81-day averaged centered index.

"""
struct _SOLFSMY_Structure
    F10::_space_indices_itp_constant{Float64,Vector{Float64}}
    F81a::_space_indices_itp_constant{Float64,Vector{Float64}}
    S10::_space_indices_itp_constant{Float64,Vector{Float64}}
    S81a::_space_indices_itp_constant{Float64,Vector{Float64}}
    M10::_space_indices_itp_constant{Float64,Vector{Float64}}
    M81a::_space_indices_itp_constant{Float64,Vector{Float64}}
    Y10::_space_indices_itp_constant{Float64,Vector{Float64}}
    Y81a::_space_indices_itp_constant{Float64,Vector{Float64}}
end

# Remote file: SOLFSMY.TXT
_solfsmy = @RemoteFile(
    "http://sol.spacenvironment.net/jb2008/indices/SOLFSMY.TXT",
    file="SOLFSMY.TXT",
    updates=:daily
)

# Optional variable that will store the `SOLFSMY.TXT` data.
@OptionalData(
    _solfsmy_data,
    _SOLFSMY_Structure,
    "Run `init_space_indices()` with `:solfsmy` in `enabled_files` array to initialize required data."
)

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
    _init_solfsmy(;force_download = false, local_path = nothing)

Initialize the data in the file `SOLFSMY.TXT` by creating `_solfsmy_data`. The
initialization process is composed of:

1. Download the file, if it is necessary;
2. Parse the file;
3. Create the interpolations and the structures.

If the keyword `force_download` is `true`, then the file will always be
downloaded.

The user can also specify a location for the file using the keyword
`local_path`. If it is `nothing`, which is the default, then the file will be
downloaded.

"""
function _init_solfsmy(;force_download = false, local_path = nothing)
    # Update the remote files if no path is given.
    if local_path == nothing
        download(_solfsmy; force = force_download, force_update = true)
        local_path = path(_solfsmy)
    end

    push!(_solfsmy_data,   _parse_solfsmy(local_path))

    return nothing
end

"""
    _parse_solfsmy(path::AbstractString)

Parse the `SOLFSMY.TXT` file in `path` and retur an instance of the structure
`_SOLFSMY_Structure` with the initialized interpolations.

The format of the file `SOLFSMY.TXT` must be:

    YYYY DDD   JulianDay  F10   F81c  S10   S81c  M10   M81c  Y10   Y81c  Ssrc

"""
function _parse_solfsmy(path::AbstractString)
    # Allocate the raw data.
    JD   = Float64[]
    F10  = Float64[]
    F81a = Float64[]
    S10  = Float64[]
    S81a = Float64[]
    M10  = Float64[]
    M81a = Float64[]
    Y10  = Float64[]
    Y81a = Float64[]

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
            push!(F81a, parse(Float64, tokens[ 5]))
            push!(S10 , parse(Float64, tokens[ 6]))
            push!(S81a, parse(Float64, tokens[ 7]))
            push!(M10 , parse(Float64, tokens[ 8]))
            push!(M81a, parse(Float64, tokens[ 9]))
            push!(Y10 , parse(Float64, tokens[10]))
            push!(Y81a, parse(Float64, tokens[11]))
        end
    end

    # Create the interpolations for each parameter.
    #
    # The parameters in `SOLFSMY.TXT` are computed at 12 UT. Hence, if we
    # interpolate by the nearest-neighbor, then it will always return the data
    # related to that day.
    knots    = (JD,)
    itp_F10  = interpolate(knots, F10 , Gridded(Constant()))
    itp_F81a = interpolate(knots, F81a, Gridded(Constant()))
    itp_S10  = interpolate(knots, S10 , Gridded(Constant()))
    itp_S81a = interpolate(knots, S81a, Gridded(Constant()))
    itp_M10  = interpolate(knots, M10 , Gridded(Constant()))
    itp_M81a = interpolate(knots, M81a, Gridded(Constant()))
    itp_Y10  = interpolate(knots, Y10 , Gridded(Constant()))
    itp_Y81a = interpolate(knots, Y81a, Gridded(Constant()))

    return _SOLFSMY_Structure(
        itp_F10,
        itp_F81a,
        itp_S10,
        itp_S81a,
        itp_M10,
        itp_M81a,
        itp_Y10,
        itp_Y81a
    )
end
