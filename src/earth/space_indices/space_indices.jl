#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   This file contains functions to retrieve space indices for many models.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export init_space_indices

################################################################################
#                              Private Structures
################################################################################

"""
Structure to store the interpolations of the data in `DTCFILE.TXT` file.

# Fields

* `DstΔTc`: Temperature variation due to Dst [K].

"""
struct _DTCFILE_Structure{T}
    DstΔTc::T
end

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

################################################################################
#                      Remote Files and Global Variables
################################################################################

# Remote file: SOLFSMY.TXT
# ==============================================================================
#
# This file contains the following indices:
#
#   F10, F81c, S10, S81c, M10, M81c, Y10, Y81c
#
# in which 81c means the 81-day averaged centered value.

solfsmy = @RemoteFile(
    "http://sol.spacenvironment.net/jb2008/indices/SOLFSMY.TXT",
    file="SOLFSMY.TXT",
    updates=:daily
   )

# Remote file: DTCFILE.TXT
# ==============================================================================
#
# This file contains the `ΔTc` variation caused by `Dst` for each hour of every
# day.

dtcfile = @RemoteFile(
    "http://sol.spacenvironment.net/jb2008/indices/DTCFILE.TXT",
    file="DTCFILE.TXT",
    updates=:daily
   )

# Global variable to store the processed data in the remote files.
@OptionalData solfsmy_data _SOLFSMY_Structure "Run `init_space_indices()` to initialize the space indices structures."
@OptionalData dtcfile_data _DTCFILE_Structure "Run `init_space_indices()` to initialize the space indices structures."

################################################################################
#                               Public Functions
################################################################################

"""
    function init_space_indices(;force_download = false)

Initialize all space indices. The initialization process is composed of:

1. Download all the files, if is necessary;
2. Parse all the files;
3. Create the interpolations and the structures.

If the keyword `force_download` is `true`, then the files will always be
downloaded.

"""
function init_space_indices(;force_download = false)
    # Update the remote files.
    download(dtcfile; force = force_download)
    download(solfsmy; force = force_download)

    # Parse each remote file.
    push!(dtcfile_data, _parse_dtcfile(path(dtcfile)))
    push!(solfsmy_data, _parse_solfsmy(path(solfsmy)))

    nothing
end

#                                   Getters
# ==============================================================================

# Create a function `get_S` for each `S` field in `_DTCFILE_Structure`.
for field in fieldnames(_DTCFILE_Structure)
    field_name = string(field)
    func_name  = "get_" * field_name
    func       = Symbol("get_" * string(field))
    qfield     = Meta.quot(field)

    @eval begin
        export $func
"""
    $($func_name)(JD::Number)

Get the value of the index `$($field_name)` at Julian Day `JD`.

This function requires the initialization of the variable `dtcfile_data`.
Otherwise, an exception will be raised. To initialize it, run
`init_space_indices()`.

"""
        ($func)(JD::Number) = getfield(get(dtcfile_data), $qfield)(JD)
    end
end

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

This function requires the initialization of the variable `solfsmy_data`.
Otherwise, an exception will be raised. To initialize it, run
`init_space_indices()`.

"""
        ($func)(JD::Number) = getfield(get(solfsmy_data), $qfield)(JD)
    end
end

################################################################################
#                              Private Functions
################################################################################

"""
    function _parse_dtcfile(path::AbstractString)

Parse the `DTCFILE.TXT` file in `path` and returning an instance of the
structure `_DTCFILE_Structure` with the initialized interpolations.

The format of the file `DTCFILE.TXT` must be:

    DTC YYYY DOY DTC_0h DTC_1h DTC_2h ... DTC_22h DTC_23h

in which `DOY` is the day of the year and `DTC_Xh` is the `ΔTc` at hour `X`.

"""
function _parse_dtcfile(path::AbstractString)
    # Allocate the raw data.
    JD     = Float64[]
    DstΔTc = Float64[]

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
            if (num_toks != 27) || (tokens[1] != "DTC")
                @warn "Line $line_num of file DTCFILE.TXT has invalid format!"
                continue
            end

            # Compute the Julian day at midnight of the day.
            year = parse(Int64, tokens[2])
            doy  = parse(Float64, tokens[3])

            JD_0h = DatetoJD(year, 1, 1, 0, 0, 0) - 1 + doy

            # Parse the data.
            for k = 1:24
                DstΔTc_h = parse(Float64, tokens[k+3])
                push!(JD,  JD_0h + (k-1)/24)
                push!(DstΔTc, DstΔTc_h)
            end

        end
    end

    # Create the interpolations for each parameter.
    knots      = (JD,)
    itp_DstΔTc = interpolate(knots, DstΔTc, Gridded(Linear()))

    _DTCFILE_Structure(itp_DstΔTc)
end

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
