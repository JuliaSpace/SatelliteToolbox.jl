# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   DTCFILE.txt
#
#   This file stores the temperature variation `ΔTc` caused by the `Dst`.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################################################################
#                       Private Structures and Variables
################################################################################

"""
    _DTCFILE_Structure

Structure to store the interpolations of the data in `DTCFILE.TXT` file.

# Fields

* `DstΔTc`: Temperature variation due to Dst [K].

"""
struct _DTCFILE_Structure
    DstΔTc::_space_indices_itp_linear{Float64}
end

# Remote file: DTCFILE.TXT
_dtcfile = @RemoteFile(
    "http://sol.spacenvironment.net/jb2008/indices/DTCFILE.TXT",
    file="DTCFILE.TXT",
    updates=:daily
)

# Optional variable that will store the `DTCFILE.TXT` data.
@OptionalData(
    _dtcfile_data,
    _DTCFILE_Structure,
    "Run `init_space_indices()` with `:dtcfile` in `enabled_files` array to initialize required data."
)

################################################################################
#                               Public Functions
################################################################################

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

This function requires the initialization of the variable `_dtcfile_data`.
Otherwise, an exception will be raised. To initialize it, run
`init_space_indices()`.

"""
        ($func)(JD::Number) = getfield(get(_dtcfile_data), $qfield)(JD)
    end
end

################################################################################
#                              Private Functions
################################################################################

"""
    _init_dctfile(;force_download = false, local_path = nothing)

Initialize the data in the file `DTCFILE.TXT` by creating `_dtcfile_data`. The
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
function _init_dtcfile(;force_download = false, local_path = nothing)
    # Update the remote files if no path is given.
    if local_path == nothing
        download(_dtcfile; force = force_download, force_update = true)
        local_path = path(_dtcfile)
    end

    push!(_dtcfile_data,   _parse_dtcfile(local_path))

    return nothing
end

"""
    _parse_dtcfile(path::AbstractString)

Parse the `DTCFILE.TXT` file in `path` and return an instance of the structure
`_DTCFILE_Structure` with the initialized interpolations.

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
            year = parse(Int,     tokens[2])
            doy  = parse(Float64, tokens[3])

            JD_0h = date_to_jd(year, 1, 1, 0, 0, 0) - 1 + doy

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

    return _DTCFILE_Structure(itp_DstΔTc)
end
