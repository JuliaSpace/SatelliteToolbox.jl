# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   fluxtable.txt
#
#   This file stores the F10.7 in different formats.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

struct _fluxtable_Structure
    F10obs::_space_indices_itp_constant{Float64}
    F10adj::_space_indices_itp_constant{Float64}
end

# Remote file: fluxtable.txt
_fluxtable = @RemoteFile(
    "ftp://ftp.seismo.nrcan.gc.ca/spaceweather/solar_flux/daily_flux_values/fluxtable.txt",
    file="fluxtable.txt",
    updates=:daily
)

# Optional variable that will store the `fluxtable.txt` data.
@OptionalData(
    _fluxtable_data,
    _fluxtable_Structure,
    "Run `init_space_indices()` with `:fluxtable` in `enabled_files` array to initialize required data."
)

################################################################################
#                              Private Functions
################################################################################

"""
    _init_fluxtable(;force_download = false, local_path = nothing)

Initialize the data in the file `fluxtable.txt` by creating `_fluxtable_data`.
The initialization process is composed of:

1. Download the file, if it is necessary;
2. Parse the file;
3. Create the interpolations and the structures.

If the keyword `force_download` is `true`, then the file will always be
downloaded.

The user can also specify a location for the file using the keyword
`local_path`. If it is `nothing`, which is the default, then the file will be
downloaded.

"""
function _init_fluxtable(;force_download = false, local_path = nothing)
    # Update the remote files if no path is given.
    if local_path == nothing
        download(_fluxtable; force = force_download, force_update = true)
        local_path = path(_fluxtable)
    end

    push!(_fluxtable_data,   _parse_fluxtable(local_path))

    return nothing
end

"""
    _parse_fluxtable(path::AbstractString)

Parse the `fluxtable.txt` file in `path` and return an instance of the structure
`_fluxtable_Structure` with the initialize interpolations.

"""
function _parse_fluxtable(path::AbstractString)
    # Allocate raw data.
    JD      = Float64[]
    F10obs  = Float64[]
    F10adj  = Float64[]

    # Store the latest processed Julian day.
    JD_k_1 = zero(Float64)

    open(path) do file
        line_num = 0

        for ln in eachline(file)
            line_num += 1

            # Jump comments at the beginning of the file.
            (line_num in [1,2]) && continue

            # Get the tokens.
            tokens   = split(ln)
            num_toks = length(tokens)

            # Check if the line is valid.
            (num_toks != 7) && @error "Error parsing the line $line_num of fluxtable.txt"

            # We will save only the flux obtained at 20h. This is what
            # NRLMSISE-00 and JB2008 uses.
            fluxtime = parse(Int, tokens[2])
            (fluxtime != 200000) && continue

            # The JD of the data will be computed at noon. Hence, we will be
            # able to use the nearest-neighbor algorithm in the interpolations.
            year  = parse(Int, tokens[1][1:4])
            month = parse(Int, tokens[1][5:6])
            day   = parse(Int, tokens[1][7:8])
            JD_k  = DatetoJD(year, month, day, 12, 0, 0)

            # If the current data is equal to the last one, it means we have a
            # duplicated data. In this case, always use the lastest one.
            if JD_k == JD_k_1
                pop!(JD)
                pop!(F10obs)
                pop!(F10adj)
            end
            JD_k_1 = JD_k

            # Get the raw data.
            push!(JD,      JD_k)
            push!(F10obs,  parse(Float64, tokens[5]))
            push!(F10adj,  parse(Float64, tokens[6]))
        end
    end

    # Create the interpolations for each parameter.
    knots = (JD,)
    itp_F10obs  = interpolate( knots, F10obs,  Gridded(Constant()) )
    itp_F10adj  = interpolate( knots, F10adj,  Gridded(Constant()) )

    return _fluxtable_Structure(itp_F10obs, itp_F10adj)
end
