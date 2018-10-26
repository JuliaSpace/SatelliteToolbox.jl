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

"""
Structure to store the interpolations of the data in WDC files.

# Fields

* `Kp`: Kp index.
* `Ap`: Ap index.

"""
struct _WDC_Structure
    Kp::AbstractExtrapolation
    Ap::AbstractExtrapolation
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

# Remote files: *.wdc
# ==============================================================================
#
# This set contains all the remote files in the directory:
#
#   ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/wdc/
#

wdcfiles = RemoteFileSet(".wdc files", Dict{Symbol,RemoteFile}())

# Global variable to store the processed data in the remote files.
@OptionalData solfsmy_data _SOLFSMY_Structure "Run `init_space_indices()` to initialize the space indices structures."
@OptionalData dtcfile_data _DTCFILE_Structure "Run `init_space_indices()` to initialize the space indices structures."
@OptionalData wdc_data     _WDC_Structure     "Run `init_space_indices()` to initialize the space indices structures."

################################################################################
#                               Public Functions
################################################################################

"""
    function init_space_indices(;force_download = false, dtcfile_path = nothing, solfsmy_path = nothing)

Initialize all space indices. The initialization process is composed of:

1. Download all the files, if is necessary;
2. Parse all the files;
3. Create the interpolations and the structures.

If the keyword `force_download` is `true`, then the files will always be
downloaded.

The user can also specify the location for each required file to retrieve the
space indices. In this case, the download will not be performed. The following
keywords can be used for this:

* `dtcfile_path`: Path to `DTCFILE.TXT`.
* `solfsmy_path`: Path to `SOLFSMY.TXT`.

For the WDC files, which contains the information about `Kp` and `Ap` indices,
the user can select what is the oldest year in which the data will be downloaded
by the keyword `wdcfiles_oldest_year`. By default, it will download the data
from 3 previous years.

"""
function init_space_indices(;force_download = false, dtcfile_path = nothing,
                             solfsmy_path = nothing, wdcfiles_dir = nothing,
                             wdcfiles_oldest_year = year(now())-3)

    # Update the remote files if no path is given.
    if dtcfile_path == nothing
        download(dtcfile; force = force_download)
        dtcfile_path = path(dtcfile)
    end

    if solfsmy_path == nothing
        download(solfsmy; force = force_download)
        solfsmy_path = path(solfsmy)
    end

    years     = Int[]
    filepaths = String[]

    if wdcfiles_dir == nothing
        _prepare_wdc_remote_files(wdcfiles_oldest_year)
        download(wdcfiles)

        # Get the files available and sort them by the year.
        for (sym,wdcfile) in wdcfiles.files
            #
            # The year must not be obtained by the data inside the file,
            # because it contains only 2 digits and will break in 2032.
            # We will obtain the year by the symbol of the remote file. The
            # symbol name is:
            #
            #       kpYYYY
            #
            # where `YYYY` is the year.
            push!(years, parse(Int, String(sym)[3:6]))
            push!(filepaths, path(wdcfile))
        end
    else
        # If the user provided a directory, check what files are available.
        # Notice that the name must be the same as the ones online.
        for (root, dirs, files) in walkdir(wdcfiles_dir)
            for file in files
                if occursin(r"^kp[1-2][0-9][0-9][0-9].wdc$", file)
                    year = parse(Int, file[3:6])

                    # Check if the year is not older than the oldest year.
                    if year >= wdcfiles_oldest_year
                        @info "Found WDC file `$file` related to the year `$year`."
                        push!(filepaths, joinpath(root, file))
                        push!(years,     year)
                    end
                end
            end
        end
    end

    p = sortperm(years)

    # Parse each remote file.
    push!(dtcfile_data, _parse_dtcfile(dtcfile_path))
    push!(solfsmy_data, _parse_solfsmy(solfsmy_path))
    push!(wdc_data,     _parse_wdcfiles(filepaths[p], years[p]))

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

"""
    function get_Kp(JD::Number)

Return the Kp index at Julian Day `JD`.

"""
function get_Kp(JD::Number)
    Kp_day = get(SatelliteToolbox.wdc_data).Kp(JD)

    # Get the hour of the day and return the appropriate Ap.
    y, m, d, h, min, sec = JDtoDate(JD)

    return Kp_day[ floor(Int, h/3) + 1 ]
end

"""
    function get_Ap(JD::Number; mean::Tuple{Int} = (), daily = false)

Return the Ap index.

If `mean` is a tuple of two integers `(hi, hf)`, then the average between `hi`
and `hf` previous hours will be computed.

If `mean` is empty and `daily` is `true`, then the day average will be computed.

If `mean` keyword is empty, and `daily` keyword is `false`, then the Ap at
Julian day `JD` will be computed.

By default, `mean` is empty and `daily` is `false`.

"""
function get_Ap(JD::Number; mean::Tuple = (), daily = false)
    # Check if we must compute the mean of previous hours.
    if isempty(mean)
        Ap_day = get(SatelliteToolbox.wdc_data).Ap(JD)

        # Check if we must compute the daily mean.
        if daily
            return sum(Ap_day)/8
        else
            # Get the hour of the day and return the appropriate Ap.
            y, m, d, h, min, sec = JDtoDate(JD)

            return Ap_day[ floor(Int, h/3) + 1 ]
        end
    else
        # Check the inputs.
        (length(mean) != 2) && @error "The keyword `mean` must be empty or a tuple with exactly 2 integers."
        hi = mean[1]
        hf = mean[2]
        (hi > hf) && @error "The first argument of the keyword `mean` must be lower than the second."

        # Assemble the vector with the previous hours that will be averaged.
        hv = hi:3:hf

        # Compute the mean.
        Ap_sum = 0
        for h in hv
            Ap_sum += get_Ap(JD - h/24; mean = (), daily = false)
        end

        return Ap_sum/length(hv)
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
            year = parse(Int,     tokens[2])
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

function _parse_wdcfiles(filepaths::Vector{String}, years::Vector{Int})
    # Allocate the raw data.
    JD = Float64[]
    Kp = Vector{Float64}[]
    Ap = Vector{Int}[]

    for (filepath, year) in zip(filepaths, years)

        open(filepath) do file
            # Read each line.
            for ln in eachline(file)
                # Get the Julian Day.
                month = parse(Int, ln[3:4])
                day   = parse(Int, ln[5:6])

                # The JD of the data will be computed at noon. Hence, we will be
                # able to use the nearest-neighbor algorithm in the
                # interpolations.
                JD_k  = DatetoJD(year, month, day, 12, 0, 0)

                # Get the vector of Kps and Aps.
                Ap_k = zeros(Int,    8)
                Kp_k = zeros(Float64,8)

                for i = 1:8
                    Kp_k[i] = parse(Int, ln[2(i-1) + 13:2(i-1) + 14])/10
                    Ap_k[i] = parse(Int, ln[3(i-1) + 32:3(i-1) + 34])
                end

                # Add data to the vector.
                push!(JD, JD_k)
                push!(Kp, Kp_k)
                push!(Ap, Ap_k)
            end
        end
    end

    # Create the interpolations for each parameter.
    knots    = (JD,)

    # Create the interpolations.
    itp_Kp = extrapolate(interpolate(knots, Kp, Gridded(Constant())), Flat())
    itp_Ap = extrapolate(interpolate(knots, Ap, Gridded(Constant())), Flat())

    _WDC_Structure(itp_Kp, itp_Ap)
end

function _prepare_wdc_remote_files(oldest_year::Number)
    # Get the current year.
    current_year = year(now())

    # If `oldest_year` is greated than current year, then consider only the
    # current year.
    (oldest_year > current_year) && (oldest_year = current_year)

    # For the current year, we must update the remote file every day. Otherwise,
    # we do not need to update at all.
    for y = oldest_year:current_year
		filename = "kp$y"
        sym = Symbol(filename)
        file_y = @RemoteFile("ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/wdc/$filename.wdc",
                             file="$filename.wdc",
                             updates= (y == current_year) ? :daily : :never)

        merge!(wdcfiles.files, Dict(sym => file_y))
    end

    nothing
end
