# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   WDCFILES.txt
#
#   The WDC files contains the Kp and Ap indices.
#
#   For more information, see:
#
#       https://www.gfz-potsdam.de/en/kp-index/
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

################################################################################
#                       Private Structures and Variables
################################################################################

"""
    _WDC_Structure

Structure to store the interpolations of the data in WDC files.

# Fields

* `Kp`: Kp index.
* `Ap`: Ap index.

"""
struct _WDC_Structure
    Kp::_space_indices_itp_constant{SVector{8,Float64}}
    Ap::_space_indices_itp_constant{SVector{8,Float64}}
end

# Remote files: *.wdc
#
# This set will be configured in the function `init_space_indices()`.
_wdcfiles = RemoteFileSet(".wdc files", Dict{Symbol,RemoteFile}())

# Optional variable that will store the WDC data.
@OptionalData(
    _wdc_data,
    _WDC_Structure,
    "Run `init_space_indices()` with `:wdcfiles` in `enabled_files` array to initialize required data."
)

################################################################################
#                               Public Functions
################################################################################

export get_Kp, get_Ap

#                                   Getters
# ==============================================================================

"""
    get_Kp(JD::Number)

Return the Kp index at Julian Day `JD`.

"""
function get_Kp(JD::Number)
    Kp_day = get(SatelliteToolbox._wdc_data).Kp(JD)

    # Get the hour of the day and return the appropriate Kp.
    y, m, d, h, min, sec = jd_to_date(JD)

    return Kp_day[ floor(Int, h/3) + 1 ]
end

"""
    get_Ap(JD::Number; mean::Tuple{Int} = (), daily = false)

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
        Ap_day = get(SatelliteToolbox._wdc_data).Ap(JD)

        # Check if we must compute the daily mean.
        if daily
            return sum(Ap_day)/8
        else
            # Get the hour of the day and return the appropriate Ap.
            y, m, d, h, min, sec = jd_to_date(JD)

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
    _init_wdcfiles(;force_download = false, local_dir = nothing, wdcfiles_oldest_year = year(now())-3)

Initialize the data in the WDC files by creating `_wdcfiles_data`. The
initialization process is composed of:

1. Download the files, if it is necessary;
2. Parse the files;
3. Create the interpolations and the structures.

If the keyword `force_download` is `true`, then the files will always be
downloaded.

The user can also specify a location for the directory with the WDC files using
the keyword `local_dir`. If it is `nothing`, which is the default, then the
file will be downloaded.

The user can select what is the oldest year in which the data will be downloaded
by the keyword `wdcfiles_oldest_year`. By default, it will download the data
from 3 previous years.

The user can select what is the newest year in which the data will be downloaded
by the keyword `wdcfiles_newest_year`. It it is `nothing`, which is the default,
then it is set to the current year.

"""
function _init_wdcfiles(;
    force_download = false,
    local_dir = nothing,
    wdcfiles_oldest_year = year(now())-3,
    wdcfiles_newest_year = nothing
)
    years     = Int[]
    filepaths = String[]

    if local_dir == nothing
        (wdcfiles_newest_year == nothing) && (wdcfiles_newest_year = year(now()))
        _prepare_wdc_remote_files(wdcfiles_oldest_year, wdcfiles_newest_year)
        download(_wdcfiles; force = force_download, force_update = true)

        # Get the files available and sort them by the year.
        for (sym, wdcfile) in _wdcfiles.files
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

            # Since we are getting those files from FTP, the timestamp is not
            # updated. Hence, we execute `touch` here to update the timestamp.
            touch(path(wdcfile))
        end
    else
        # If the user provided a directory, check what files are available.
        # Notice that the name must be the same as the ones online.
        for (root, dirs, files) in walkdir(local_dir)
            for file in files
                if occursin(r"^kp[1-2][0-9][0-9][0-9].wdc$", file)
                    year = parse(Int, file[3:6])

                    # Check if the year is not older than the oldest year.
                    if year >= wdcfiles_oldest_year
                        @info "Found WDC file `$file` related to the year `$year`."
                        push!(filepaths, joinpath(root, file))
                        push!(years, year)
                    end
                end
            end
        end
    end

    p = sortperm(years)

    push!(_wdc_data, _parse_wdcfiles(filepaths[p], years[p]))

    return nothing
end

"""
    _parse_wdcfiles(filepaths::Vector{String}, years::Vector{Int})

Parse the WDC files with paths in `filepaths` related to the years in `years`.

**Notice that the files must be sorted by the year!**

"""
function _parse_wdcfiles(filepaths::Vector{String}, years::Vector{Int})
    # Allocate the raw data.
    JD = Float64[]
    Kp = SVector{8,Float64}[]
    Ap = SVector{8,Float64}[]

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
                JD_k  = date_to_jd(year, month, day, 12, 0, 0)

                # Get the vector of Kps and Aps.
                Ap_k = zeros(Float64,8)
                Kp_k = zeros(Float64,8)

                for i = 1:8
                    Kp_k[i] = parse(Int, ln[2(i-1) + 13:2(i-1) + 14])/10
                    Ap_k[i] = parse(Int, ln[3(i-1) + 32:3(i-1) + 34])
                end

                # Add data to the vector.
                push!(JD, JD_k)
                push!(Kp, SVector{8,Float64}(Kp_k))
                push!(Ap, SVector{8,Float64}(Ap_k))
            end
        end
    end

    # Create the interpolations for each parameter.
    knots    = (JD,)

    # Create the interpolations.
    itp_Kp = interpolate(knots, Kp, Gridded(Constant()))
    itp_Ap = interpolate(knots, Ap, Gridded(Constant()))

    return _WDC_Structure(itp_Kp, itp_Ap)
end

"""
    _prepare_wdc_remote_files(oldest_year::Number, newest_year::Number)

Configure all the WDC remote files between `newest_year` and `oldest_year`.
Notice that previous years will never be updated whereas the current year will
be updated daily.

If `oldest_year` is greater than current year, then only the files from the
current year will be downloaded.

If `newest_year` is smaller than `oldest_year`, then only the files from the
`oldest_year` will be downloaded.

This function modifies the global variable `_wdcfiles`.

"""
function _prepare_wdc_remote_files(oldest_year::Number, newest_year::Number)
    # Get the current year.
    current_year = year(now())

    # If `oldest_year` is greater than current year, then consider only the
    # current year.
    (oldest_year > current_year) && (oldest_year = current_year)
    (newest_year < oldest_year)  && (newest_year = oldest_year)
    (newest_year > current_year) && (newest_year = current_year)

    # For the current year, we must update the remote file every day. Otherwise,
    # we do not need to update at all.
    for y = oldest_year:newest_year
		filename = "kp$y"
        sym = Symbol(filename)
        file_y = @RemoteFile("ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/wdc/yearly/$filename.wdc",
                             file="$filename.wdc",
                             updates= (y == current_year) ? :daily : :never)

        merge!(_wdcfiles.files, Dict(sym => file_y))
    end

    return nothing
end
