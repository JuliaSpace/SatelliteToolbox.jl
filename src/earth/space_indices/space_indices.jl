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

include("./dtcfile.jl")
include("./solfsmy.jl")
include("./wdcfiles.jl")

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
* `wdcfiles_dir`: Path to the directory containing the WDC files.

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
        download(_dtcfile; force = force_download)
        dtcfile_path = path(_dtcfile)
    end

    if solfsmy_path == nothing
        download(_solfsmy; force = force_download)
        solfsmy_path = path(_solfsmy)
    end

    years     = Int[]
    filepaths = String[]

    if wdcfiles_dir == nothing
        _prepare_wdc_remote_files(wdcfiles_oldest_year)
        download(_wdcfiles; force = force_download)

        # Get the files available and sort them by the year.
        for (sym,wdcfile) in _wdcfiles.files
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
    push!(_dtcfile_data, _parse_dtcfile(dtcfile_path))
    push!(_solfsmy_data, _parse_solfsmy(solfsmy_path))
    push!(_wdc_data,     _parse_wdcfiles(filepaths[p], years[p]))

    nothing
end
