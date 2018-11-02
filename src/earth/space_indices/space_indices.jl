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
include("./fluxtable.jl")
include("./solfsmy.jl")
include("./wdcfiles.jl")

################################################################################
#                               Public Functions
################################################################################

"""
    function init_space_indices(...)

Initialize all space indices. The files that will be initialized must be
indicated by the array of symbols passed to the keyword argument
`enabled_files`. If this is `nothing`, which is the default, then all files will
be initialized. The symbol related to each file is described next.

Notice that the initialization process can be changed by a set of keywords as
described next.

## DTCFILE

**Symbol**: `:dtcfile`

This file contains the exospheric temperature variation caused by the Dst index.
This is used for the JB2008 atmospheric model.

### Keywords

* `dtcfile_path`: Path for the file `DTCFILE.TXT`. If `nothing`, then it will be
                  downloaded. (**Default** = `nothing`)
* `dtcfile_force_download`: If `true`, then the file will always be downloaded
                            if the path is not specified. (**Default** =
                            `false`).

## fluxtable

**Symbol**: `:fluxtable`

This file contains the F10.7 flux data in different formats.

### Keywords

* `fluxtable_path`: Path for the file `fluxtable.txt`. If `nothing`, then it
                    will be downloaded. (**Default** = `nothing`)
* `fluxtable_force_download`: If `true`, then the file will always be downloaded
                              if the path is not specified.
                              (**Default** = `false`).

## SOLFSMY

**Symbol**: `:solfsmy`

This files contains the indices necessary for the JB2008 atmospheric model.

### Keywords

* `solfsmy_path`: Path for the file `SOLFSMY.TXT`. If `nothing`, then it will be
                  downloaded. (**Default** = `nothing`)
* `solfsmy_force_download`: If `true`, then the file will always be downloaded
                            if the path is not specified. (**Default** =
                            `false`).

## WDC Files

**Symbol**: `:wdcfiles`

This set of files contain the Kp and Ap indices.

### Keywords

* `wdcfiles_path`: Path for the directory with the WDC files. If `nothing`, then
                   they will be downloaded. (**Default** = `nothing`)
* `wdcfiles_force_download`: If `true`, then the files will always be downloaded
                            if the path is not specified. (**Default** =
                            `false`).
* `wdcfiles_oldest_year`: Oldest year in which the WDC file will be obtained.
                          (**Default** = past 3 years).

"""
function init_space_indices(;enabled_files = nothing,
                             dtcfile_path = nothing,
                             dtcfile_force_download = false,
                             fluxtable_path = nothing,
                             fluxtable_force_download = false,
                             solfsmy_path = nothing,
                             solfsmy_force_download = false,
                             wdcfiles_dir = nothing,
                             wdcfiles_force_download = false,
                             wdcfiles_oldest_year = year(now())-3)

    dtcfile   = (enabled_files == nothing) || (:dtcfile in enabled_files)
    fluxtable = (enabled_files == nothing) || (:fluxtable in enabled_files)
    solfsmy   = (enabled_files == nothing) || (:solfsmy in enabled_files)
    wdcfiles  = (enabled_files == nothing) || (:wdcfiles in enabled_files)

    dtcfile && _init_dtcfile(local_path = dtcfile_path,
                             force_download = dtcfile_force_download)

    fluxtable && _init_fluxtable(local_path = fluxtable_path,
                                 force_download = fluxtable_force_download)

    solfsmy && _init_solfsmy(local_path = solfsmy_path,
                             force_download = solfsmy_force_download)

    wdcfiles && _init_wdcfiles(local_dir = wdcfiles_dir,
                               force_download = wdcfiles_force_download,
                               wdcfiles_oldest_year = wdcfiles_oldest_year)

    nothing
end
