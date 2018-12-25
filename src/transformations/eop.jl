#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions related with the IERS (International Earth Orientation and
#   Reference Systems Service) EOP (Earth Orientation Parameters) data.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export EOPData_IAU1980, EOPData_IAU2000A
export get_iers_eop, get_iers_eop_iau_1980, get_iers_eop_iau_2000A
export read_iers_eop

################################################################################
#                       Private Structures and Variables
################################################################################

"""
    function get_iers_eop(data_type::Symbol = :IAU1980; force_download = false)

Download and parse the IERS EOP C04 data. The data type is specified by
`data_type` symbol. Supported values are:

* `IAU1980`: Get IERS EOP C04 IAU1980 data.
* `IAU2000A`: Get IERS EOP C04 IAU2000A data.

If `data_type` is omitted, then it defaults to `IAU1980`.

The files are downloaded using the `RemoteFile` package with daily updates.
Hence, if one desires to force a download before the scheduled time, then set
the keyword `force_download` to `true`.

# Returns

A structure (`EOPData_IAU1980` or `EOPData_IAU2000A`, depending on `data_type`)
with the interpolations of the EOP parameters. Notice that the interpolation
indexing is set to the Julian Day.

"""
function get_iers_eop(data_type::Symbol = :IAU1980; force_download = false)
    if data_type == :IAU1980
        return get_iers_eop_iau_1980(force_download = force_download)
    elseif data_type == :IAU2000A
        return get_iers_eop_iau_2000A(force_download = force_download)
    else
        throw(ArgumentError("Unknow EOP type. It must be :IAU1980 or :IAU2000A."))
    end
end

"""
    function get_iers_eop_iau_1980(url::String = "https://datacenter.iers.org/data/latestVersion/223_EOP_C04_14.62-NOW.IAU1980223.txt")

Get the IERS EOP C04 IAU1980 data from the URL `url`. If `url` is omitted, then
it defaults to https://datacenter.iers.org/data/latestVersion/223_EOP_C04_14.62-NOW.IAU1980223.txt

The file is downloaded using the `RemoteFile` package with daily updates. Hence,
if one desires to force a download before the scheduled time, then set the
keyword `force_download` to `true`.

# Returns

The structure `EOPData_IAU1980` with the interpolations of the EOP parameters.
Notice that the interpolation indexing is set to the Julian Day.

# Remarks

For every field in `EOPData_IAU1980` to interpolation between two points in the
grid is linear. If extrapolation is needed, then if will use the nearest value
(flat extrapolation).

"""
function get_iers_eop_iau_1980(url::String = "https://datacenter.iers.org/data/latestVersion/223_EOP_C04_14.62-NOW.IAU1980223.txt";
                               force_download = false)

    _eop_iau1980 = @RemoteFile(@eval($url),
                               file="EOP_IAU1980.TXT",
                               updates=:daily)
    download(_eop_iau1980; force = force_download)

    # Parse the data removing the header.
    eop = readdlm(path(_eop_iau1980); skipblanks=true, skipstart=14)

    # Return the parsed EOP data.
    parse_iers_eop_iau_1980(eop)
end


"""
    function get_iers_eop_iau_2000A(url::String = "https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt"; force_download = false)

Get the IERS EOP C04 IAU2000A data from the URL `url`. If `url` is omitted, then
it defaults to https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt

The file is downloaded using the `RemoteFile` package with daily updates. Hence,
if one desires to force a download before the scheduled time, then set the
keyword `force_download` to `true`.

# Returns

The structure `EOPData_IAU2000A` with the interpolations of the EOP parameters.
Notice that the interpolation indexing is set to the Julian Day.

# Remarks

For every field in `EOPData_IAU2000A` to interpolation between two points in the
grid is linear. If extrapolation is needed, then if will use the nearest value
(flat extrapolation).

"""
function get_iers_eop_iau_2000A(url::String = "https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt";
                                force_download = false)

    _eop_iau2000A = @RemoteFile(@eval($url),
                                file="EOP_IAU2000A.TXT",
                                updates=:daily)
    download(_eop_iau2000A; force = force_download)

    # Parse the data removing the header.
    eop = readdlm(path(_eop_iau2000A); skipblanks=true, skipstart=14)

    # Return the parsed EOP data.
    parse_iers_eop_iau_2000A(eop)
end

"""
    function read_iers_eop(filename::String, data_type::Symbol = :IAU1980)

Read IERS EOP Data from the file `filename`. The user must specify if the data
is related to the model IAU 1980 (`data_type = :IAU1980`), which is the default,
or to the model IAU 2000A (`data_type = :IAU2000A`).

# Returns

A structure (`EOPData_IAU1980` or `EOPData_IAU2000A`, depending on `data_type`)
with the interpolations of the EOP parameters. Notice that the interpolation
indexing is set to the Julian Day.

# Remarks

The input file **must be exactly the same** as provided by IERS. One can
download it using the following commands:

* IAU 1980

        curl -O https://datacenter.iers.org/data/latestVersion/223_EOP_C04_14.62-NOW.IAU1980223.txt
        wget https://datacenter.iers.org/data/latestVersion/223_EOP_C04_14.62-NOW.IAU1980223.txt

* IAU 2000A

        curl -O https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt
        wget https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt

"""
function read_iers_eop(filename::String, data_type::Symbol = :IAU1980)
    # Parse the data removing the header.
    eop = readdlm(filename; skipblanks=true, skipstart=14)

    if data_type == :IAU1980
        parse_iers_eop_iau_1980(eop)
    elseif data_type == :IAU2000A
        parse_iers_eop_iau_2000A(eop)
    else
        throw(ArgumentError("Unknow EOP type. It must be :IAU1980 or :IAU2000A."))
    end
end

################################################################################
#                              Private Functions
################################################################################

function parse_iers_eop_iau_1980(eop::Matrix)
    # Check if the dimension seems correct.
    if size(eop,2) != 16
        error("The input data does not have the correct dimension.")
    end

    # Create the EOP Data structure by creating the interpolations.
    #
    # The interpolation will be linear between two points in the grid. The
    # extrapolation will be flat, considering the nearest point.
	knots = (eop[:,4] .+ 2400000.5,)

    EOPData_IAU1980(
        extrapolate(interpolate(knots, eop[:, 5], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:, 6], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:, 7], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:, 8], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:, 9], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,10], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,11], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,12], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,13], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,14], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,15], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,16], Gridded(Linear())), Flat())
    )
end

function parse_iers_eop_iau_2000A(eop::Matrix)
    # Check if the dimension seems correct.
    if size(eop,2) != 16
        error("The input data does not have the correct dimension.")
    end

    # Create the EOP Data structure by creating the interpolations.
    #
    # The interpolation will be linear between two points in the grid. The
    # extrapolation will be flat, considering the nearest point.
	knots = (eop[:,4] .+ 2400000.5,)

    EOPData_IAU2000A(
        extrapolate(interpolate(knots, eop[:, 5], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:, 6], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:, 7], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:, 8], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:, 9], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,10], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,11], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,12], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,13], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,14], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,15], Gridded(Linear())), Flat()),
        extrapolate(interpolate(knots, eop[:,16], Gridded(Linear())), Flat())
    )
end
