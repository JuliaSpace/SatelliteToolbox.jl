# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions related with the IERS (International Earth Orientation and
#   Reference Systems Service) EOP (Earth Orientation Parameters) data.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] IERS (2010). Transformation between the International Terrestrial
#       Reference System and the Geocentric Celestial Reference System. IERS
#       Technical Note No. 36, Chapter 5.
#
#   [3] ftp://hpiers.obspm.fr/eop-pc/models/uai2000.package
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export EOPData_IAU1980, EOPData_IAU2000A
export get_iers_eop, get_iers_eop_iau_1980, get_iers_eop_iau_2000A
export read_iers_eop
export deps_dpsi

################################################################################
#                       Private Structures and Variables
################################################################################

# Interpolation types used in EOP data.
const _eop_itp_linear{T} =
    Interpolations.GriddedInterpolation{T,1,T, Gridded{Linear}, Tuple{Array{T,1}}}

const _eop_etp_linear{T} =
    Interpolations.Extrapolation{T, 1, _eop_itp_linear{T}, Gridded{Linear}, Flat{Nothing}}

################################################################################
#                                  Functions
################################################################################

"""
    get_iers_eop(data_type::Symbol = :IAU1980; force_download = false)

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
    get_iers_eop_iau_1980(url::String = "https://datacenter.iers.org/data/latestVersion/223_EOP_C04_14.62-NOW.IAU1980223.txt")

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
    get_iers_eop_iau_2000A(url::String = "https://datacenter.iers.org/data/latestVersion/224_EOP_C04_14.62-NOW.IAU2000A224.txt"; force_download = false)

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
    read_iers_eop(filename::String, data_type::Symbol = :IAU1980)

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

    EOPData_IAU1980{_eop_etp_linear{Float64}}(
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

    EOPData_IAU2000A{_eop_etp_linear{Float64}}(
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

################################################################################
#                                   Helpers
################################################################################

"""
    deps_dpsi(eop_iau2000a::EOPData_IAU2000A, JD::Number)

Returns the celestial pole offsets in obliquity (δϵ_2000) and longitude
(δΨ_2000) [arcsec]. This function obtains those values by converting the
celestial pole offsets with respect to the GCRS (`dX` and `dY`). These values
are necessary in the equinox-based IAU-2006 theory.

The algorithm was obtained from [2, eq. 5.25] and [3, `DPSIDEPS2000_DXDY2000`].

"""
function deps_dpsi(eop_iau2000a::EOPData_IAU2000A, JD::Number)
    # Constants.
    d2r = π / 180
    a2d = 1 / 3600
    a2r = a2d * d2r

    # Obtain the parameters `dX` and `dY` [arcseg].
    dX = eop_iau2000a.dX(JD)
    dY = eop_iau2000a.dY(JD)

    # Compute the Julian century.
    T_TT = (JD - JD_J2000) / 36525

    # Luni-solar precession [rad].
    Ψ_a = @evalpoly(T_TT, 0, +5038.47875, -1.07259, -0.001147) * a2r

    # Planetary precession [rad].
    χ_a = @evalpoly(T_TT, 0, +10.5526, -2.38064, -0.001125) * a2r

    sϵ₀, cϵ₀ = sincos(84381.406 * a2r)

    aux = Ψ_a * cϵ₀ - χ_a
    den = aux^2 * sϵ₀ - sϵ₀
    δϵ  = (aux * sϵ₀ * dX - sϵ₀ * dY)/den
    δΨ  = (dX - aux * dY)/den

    return δϵ, δΨ
end
