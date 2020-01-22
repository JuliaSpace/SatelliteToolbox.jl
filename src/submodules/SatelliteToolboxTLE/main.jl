# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Main file for TLE.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export @tle_str, @tlenc_str
export print_tle, read_tle, read_tle_from_string, tle_to_str

################################################################################
#                                    Macros
################################################################################

"""
    macro tle_str(str)

Parse a set of TLEs in the string `str` and return as an array of `TLE`. This
version verifies the checksum of the TLE. If the checksum verification is not
desired, see `@tlenc_str`.

# Example

```julia-repl
julia> tles = tle\"""
       CBERS 4
       1 40336U 14079A   18166.15595376 -.00000014  00000-0  10174-4 0  9993
       2 40336  98.4141 237.7928 0001694  75.7582 284.3804 14.35485112184485
       SCD 1
       1 22490U 93009B   18165.62596833  .00000225  00000-0  11410-4 0  9991
       2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338066
       SCD 2
       1 25504U 98060A   18165.15074951  .00000201  00000-0  55356-5 0  9994
       2 25504  24.9961  80.1303 0017060 224.4822 286.6438 14.44043397 37312
       \"""
```
"""
macro tle_str(str)
    read_tle_from_string(str, true)
end

"""
    macro tlenc_str(str)

Parse a set of TLEs in the string `str` and return as an array of `TLE`. This
version **does not** verify the checksum of the TLE. If the checksum
verification is required, see `@tle_str`.

# Example

```julia-repl
julia> tles = tlenc\"""
       CBERS 4
       1 40336U 14079A   18166.15595376 -.00000014  00000-0  10174-4 0  9993
       2 40336  98.4141 237.7928 0001694  75.7582 284.3804 14.35485112184485
       SCD 1
       1 22490U 93009B   18165.62596833  .00000225  00000-0  11410-4 0  9991
       2 22490  24.9690 231.7852 0042844 200.7311 292.7198 14.44524498338066
       SCD 2
       1 25504U 98060A   18165.15074951  .00000201  00000-0  55356-5 0  9994
       2 25504  24.9961  80.1303 0017060 224.4822 286.6438 14.44043397 37312
       \"""
```
"""
macro tlenc_str(str)
    read_tle_from_string(str, false)
end

################################################################################
#                                  Overloads
################################################################################

function show(io::IO, tle::TLE)
    print(io, "TLE: ", tle.name, " (Epoch = ", julian2datetime(tle.epoch), ")")
    return nothing
end

function show(io::IO, mime::MIME"text/plain", tle::TLE)
    color = get(io, :color, false)
    _show_tle(io, tle, color)
end

################################################################################
#                                  Functions
################################################################################

"""
    function print_tle(io::IO, tle; kwargs...)

Print the TLE `tle` to the IO `io`. If `io` is omited, then `stdout` is used.

The keywords of this function are the same that can be used in `tle_to_str`.

"""
print_tle(tle; kwargs...) = print_tle(stdout, tle; kwargs...)

function print_tle(io::IO, tle; kwargs...)
    str = tle_to_str(tle; kwargs...)
    print(io, str)

    return nothing
end

"""
    function read_tle(tle_filename::String, verify_checksum::Bool = true)

Read the TLEs in the file `tle_filename` and return an array of `TLE` with the
parsed TLEs.

If `verify_checksum` if `true`, then the checksum of both TLE lines will be
verified. Otherwise, the checksum will not be checked. If `verify_checksum` is
omitted, then it defaults to `true`.

"""
@inline function read_tle(tle_filename::String, verify_checksum::Bool = true)
    # Open the file in read mode.
    file = open(tle_filename, "r")

    _parse_tle(file, verify_checksum)
end

"""
    function read_tle_from_string(tles::String, verify_checksum::Bool = true)
    function read_tle_from_string(tle_l1::String, tle_l2::String, verify_checksum::Bool = false)

Parse a set of TLEs in the string `tles` or one TLE with first line `tle_l1` and
second line `tle_l2`. This function returns an array of `TLE` with the parsed
TLEs.

If `verify_checksum` if `true`, then the checksum of both TLE lines will be
verified. Otherwise, the checksum will not be checked. If `verify_checksum` is
omitted, then it defaults to `true`.

"""
@inline function read_tle_from_string(tles::String,
                                      verify_checksum::Bool = true)
    # Convert the string to an IOBuffer and call the function to parse it.
    _parse_tle(IOBuffer(tles), verify_checksum)
end

@inline function read_tle_from_string(tle_l1::String,
                                      tle_l2::String,
                                      verify_checksum::Bool = false)
    # Assemble the TLE into one string and call the function to parse it.
    read_tle_from_string(tle_l1 * '\n' * tle_l2, verify_checksum)
end


"""
    function tle_to_str(tle::TLE; recompute_checksum = true, bstar_exp_le = true)

Convert the TLE `tle` to a string. If `recompute_checksum` is `true`, then the
checksums in `tle` will be ignored and they will be computed considering the TLE
data.

The keyword `bstar_exp_le` selects if the BSTAR exponent signal will be `+` or
`-` when BSTAR is zero. This is required for the tests because it is not
standardized in TLE generation. If it is `true`, then the exponent signal will
be `-` when BSTAR is zero.

"""
function tle_to_str(tle::TLE; recompute_checksum = true, bstar_exp_le = true)
    @unpack_TLE tle

    # Name
    # ==========================================================================

    str_name = @sprintf("%-24s",name)

    # First line
    # ==========================================================================

    # Line ID
    str_l1  = "1 "

    # Satellite number and classification
    str_sat_num = @sprintf("%05d", sat_num)[1:5]
    str_l1 *= str_sat_num
    str_l1 *= classification * " "

    # International designator
    if length(int_designator) > 8
        str_l1 *= int_designator[1:8]
    else
        str_l1 *= int_designator * " "^(8-length(int_designator))
    end
    str_l1 *= " "

    # Epoch (year)
    str_l1 *= @sprintf("%-2d", epoch_year)[1:2]

    # Epoch (day + fraction of the day)
    i_epoch_day = floor(Int,epoch_day)
    f_epoch_day = epoch_day - i_epoch_day
    str_l1 *= @sprintf("%03d", i_epoch_day)[1:3] * "." * @sprintf("%-10.8f", f_epoch_day)[3:10]
    str_l1 *= " "

    # First time derivative of the mean motion
    str_l1 *= dn_o2 < 0 ? "-." : " ."
    str_l1 *= @sprintf("%-10.8f", abs(dn_o2))[3:10]
    str_l1 *= " "

    # Second time derivative of the mean motion
    mant, exp = _get_mant_exp(abs(ddn_o6))
    str_l1 *= ddn_o6 < 0 ? "-" : " "
    str_l1 *= @sprintf("%-7.5f", mant)[3:7]
    str_l1 *= exp ≤ 0 ? "-" : "+"
    str_l1 *= @sprintf("%-2d", abs(exp))[1]
    str_l1 *= " "

    # BSTAR drag term
    mant, exp = _get_mant_exp(abs(bstar))
    str_l1 *= bstar < 0 ? "-" : " "
    str_l1 *= @sprintf("%-7.5f", mant)[3:7]
    if bstar_exp_le
        str_l1 *= exp ≤ 0 ? "-" : "+"
    else
        str_l1 *= exp < 0 ? "-" : "+"
    end
    str_l1 *= @sprintf("%-2d", abs(exp))[1]
    str_l1 *= " "

    # Ephemeris type
    str_l1 *= "0 "

    # Element number
    str_l1 *= @sprintf("%4d", elem_set_number)[1:4]

    # Checksum
    recompute_checksum && (checksum_l1 = compute_checksum(str_l1))
    str_l1 *= @sprintf("%d", checksum_l1)[1]

    # Second line
    # ==========================================================================

    # ID
    str_l2 = "2 "

    # Satellite number
    str_l2 *= str_sat_num
    str_l2 *= " "

    # Inclination [°]
    str_l2 *= @sprintf("%8.4f", i)
    str_l2 *= " "

    # Right ascension of the ascending node [°]
    str_l2 *= @sprintf("%8.4f", Ω)
    str_l2 *= " "

    # Eccentricity [°]
    str_l2 *= @sprintf("%-9.7f", e)[3:9]
    str_l2 *= " "

    # Argument of perigee [°]
    str_l2 *= @sprintf("%8.4f", ω)
    str_l2 *= " "

    # Mean anomaly [°]
    str_l2 *= @sprintf("%8.4f", M)
    str_l2 *= " "

    # Mean motion [rev/day]
    str_l2 *= @sprintf("%11.8f", n)

    # Revolution number at epoch [revs]
    str_l2 *= @sprintf("%5d", rev_num)[1:5]

    # Checksum
    recompute_checksum && (checksum_l2 = compute_checksum(str_l2))
    str_l2 *= @sprintf("%d", checksum_l2)[1]

    # Assemble and return the TLE.
    # ==========================================================================
    return str_name * '\n' * str_l1 * '\n' * str_l2
end
