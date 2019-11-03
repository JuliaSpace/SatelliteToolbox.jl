# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Main file for TLE.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export @tle_str, @tlenc_str
export read_tle, read_tle_from_string

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
    print_tle(io, tle, color)
end

################################################################################
#                                  Functions
################################################################################

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
