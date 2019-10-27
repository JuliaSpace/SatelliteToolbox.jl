#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Function related to TLE (Two line elements).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export read_tle, read_tle_from_string
export @tle_str, @tlenc_str

"""
### macro parse_value(T, str, line_num)

Parse the string `str` using the type `T`. If it is not succeeded, then throw an
error indicating the line `line_num` with the problem.

"""
macro parse_value(T, str, line_num)
    quote
        local _T = $(esc(T))
        local _str = $(esc(str))
        local _line_num = $(esc(line_num))

        try
            parse(_T, _str)
        catch e
            throw(ErrorException("The TLE file is not valid (error in line $_line_num): $(e.msg)."))
        end
    end
end

"""
    function compute_checksum(str::AbstractString)

Compute the checksum of the line `str` modulo 10.

The algorithm is simple: add all the numbers in the line, ignoring letters,
spaces, periods, and plus signs, but assigning +1 to the minus signs. The
checksum is the remainder of the division by 10.

"""
function compute_checksum(str::AbstractString)
    checksum = 0

    for c in str
        # Check if `c` is a number.
        if isnumeric(c)
            checksum += parse(Int, c)

        # Check if `c` is a minus sign, which has value 1.
        elseif c == '-'
            checksum += 1
        end

        # Otherwise, the value of the character is 0.
    end

    # Return the checksum modulo 10.
    checksum % 10
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
#                              Private Functions
################################################################################

function _parse_tle(io::IO, verify_checksum::Bool = true)
    # State machine to read the TLE. It has three possible states:
    #
    #   :name -> Satellite name.
    #   :l1   -> Line 1.
    #   :l2   -> Line 2.
    #
    # The transitions are:
    #
    #   :name --> :l1 --> :l2
    #     ^                |
    #     |                |
    #      ----------------
    #
    state = :name

    # Auxiliary variables that will store the TLE.
    name = ""
    sat_num = 0
    classification = Char(0)
    int_designator = ""
    epoch_year = 0
    epoch_day = 0.0
    epoch = 0.0
    dn_o2 = 0.0
    ddn_o6 = 0.0
    bstar = 0.0
    elem_set_number = 0
    checksum_l1 = 0
    i = 0.0
    Ω = 0.0
    e = 0.0
    ω = 0.0
    M = 0.0
    n = 0.0
    rev_num = 0
    checksum_l2 = 0

    # Output array with the TLEs found in the file.
    tles = Array{TLE}(undef,0)

    line_num = 0

    skip_line_read = false

    line = nothing

    while !eof(io)
        # Read the current line, strip white spaces, and skip if it is blank or
        # is a comment.
        if !skip_line_read
            line = strip(readline(io))
            line_num += 1
            (isempty(line))  && continue
            (line[1] == '#') && continue
        else
            skip_line_read = false
        end

        # Check the state of the reading.
        if state == :name
            # Check if the line seems to be the first line of the TLE. In this
            # case, maybe the user has not provided a name. Then, change the
            # state to `:l1`, and skip the line reading in the next loop.
            if (line[1:2] == "1 ") && (length(line) == 69)
                name = "UNDEFINED"
                state = :l1
                skip_line_read  = true
                continue
            end

            # Otherwise, if the line is not blank, then it must be the name of
            # the satellite.
            #
            # NOTE: The name should not be bigger than 24 characters. However,
            # we will not check this here.

            # Copy the name and change to state to wait for the 1st line.
            name  = line
            state = :l1

        # TLE Line 1
        # ==========
        elseif state == :l1
            # The next non-blank line must be the first line of the TLE.
            # Otherwise, the file is not valid.

            # The first line must start with "1 " and have 69 characters.
            if (line[1:2] != "1 ") || (length(line) != 69)
                throw(ErrorException("The TLE file is not valid (error in line $line_num): This is not a valid 1st line."))
            end

            # Verify the checksum.
            # ====================

            checksum_l1   = @parse_value(Int, line[69], line_num)
            checksum_line = compute_checksum(line[1:end-1])

            if (verify_checksum) && (checksum_l1 != checksum_line)
                throw(ErrorException("The TLE file is not valid (error in line $line_num): Expected checksum: $checksum_line, line checksum: $checksum_l1."))
            end

            # Satellite number and classification
            # ===================================
            sat_num        = @parse_value(Int, line[3:7], line_num)
            classification = Char(line[8])

            # International designator
            # ========================
            int_designator = line[10:17]

            # Epoch
            # =====
            epoch_year = @parse_value(Int,     line[19:20], line_num)
            epoch_day  = @parse_value(Float64, line[21:32], line_num)

            # Convert the TLE year and date to JD.
            # ------------------------------------

            # If the two-number year is higher than the current one, then
            # consider it in the past (e.g. if today is 2018, then 18 = 2018 but
            # 19 = 1919.
            aux = Dates.year(now()) - 2000

            if epoch_year > aux
                epoch = Dates.datetime2julian(DateTime(1900 + epoch_year, 1, 1, 0, 0, 0)) - 1 +
                        epoch_day
            else
                epoch = Dates.datetime2julian(DateTime(2000 + epoch_year, 1, 1, 0, 0, 0)) - 1 +
                        epoch_day
            end

            # Mean motion derivatives
            # =======================
            dn_o2  = @parse_value(Float64, line[34:43], line_num)

            aux = ( (line[45] == ' ') ? "+." : line[45:45] * "." )*line[46:50]

            ddn_o6_dec = @parse_value(Float64, aux, line_num)
            ddn_o6_exp = @parse_value(Float64, line[51:52], line_num)
            ddn_o6 = ddn_o6_dec*10^ddn_o6_exp

            # BSTAR
            # =====

            aux = ( (line[54] == ' ') ? "+." : line[54:54] * "." )*line[55:59]

            bstar_dec = @parse_value(Float64, aux, line_num)
            bstar_exp = @parse_value(Float64, line[60:61], line_num)
            bstar = bstar_dec*10^bstar_exp

            # Ephemeris type
            # ==============
            (line[63] != '0' && line[63] != ' ') &&
                warn("Warning in TLE file (line $line_num): Ephemeris type should be 0!")

            # Element number
            # ==============
            elem_set_number = @parse_value(Int, line[65:68], line_num)

            # Change the state to wait for the line 2.
            state = :l2

        # TLE Line 2
        # ==========
        elseif state == :l2
            # The next non-blank line must be the second line of the TLE.
            # Otherwise, the file is not valid.

            # The second line must start with "2 " and have 69 characters.
            if (line[1:2] != "2 ") || (length(line) != 69)
                throw(ErrorException("The TLE file is not valid (error in line $line_num): This is not a valid 2st line."))
            end

            # Verify the checksum.
            # ====================

            checksum_l2   = @parse_value(Int, line[69], line_num)
            checksum_line = compute_checksum(line[1:end-1])

            if (verify_checksum) && (checksum_l2 != checksum_line)
                throw(ErrorException("The TLE file is not valid (error in line $line_num): Expected checksum: $checksum_line, line checksum: $checksum_l2."))
            end

            # Check satellite number with the one in the first line.
            # ======================================================
            sat_num_l2 = @parse_value(Int, line[3:7], line_num)

            if sat_num_l2 != sat_num
                throw(ErrorException("The TLE file is not valid (error in line $line_num): Satellite number in line 2 is not equal to that in line 1."))
            end

            # Inclination
            # ===========
            i = @parse_value(Float64, line[9:16], line_num)

            # RAAN
            # ====
            Ω = @parse_value(Float64, line[18:25], line_num)

            # Eccentricity
            # ============
            e = @parse_value(Float64, "." * line[27:33], line_num)

            # Argument of perigee
            # ===================
            ω = @parse_value(Float64, line[35:42], line_num)

            # Mean anomaly
            # ============
            M = @parse_value(Float64, line[44:51], line_num)

            # Mean motion
            # ===========
            n = @parse_value(Float64, line[53:63], line_num)

            # Revolution number at epoch
            # ==========================
            rev_num = @parse_value(Int, line[64:68], line_num)

            # Now that we have all the information, we can store our TLE in the
            # output array.
            push!(tles,
                  TLE(name,
                      sat_num,
                      classification,
                      int_designator,
                      epoch_year,
                      epoch_day,
                      epoch,
                      dn_o2,
                      ddn_o6,
                      bstar,
                      elem_set_number,
                      checksum_l1,
                      i,
                      Ω,
                      e,
                      ω,
                      M,
                      n,
                      rev_num,
                      checksum_l2))

            # Change the state to wait for another TLE.
            state = :name
        end
    end

    # If the final state is not :name, then we have an incomplete TLE. Thus,
    # throw and exception because the file is not valid.
    if state != :name
        throw(ErrorException("The TLE file is not valid.\nThere is an incomplete TLE."))
    end

    tles
end

"""
    function print_tle(io::IO, tle::TLE, color::Bool = true)

Print the TLE `tle` in the IO `io`.

If `color` is `true`, then the text will be printed using colors. If `color` is
omitted, then it defaults to `true`.

"""
function print_tle(io::IO, tle::TLE, color::Bool = true)
    # Colors will be printed only for STDOUT.
    b = (color) ? _b : ""
    d = (color) ? _d : ""
    g = (color) ? _g : ""
    y = (color) ? _y : ""

    # Print the TLE information.
    print(io, "                             $(g)TLE$(d)\n")
    print(io, "$(y)    ==========================================================$(d)\n")
    print(io, "$(b)                            Name: $(d)"); @printf(io, "%s\n", tle.name)
    print(io, "$(b)                Satellite number: $(d)"); @printf(io, "%d\n", tle.sat_num)
    print(io, "$(b)        International designator: $(d)"); @printf(io, "%s\n", tle.int_designator)
    print(io, "$(b)                    Epoch (Year): $(d)"); @printf(io, "%d\n", tle.epoch_year)
    print(io, "$(b)                     Epoch (Day): $(d)"); @printf(io, "%12.8f\n", tle.epoch_day)
    print(io, "$(b)              Epoch (Julian Day): $(d)"); @printf(io, "%12.5f\n", tle.epoch)
    print(io, "$(b)              Element set number: $(d)"); @printf(io, "%d\n", tle.elem_set_number)
    print(io, "$(b)                    Eccentricity: $(d)"); @printf(io, "%12.8f deg\n", tle.e)
    print(io, "$(b)                     Inclination: $(d)"); @printf(io, "%12.8f deg\n", tle.i)
    print(io, "$(b)                            RAAN: $(d)"); @printf(io, "%12.8f deg\n", tle.Ω)
    print(io, "$(b)             Argument of perigee: $(d)"); @printf(io, "%12.8f deg\n", tle.ω)
    print(io, "$(b)                    Mean anomaly: $(d)"); @printf(io, "%12.8f deg\n", tle.M)
    print(io, "$(b)                 Mean motion (n): $(d)"); @printf(io, "%12.8f revs/day\n", tle.n)
    print(io, "$(b)               Revolution number: $(d)"); @printf(io, "%d\n", tle.rev_num)
    print(io, "\n")
    print(io, "$(b)                              B*: $(d)"); @printf(io, "%f 1/[er]\n", tle.bstar)
    print(io, "\n")
    print(io, "$(b)                        1   d\n$(d)")
    print(io, "$(b)                       ---.--- n: $(d)"); @printf(io, "%f rev/day²\n", tle.dn_o2)
    print(io, "$(b)                        2  dt\n$(d)")
    print(io, "\n")
    print(io, "$(b)                        1   d²\n$(d)")
    print(io, "$(b)                       ---.--- n: $(d)"); @printf(io, "%f rev/day³\n", tle.ddn_o6)
    print(io, "$(b)                        6  dt²\n$(d)")
    print(io, "$(y)    ==========================================================$(d)")

    nothing
end

function show(io::IO, tle::TLE)
    print(io, "TLE: ", tle.name, " (Epoch = ", JDtoDate(DateTime,tle.epoch), ")")
    return nothing
end

function show(io::IO, mime::MIME"text/plain", tle::TLE)
    color = get(io, :color, false)
    print_tle(io, tle, color)
end
