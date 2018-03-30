#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Function related to TLE (Two line elements).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-03-28: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export read_tle

"""
### macro parse_value(T, str, line_num)

Parse the string `str` using the type `T`. If it is not succeeded, then throw an
error indicating the line `line_num` with the problem.

##### Args

* T: Type of the output.
* str: Input string.
* line_num: Line number that is being analyzed.

##### Returns

The `str` parsed to the type `T`. If an error occurred, then an exception is
thrown.

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
### function compute_checksum(str::String)

Compute the checksum of the line `str` modulo 10. The algorithm is simple: add
all the numbers in the line, ignoring letters, spaces, periods, and plus signs,
but assigning +1 to the minus signs. The checksum is the remainder of the
division by 10.

##### Args

* str: String to be checked, the checksum **must not** be in this string.

##### Returns

The computed checksum.

"""

function compute_checksum(str::String)
    checksum = 0

    for c in str
        # Check if `c` is a number.
        if isnumber(c)
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
### function read_tle(tle_filename::String)

Read the TLEs in the file `tle_filename`.

##### Args

* tle_filename: TLE file name.
* verify_checksum: (OPTIONAL) If false, then the checksum will not be verified
  (**DEFAULT** = true).

##### Returns

An array with all the TLEs that were parsed.

"""

function read_tle(tle_filename::String, verify_checksum::Bool = true)
    # Open the file in read mode.
    file = open(tle_filename, "r")

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
    tles = Array{TLE}(0)

    line_num = 0

    while !eof(file)
        # Read the current line, strip white spaces, and skip if it is blank.
        line = strip(readline(file))
        line_num += 1
        (isempty(line)) && continue

        # Check the state of the reading.
        if state == :name
            # If the line is not blank, then it must be the name of the
            # satellite.
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
            (line[63] != '0') && warn("Ephemeris type should be 0!")

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
