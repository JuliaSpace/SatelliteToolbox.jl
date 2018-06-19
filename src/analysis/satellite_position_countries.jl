#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions to verify if the satellite is above some countries.
#
#   Initial version based on the code sent by Renato Branco.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export satellite_check_Brazil

"""
    function satellite_check_Brazil(lat::Number, lon::Number)

Verify if a point described by latitude `lat` and longitude `lon` is inside
Brazil.

# Args

* `lat`: Latitude [rad].
* `lon`: Longitude [rad].

# Returns

* `TRUE`: The point is inside Brazil.
* `FALSE`: The point is not inside Brazil.

# Remarks

This function was based on the algorithm sent by Renato Branco to Ronan Arraes
by e-mail at 2016-02-16.

"""
function satellite_check_Brazil(lat::Number, lon::Number)
    # Notes
    # ====
    #
    # SB: Superior bound.
    # IB: Inferior bound.

    # Convert latitude and longitude to degrees.
    lat *= 180.0/pi
    lon *= 180.0/pi

    # Rough verification
    # ==================
    if ( (lon <= -32.5) && (lon >= -75.0) && (lat <= 5.0) && (lat >= -35.0) )
        # Accurate verification
        # =====================

        # Zone 01.
        if (lon > -40.0)
            SB = -0.4*lon - 18.0
            IB =  2.0*lon + 60.0
        # Zone 02.
        elseif (lon > -45.0)
            SB = -0.4*lon   - 18.0
            IB =  0.667*lon + 6.667
        # Zone 03.
        elseif (lon > -47.5)
            SB = 0.0
            IB = 0.667*lon + 6.667
        # Zone 04.
        elseif (lon > -50.0)
            SB = 0.0
            IB = 2.0*lon + 70.0
        # Zone 05.
        elseif (lon > -52.5)
            SB = -2.0*lon - 100.0
            IB =  2.0*lon +  70.0
        # Zone 06.
        elseif (lon > -55.0)
            SB =  1.0*lon + 57.5
            IB = -1.0*lon - 87.5
        # Zone 07.
        elseif (lon > -57.5)
            SB = -0.5*lon - 25.0
            IB = -1.0*lon - 87.5
        # Zone 08.
        elseif (lon > -60.0)
            SB = -0.5*lon -  25.0
            IB = -5.0*lon - 317.5
        # Zone 09.
        elseif (lon > -67.5)
            SB =  0.25*lon + 20.0
            IB = -1.00*lon - 77.5
        # Zone 10.
        elseif (lon > -70.0)
            SB = 0.25*lon + 20.0
            IB = -10.0
        # Zone 11.
        elseif (lon > -72.5)
            SB = 1.5*lon + 107.5
            IB = -10.0
        # Zone 12.
        else
            SB =  1.5*lon + 107.5
            IB = -2.0*lon - 155.0
        end

        return ( (lat <= SB) && (lat >= IB) )
    end

    return false
end
