#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions to verify if the satellite is within the station visibility
#   circle.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export satellite_check_station

"""
    function satellite_check_station(r_e::Vector, rs_e::Vector, minElev::Number)

Check if the satellite with position vector `r_e` (ECEF) is inside the
visibility circle of a ground station with position vector `rs_e` (ECEF) and a
minimum elevation angle of `minElev` [rad].

Notice that `r_e` and `rs_e` must be represented in the same ECEF frame, and
must have the same unit.

Returns `true` if the satellite is inside the visibility circle, or `false`
otherwise.

"""
function satellite_check_station(r_e::Vector, rs_e::Vector, minElev::Number)
    # Check if the satellite is within the visibility circle of the station.
    dr_e = r_e - rs_e
    cos_beta = ( ( dr_e/norm(dr_e) )'*( rs_e/norm(rs_e) ) )[1]

    if (cos_beta > cos(pi/2-minElev))
        return true
    else
        return false
    end
end


"""
    function satellite_check_station(r_e::Vector, lat_s::Number, lon_s::Number, h_s::Number, minElev::Number)

Check if the satellite with position vector `r_e` (ECEF) is inside the
visibility circle of a ground station with latitude `lat_s` [rad], longitude
`lon_s` [rad], altitude `h_s` (WGS-84), and a minimum elevation angle of
`minElev` [rad].

Notice that the units of `r_e` and `h_s` must be the same.

Returns `true` if the satellite is inside the visibility circle, or `false`
otherwise.

"""
function satellite_check_station(r_e::Vector,
                                 lat_s::Number,
                                 lon_s::Number,
                                 h_s::Number,
                                 minElev::Number)
    # Convert the ground station LLA to the ECEF frame.
    rs_e = GeodetictoECEF(lat_s, lon_s, h_s)

    satellite_check_station(r_e, rs_e, minElev)
end
