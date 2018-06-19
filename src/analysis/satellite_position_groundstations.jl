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
minimum elevation angle of `minElev`.

# Args

* `r_e`: Satellite position represented in the ECEF reference frame.
* `rs_e`: Ground station position represented in the ECEF reference frame.
* `minElev`: Minimum elevation angle in which the station can see the satellite
             [rad].

# Returns

* `TRUE`: The satellite is inside the visibility circle of the ground station.
* `FALSE`: The satellite is not inside the visibility circle of the ground
           station.

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
visibility circle of a ground station with latitude `lat_s`, longitude `lon_s`,
altitude `h_s`, and a minimum elevation angle of `minElev`.

# Args

* `r_e`: Satellite position represented in the ECEF reference frame.
* `lat_s`: Latitude of the ground station [rad].
* `lon_s`: Longitude of the ground station [rad].
* `h_s`: Altitude of the ground station (WGS-84) [m].
* `minElev`: Minimum elevation angle in which the station can see the satellite
           [rad].

# Returns

* `TRUE`: The satellite is inside the visibility circle of the ground station.
* `FALSE`: The satellite is not inside the visibility circle of the ground
           station.

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
