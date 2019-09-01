#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions to verify if the satellite is within the station visibility
#   circle.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export ground_station_visible

"""
    function ground_station_visible(r_e::AbstractVector, rs_e::AbstractVector, θ::Number)

Check if the satellite with position vector `r_e` (ECEF) is inside the
visibility circle of a ground station with position vector `rs_e` (ECEF) and a
minimum elevation angle of `θ` [rad].

Notice that `r_e` and `rs_e` must be represented in the same ECEF frame, and
must have the same unit.

Returns `true` if the satellite is inside the visibility circle, or `false`
otherwise.

"""
function ground_station_visible(r_e::AbstractVector, rs_e::AbstractVector,
                                θ::Number)
    # Check if the satellite is within the visibility circle of the station.
    dr_e = r_e - rs_e
    cos_beta = dot( dr_e/norm(dr_e), rs_e/norm(rs_e) )

    return cos_beta > cos(π/2-θ)
end


"""
    function ground_station_visible(r_e::AbstractVector, lat_s::Number, lon_s::Number, h_s::Number, θ::Number)

Check if the satellite with position vector `r_e` (ECEF) is inside the
visibility circle of a ground station with latitude `lat_s` [rad], longitude
`lon_s` [rad], altitude `h_s` (WGS-84), and a minimum elevation angle of
`θ` [rad].

Notice that the units of `r_e` and `h_s` must be the same.

Returns `true` if the satellite is inside the visibility circle, or `false`
otherwise.

"""
function ground_station_visible(r_e::AbstractVector, lat_s::Number,
                                lon_s::Number, h_s::Number, θ::Number)
    # Convert the ground station LLA to the ECEF frame.
    rs_e = GeodetictoECEF(lat_s, lon_s, h_s)

    return ground_station_visible(r_e, rs_e, θ)
end
