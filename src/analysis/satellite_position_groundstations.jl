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
#   Functions to verify if the satellite is within the station visibility
#   circle.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2016-03-10: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export satellite_check_station

"""
### function satellite_check_station(r_e::Vector{Float64}, rs_e::Vector{Float64}, minElev::Float64)

Check if the satellite is inside the visibility circle of a ground station.

##### Args

* r_e: Satellite position represented in the ECEF reference frame.
* rs_e: Ground station position represented in the ECEF reference frame.
* minElev: Minimum elevation angle in which the station can see the satellite
[rad].

##### Returns

* **TRUE**: The satellite is inside the visibility circle of the ground station.
* **FALSE**: The satellite is not inside the visibility circle of the ground
station.

"""

function satellite_check_station(r_e::Vector{Float64},
                                 rs_e::Vector{Float64},
                                 minElev::Float64)
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
### function satellite_check_station(r_e::Vector{Float64}, lat_s::Float64, lon_s::Float64, h_s::Float64, minElev::Float64)

Check if the satellite is inside the visibility circle of a ground station.

##### Args

* r_e: Satellite position represented in the ECEF reference frame.
* lat_s: Latitude of the ground station [rad].
* lon_s: Longitude of the ground station [rad].
* h_s: Altitude of the ground station (WGS-84) [m].
* minElev: Minimum elevation angle in which the station can see the satellite
[rad].

##### Returns

* **TRUE**: The satellite is inside the visibility circle of the ground station.
* **FALSE**: The satellite is not inside the visibility circle of the ground
station.

"""

function satellite_check_station(r_e::Vector{Float64},
                                 lat_s::Float64,
                                 lon_s::Float64,
                                 h_s::Float64,
                                 minElev::Float64)
    # Convert the ground station LLA to the ECEF frame.
    rs_e = LLAtoECEF(lat_s, lon_s, h_s)

    satellite_check_station(r_e, rs_e, minElev)
end
