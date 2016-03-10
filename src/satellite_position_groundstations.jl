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
### function satellite_check_station(r_e::Vector{Float64}, lat_s::Float64, lon_s::Float64, h_s::Float64, minElev::Float64)

Check if the satellite is inside the visibility circle of a ground station.

##### Args

* r_e: Satellite position represented in the ECEF reference frame.
* lat_s: Latitude of the ground station [rad].
* lon_s: Longitude of the ground station [rad].
* h_s: Altitude of the ground station (

##### Returns

* Greenwich mean sideral time [rad].

##### Remarks

Based on algorithm in [2] (http://www.navipedia.net/index.php/CEP_to_ITRF),
accessed at 2015-12-01.

"""

function satellite_check_station(r_e::Vector{Float64},
                                 lat_s::Float64,
                                 lon_s::Float64,
                                 h_s::Float64,
                                 minElev::Float64)
    # Convert the station to the ECEF frame.
    Des  = angle2dcm(lat_s, -lon_s, 0.0, "YZX")
    rs_e = Des*[R0+h_s; 0.0; 0.0]

    # Check if the satellite is within the visibility circle of the station.
    dr_e = r_e - rs_e
    cos_beta = ( ( dr_e/norm(dr_e) )'*( rs_e/norm(rs_e) ) )[1]

    display(cos_beta)

    if (cos_beta > cos(pi/2-minElev))
        return true
    else
        return false
    end
end
