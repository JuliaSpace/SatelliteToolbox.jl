#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    Coordinate transformations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References:
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] ESA Navipedia: http://www.navipedia.net/
#
#   [3] mu-blox ag (1999). Datum Transformations of GPS Positions. Application
#       Note.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2016-03-10: Ronan Arraes Jardim Chagsa <ronan.arraes@inpe.br>
#   Add code to convert from LLA (WGS-84) to ECEF reference frame.
x
# 2015-11-23: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export ECEFtoLLA, LLAtoECEF, J2000toGMST, JDtoGMST

"""
### function ECEFtoLLA(r_e::Vector{Float64})

Convert the ECEF reference frame into LLA (WGS-84).

##### Inputs

* r_e: Position vector represented in the ECEF reference frame.

##### Outputs

* The position vector represented in the ECEF reference frame.

##### Remarks

Based on algorithm in [3].

"""

function ECEFtoLLA(r_e::Vector{Float64})
    # Auxiliary variables.
    X = r_e[1]
    Y = r_e[2]
    Z = r_e[3]

    p = sqrt(X^2 + Y^2)
    θ = atan2(Z*a_wgs84, p*b_wgs84)


    # Compute LLA.
    lon = atan2(Y, X)
    lat = atan2(Z + el_wgs84^2*b_wgs84*sin(θ)^3,
                p -  e_wgs84^2*a_wgs84*cos(θ)^3)
    N   = a_wgs84/sqrt(1 - e_wgs84^2*sin(lat)^2 )
    h   = p/cos(lat) - N

    (lat, lon, h)
end

"""
### function LLAtoECEF(lat::Real, lon::Real, h::Real)

Convert LLA (WGS-84) into the ECEF reference frame.

##### Inputs

* lat: Latitude in WGS-84 [rad].
* lon: Longitude in WGS-84 [rad].
* h: Altitude in WGS-84 [rad].

##### Outputs

* The position vector represented in the ECEF reference frame.

##### Remarks

Based on algorithm in [3].

"""

function LLAtoECEF(lat::Real, lon::Real, h::Real)
    # Radius of curvature [m].
    N = a_wgs85/sqrt(1 - e_wgs84^2*sin(lat)^2 )

    # Compute the position in ECEF frame.
    [ (         N + h)*cos(lat)*cos(lon);
      (         N + h)*cos(lat)*sin(lon);
      ( (b/a)^2*N + h)*sin(lon);]
end

"""
### function J2000toGMST(J2000::Real)

Compute the Greenwich Mean Sideral Time (GMST).

##### Args

* J2000: Date in J2000.0 reference (UT1).

##### Returns

* Greenwich mean sideral time [rad].

##### Remarks

Based on algorithm in [2] (http://www.navipedia.net/index.php/CEP_to_ITRF),
accessed at 2015-12-01.

"""

function J2000toGMST(J2000::Real)
    # Julian UT1 Date at 0h.
    JD_aux    = floor(J2000-0.5)
    JD_UT1_0h = (JD_aux <= 0) ? JD_aux + 0.5 : JD_aux - 0.5

    # UT1 in the Julian day [s].
    UT1 = (J2000 - JD_UT1_0h)*86400

	# Julian centuries elapsed from the epoch J2000.0.
	T_UT1 = JD_UT1_0h/36525;

	# Greenwich Mean Sideral Time at 0h [s].
    GMST_0h = 24110.54841 + 8640184.812866*T_UT1
		  	              + 0.093104*T_UT1^2
						  - 6.2e-6*T_UT1^3;

    # Greenwich Mean Sideral Time at T_UT1 [s].
    GMST = GMST_0h + 1.002737909350795*UT1

	# Greenwich Mean Sideral Time at provided date [rad].
    mod(GMST*pi/43200, 2*pi);
end

"""
### function JDtoGMST(JD::Real)

Compute the Greenwich Mean Sideral Time (GMST).

##### Inputs

* JD: Julian day.

##### Outputs

* Greenwich mean sideral time [rad].

##### Remarks

Based on algorithm in [1, pp. 188].

"""
function JDtoGMST(JD::Real)
	J2000toGMST(JD - JD_J2000);
end


