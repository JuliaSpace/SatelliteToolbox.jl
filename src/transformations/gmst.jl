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
#    Compute the Greenwich Mean Sideral Time (GMST).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-03-27: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version, but the functions were copied from another file. The change
#   occurred to improve the source code organization.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export J2000toGMST, JDtoGMST

"""
### function J2000toGMST(J2000::Number)

Compute the Greenwich Mean Sideral Time (GMST) given the instant `J2000`
represented in J2000.0 reference frame (UT1).

##### Args

* J2000: Instant in J2000.0 reference (UT1).

##### Returns

* Greenwich mean sideral time [rad].

##### Remarks

Based on algorithm in [2] (http://www.navipedia.net/index.php/CEP_to_ITRF),
accessed at 2015-12-01.

"""

function J2000toGMST(J2000::Number)
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
### function JDtoGMST(JD::Number)

Compute the Greenwich Mean Sideral Time (GMST) for a Julian Day `JD`.

##### Args

* JD: Julian day.

##### Returns

* Greenwich mean sideral time [rad].

##### Remarks

Based on algorithm in [1, pp. 188].

"""
function JDtoGMST(JD::Number)
	J2000toGMST(JD - JD_J2000);
end

