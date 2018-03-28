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
#   Definition of constants.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorne, CA.
#
#   [2] http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html
#       Accessed on 2017-08-07.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2017-08-07: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version. Constants were moved from `SatToolbox.jl` file.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# Julian Day of J2000.0 epoch.
const JD_J2000 = 2451545.0

# Earth Equatorial radius [m].
const R0 = 6378137.0

# Earth mean radius [m].
const Rm = 6371009.0

# Standard gravitational parameter for Earth [m^3/s^2].
const m0 = 3.986004418e14

# Perturbation terms based on EGM-08 standard gravitational model [1, pp. 1039].
const J2 = 1.08262617385222e-3
const J3 = 2.53241051856772e-6
const J4 = 1.61989759991697e-6

# Sun radius [m].
const Rs = 6.963e8

# Earth's orbit mean motion [rad/s]
const ne = (360.0/365.2421897)*pi/180/86400

# Conversion factor from AU to m.
const au2m = 149597870700.0

# Sun radiation emitted [J/sec].
const sunRad = 3.826e26

# WGS-84 Data [2].
const a_wgs84  = 6378137.0
const f_wgs84  = 1/298.257223563
const b_wgs84  = a_wgs84*(1-f_wgs84)
const e_wgs84  = sqrt( (a_wgs84^2-b_wgs84^2)/a_wgs84^2 )
const el_wgs84 = sqrt( (a_wgs84^2-b_wgs84^2)/b_wgs84^2 )
