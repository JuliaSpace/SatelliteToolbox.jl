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
#   Auxiliary functions related to coordinate transformations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-05-25: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Add support for TEME.
#
# 2018-05-13: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export FK5, ITRF, PEF, TOD, MOD, GCRF, J2000, TEME

# Auxiliary functions to define the conversion models.
@inline FK5()   = Val{:FK5}

# Auxiliary functions to define the reference frames.
@inline ITRF()  = Val{:ITRF}
@inline PEF()   = Val{:PEF}
@inline TOD()   = Val{:TOD}
@inline MOD()   = Val{:MOD}
@inline GCRF()  = Val{:GCRF}
@inline J2000() = Val{:J2000}
@inline TEME()  = Val{:TEME}
