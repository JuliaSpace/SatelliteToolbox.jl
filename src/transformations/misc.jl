#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Auxiliary functions related to coordinate transformations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export FK5, ITRF, PEF, TOD, MOD, GCRF, J2000, TEME, CIRS, TIRS, ERS

# Auxiliary functions to define the conversion models.
@inline FK5()   = Val(:FK5)

# Auxiliary functions to define the reference frames.
@inline ITRF()  = Val(:ITRF)
@inline PEF()   = Val(:PEF)
@inline TOD()   = Val(:TOD)
@inline MOD()   = Val(:MOD)
@inline GCRF()  = Val(:GCRF)
@inline J2000() = Val(:J2000)
@inline TEME()  = Val(:TEME)
@inline CIRS()  = Val(:CIRS)
@inline TIRS()  = Val(:TIRS)
@inline ERS()   = Val(:ERS)
