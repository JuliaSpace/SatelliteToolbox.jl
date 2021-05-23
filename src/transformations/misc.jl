# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Auxiliary functions related to coordinate transformations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export ITRF, PEF, TOD, MOD, GCRF, J2000, TEME, CIRS, TIRS, ERS, MOD06, MJ2000

# Auxiliary functions to define the reference frames.
@inline ITRF()   = Val(:ITRF)
@inline PEF()    = Val(:PEF)
@inline TOD()    = Val(:TOD)
@inline MOD()    = Val(:MOD)
@inline GCRF()   = Val(:GCRF)
@inline J2000()  = Val(:J2000)
@inline TEME()   = Val(:TEME)
@inline CIRS()   = Val(:CIRS)
@inline TIRS()   = Val(:TIRS)
@inline ERS()    = Val(:ERS)
@inline MOD06()  = Val(:MOD06)
@inline MJ2000() = Val(:MJ2000)
