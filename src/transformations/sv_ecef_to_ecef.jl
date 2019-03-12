#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Convert a satellite state vector from an Earth-Centered, Earth-Fixed (ECEF)
#   reference frame to another ECEF reference frame.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export svECEFtoECEF

"""
    function svECEFtoECEF(sv::SatelliteStateVector, args...)

Convert the satellite state vector `sv` from an ECEF frame to another ECEF
frame. The arguments `args...` must match those of the function `rECEFtoECEF`
**wihtout** the rotation representation.

"""
function svECEFtoECEF(sv::SatelliteStateVector, args...)
    D = rECEFtoECEF(DCM, args...)

    # Since both frames does not have a significant angular velocity between
    # them, then we just need to convert the representations.
    return SatelliteStateVector(t = sv.t, r = D*sv.r, v = D*sv.v, a = D*sv.a)
end
