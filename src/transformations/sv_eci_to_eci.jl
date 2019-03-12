#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Convert a satellite state vector from an Earth-Centered Inertial (ECI)
#   reference frame to another ECI reference frame.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export svECItoECI

"""
    function svECItoECI(sv::SatelliteStateVector, args...)

Convert the satellite state vector `sv` from an ECI frame to another ECI frame.
The arguments `args...` must match those of the function `rECItoECI` **wihtout**
the rotation representation.

"""
function svECItoECI(sv::SatelliteStateVector, args...)
    D = rECItoECI(DCM, args...)

    # Since both frames does not have a significant angular velocity between
    # them, then we just need to convert the representations.
    return SatelliteStateVector(t = sv.t, r = D*sv.r, v = D*sv.v, a = D*sv.a)
end
