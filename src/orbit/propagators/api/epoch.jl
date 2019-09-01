#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    SatelliteToolbox orbit propagator API: epoch
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export epoch

"""
    function epoch(orbp)

Return the epoch of the propagator `orbp` [JD].

"""
epoch(orbp::OrbitPropagatorJ2)      = orbp.j2d.epoch
epoch(orbp::OrbitPropagatorJ4)      = orbp.j4d.epoch
epoch(orbp::OrbitPropagatorSGP4)    = orbp.sgp4d.epoch
epoch(orbp::OrbitPropagatorTwoBody) = orbp.tbd.epoch
