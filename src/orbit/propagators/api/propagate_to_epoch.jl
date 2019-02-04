#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    SatelliteToolbox orbit propagator API: propagate_to_epoch!
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export propagate_to_epoch!

"""
    function propagate_to_epoch!(orbp, JD::Number) where T
    function propagate_to_epoch!(orbp, JD::AbstractVector) where T

If `t` is a number, then propagate `orbp` until the epoch `JD` [Julian Day].
Otherwise, if `JD` is an array, then propagate the orbit in `orbp` using the
epochs defined in the vector `t` [Julian Day].

In both cases, the orbit propagator algorithm is the one related to the
structure `orbp`.

The structure `orbp` will contain the elements at the last propagation instant.

# Returns

* The mean Keplerian elements represented in inertial frame in each time instant
  (see `Orbit`) [SI units].
* The position vector represented in inertial frame in each time instant [m].
* The velocity vector represented in inertial frame in each time instant [m].

If `JD` is an array, then those values will be an array containing the
information related to each epoch in `JD`.

# Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. If the orbit parameters are obtained
from a TLE, then the inertial frame will be TEME. Notice, however, that the
perturbation theory requires an inertial frame with true equator.

"""
function propagate_to_epoch!(orbp::OrbitPropagatorJ2{T},
                             JD::Union{Number,AbstractVector}) where T
    propagate!(orbp, (JD .- orbp.j2d.epoch)*86400)
end

function propagate_to_epoch!(orbp::OrbitPropagatorSGP4{T},
                             JD::Union{Number,AbstractVector}) where T
    propagate!(orbp, (JD .- orbp.sgp4d.epoch)*86400)
end

function propagate_to_epoch!(orbp::OrbitPropagatorTwoBody{T},
                             JD::Union{Number,AbstractVector}) where T
    propagate!(orbp, (JD .- orbp.tbd.epoch)*86400)
end
