#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    SatelliteToolbox orbit propagator API: step
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export step!

"""
    function step!(orbp, Δt::Number)

Propagate the orbit in `orbp` by `Δt` s using the algorithm of `orbp`. The new
parameters will be written in `orbp`.

# Args

* `orbp`: Propagator structure.
* `Δt`: Step time [s].

# Returns

* The Keplerian elements represented in the inertial frame after the step (see
  `Orbit`) [SI units].
* The position vector represented in the inertial frame after the step [m].
* The velocity vector represented in the inertial frame after the step [m].

# Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. If the orbit parameters are obtained
from a TLE, then the inertial frame will be TEME. Notice, however, that the
perturbation theory requires an inertial frame with true equator.

"""
function step!(orbp::OrbitPropagatorJ2, Δt::Number)
    # Auxiliary variables.
    orb = orbp.orb
    j2d = orbp.j2d

    # Propagate the orbit.
    (r_i, v_i) = j2!(j2d, j2d.Δt + Δt)

    # Update the elements in the `orb` structure.
    orb.t += Δt/86400
    orb.a  = j2d.a_k*j2d.j2_gc.R0
    orb.e  = j2d.e_k
    orb.i  = j2d.i_k
    orb.Ω  = j2d.Ω_k
    orb.ω  = j2d.ω_k
    orb.f  = j2d.f_k

    # Return the information about the step.
    (copy(orbp.orb), r_i, v_i)
end

function step!(orbp::OrbitPropagatorSGP4{T}, Δt::Number) where T
    # Auxiliary variables.
    orb     = orbp.orb
    sgp4d   = orbp.sgp4d
    sgp4_gc = orbp.sgp4_gc

    # Propagate the orbit.
    (r_teme, v_teme) = sgp4!(sgp4d, sgp4d.Δt + Δt/60)

    # Convert km to m.
    r_teme *= 1000
    v_teme *= 1000

    # Update the elements in the `orb` structure.
    orb.t = sgp4d.epoch + sgp4d.Δt*60/86400
    orb.a = sgp4d.a_k*sgp4_gc.R0*1000
    orb.e = sgp4d.e_k
    orb.i = sgp4d.i_k
    orb.Ω = sgp4d.Ω_k
    orb.ω = sgp4d.ω_k
    orb.f = M_to_f(sgp4d.e_k, sgp4d.M_k)

    # Return the information about the step.
    (copy(orbp.orb), r_teme, v_teme)
end

function step!(orbp::OrbitPropagatorTwoBody, Δt::Number)
    # Auxiliary variables.
    orb = orbp.orb
    tbd = orbp.tbd

    # Propagate the orbit.
    (r_i, v_i) = twobody!(tbd, tbd.Δt + Δt)

    # Update the elements in the `orb` structure.
    orb.t += Δt/86400
    orb.f  = tbd.f_k

    # Return the information about the step.
    (copy(orbp.orb), r_i, v_i)
end
