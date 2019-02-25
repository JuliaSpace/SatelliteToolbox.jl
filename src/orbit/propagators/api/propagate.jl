#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    SatelliteToolbox orbit propagator API: propagate!
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export propagate!

"""
    function propagate!(orbp, t::Number) where T
    function propagate!(orbp, t::AbstractVector) where T

If `t` is a number, then propagate `orbp` by `t` [s] from the orbit epoch.
Otherwise, if `t` is an array, then propagate the orbit in `orbp` using the time
instants defined in the vector `t` [s].

In both cases, the orbit propagator algorithm is the one related to the
structure `orbp`.

The structure `orbp` will contain the elements at the last propagation instant.

# Returns

* The mean Keplerian elements represented in inertial frame in each time instant
  (see `Orbit`) [SI units].
* The position vector represented in inertial frame in each time instant [m].
* The velocity vector represented in inertial frame in each time instant [m].

If `t` is an array, then those values will be an array containing the
information related to each epoch in `t`.

# Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. If the orbit parameters are obtained
from a TLE, then the inertial frame will be TEME. Notice, however, that the
perturbation theory requires an inertial frame with true equator.

"""
function propagate!(orbp::OrbitPropagatorJ2{T}, t::Number) where T
    # Auxiliary variables.
    orb = orbp.orb
    j2d = orbp.j2d

    # Propagate the orbit.
    (r_i, v_i) = j2!(j2d, t)

    # Update the elements in the `orb` structure.
    orb.t = j2d.epoch + t/86400
    orb.a = j2d.al_k*j2d.j2_gc.R0
    orb.e = j2d.e_k
    orb.i = j2d.i_k
    orb.Ω = j2d.Ω_k
    orb.ω = j2d.ω_k
    orb.f = j2d.f_k

    # Return.
    (copy(orb), r_i, v_i)
end

function propagate!(orbp::OrbitPropagatorJ4{T}, t::Number) where T
    # Auxiliary variables.
    orb = orbp.orb
    j4d = orbp.j4d

    # Propagate the orbit.
    (r_i, v_i) = j4!(j4d, t)

    # Update the elements in the `orb` structure.
    orb.t = j4d.epoch + t/86400
    orb.a = j4d.al_k*j4d.j4_gc.R0
    orb.e = j4d.e_k
    orb.i = j4d.i_k
    orb.Ω = j4d.Ω_k
    orb.ω = j4d.ω_k
    orb.f = j4d.f_k

    # Return.
    (copy(orb), r_i, v_i)
end

function propagate!(orbp::OrbitPropagatorSGP4{T}, t::Number) where T
    # Auxiliary variables.
    orb     = orbp.orb
    sgp4d   = orbp.sgp4d
    sgp4_gc = orbp.sgp4_gc

    # Propagate the orbit.
    (r_teme, v_teme) = sgp4!(sgp4d, t/60)

    # Convert km to m.
    r_teme *= 1000
    v_teme *= 1000

    # Update the elements in the `orb` structure.
    orb.t = sgp4d.epoch + t/86400
    orb.a = sgp4d.a_k*sgp4_gc.R0*1000
    orb.e = sgp4d.e_k
    orb.i = sgp4d.i_k
    orb.Ω = sgp4d.Ω_k
    orb.ω = sgp4d.ω_k
    orb.f = M_to_f(sgp4d.e_k, sgp4d.M_k)

    # Return.
    (copy(orb), r_teme, v_teme)
end

function propagate!(orbp::OrbitPropagatorTwoBody{T}, t::Number) where T
    # Auxiliary variables.
    orb = orbp.orb
    tbd = orbp.tbd

    # Propagate the orbit.
    (r_i, v_i) = twobody!(tbd, t)

    # Update the elements in the `orb` structure.
    orb.t = tbd.epoch + t/86400
    orb.f = tbd.f_k

    (copy(orb), r_i, v_i)
end

function propagate!(orbp::Union{OrbitPropagatorJ2{T},
                                OrbitPropagatorJ4{T},
                                OrbitPropagatorSGP4{T},
                                OrbitPropagatorTwoBody{T}},
                    t::AbstractVector) where T
    # Output.
    num_elems  = length(t)
    result_orb = Array{Orbit{T,T,T,T,T,T,T}}(undef, num_elems)
    result_r   = Array{SVector{3,T}}(undef, num_elems)
    result_v   = Array{SVector{3,T}}(undef, num_elems)

    for k = 1:num_elems
        # Propagate the orbit.
        (orb, r_i_k, v_i_k) = propagate!(orbp, t[k])

        result_orb[k] = copy(orb)
        result_r[k]   = r_i_k
        result_v[k]   = v_i_k
    end

    (result_orb, result_r, result_v)
end
