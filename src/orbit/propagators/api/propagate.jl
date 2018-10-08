#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    SatelliteToolbox orbit propagator API: propagate!
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export propagate!

"""
    function propagate!(orbp, t::AbstractVector) where T

Propagate the orbit in `orbp` using the time instants defined in the vector `t`
using the algorithm related to the structure `orbp`. The structure `orbp` will
contain the elements at the last propagation instant.

# Args

* `orbp`: Propagator structure.
* `t`: Time instants from orbit epoch in which the orbit will be propagated
       [s].

# Returns

* An array with the mean Keplerian elements represented in inertial frame in
  each time instant (see `Orbit`) [SI units].
* An array with the position vector represented in inertial frame in each time
  instant [m].
* An array with the velocity vector represented in inertial frame in each time
  instant [m].

# Remarks

The inertial frame in which the output is represented depends on which frame it
was used to generate the orbit parameters. If the orbit parameters are obtained
from a TLE, then the inertial frame will be TEME. Notice, however, that the
perturbation theory requires an inertial frame with true equator.

"""
function propagate!(orbp::OrbitPropagatorJ2{T}, t::AbstractVector) where T
    # Auxiliary variables.
    orb = orbp.orb
    j2d = orbp.j2d

    # Output.
    result_orb = Array{Orbit{T,T,T,T,T,T,T}}(undef, 0)
    result_r   = Array{Vector{T}}(undef, 0)
    result_v   = Array{Vector{T}}(undef, 0)

    for k in t
        # Propagate the orbit.
        (r_i_k, v_i_k) = j2!(j2d, k)

        # Update the elements in the `orb` structure.
        orb.t = j2d.epoch + k/86400
        orb.a = j2d.a_k*j2d.j2_gc.R0
        orb.e = j2d.e_k
        orb.i = j2d.i_k
        orb.Ω = j2d.Ω_k
        orb.ω = j2d.ω_k
        orb.f = j2d.f_k

        push!(result_orb, copy(orb))
        push!(result_r,   r_i_k)
        push!(result_v,   v_i_k)
    end

    (result_orb, result_r, result_v)
end

function propagate!(orbp::OrbitPropagatorSGP4{T}, t::AbstractVector) where T
    # Auxiliary variables.
    orb     = orbp.orb
    sgp4d   = orbp.sgp4d
    sgp4_gc = orbp.sgp4_gc

    # Output.
    result_orb = Array{Orbit{T,T,T,T,T,T,T}}(undef, 0)
    result_r   = Array{Vector{T}}(undef, 0)
    result_v   = Array{Vector{T}}(undef, 0)

    for k in t
        # Propagate the orbit.
        (r_teme_k, v_teme_k) = sgp4!(sgp4d, k/60)

        # Convert km to m.
        r_teme_k *= 1000
        v_teme_k *= 1000

        # Update the elements in the `orb` structure.
        orb.t = sgp4d.epoch + k/86400
        orb.a = sgp4d.a_k*sgp4_gc.R0*1000
        orb.e = sgp4d.e_k
        orb.i = sgp4d.i_k
        orb.Ω = sgp4d.Ω_k
        orb.ω = sgp4d.ω_k
        orb.f = M_to_f(sgp4d.e_k, sgp4d.M_k)

        push!(result_orb, copy(orb))
        push!(result_r,   r_teme_k)
        push!(result_v,   v_teme_k)
    end

    (result_orb, result_r, result_v)
end

function propagate!(orbp::OrbitPropagatorTwoBody{T}, t::AbstractVector) where T
    # Auxiliary variables.
    orb = orbp.orb
    tbd = orbp.tbd

    # Output.
    result_orb = Array{Orbit{T,T,T,T,T,T,T}}(undef, 0)
    result_r   = Array{Vector{T}}(undef, 0)
    result_v   = Array{Vector{T}}(undef, 0)

    for k in t
        # Propagate the orbit.
        (r_i_k, v_i_k) = twobody!(tbd, k)

        # Update the elements in the `orb` structure.
        orb.t = tbd.epoch + k/86400
        orb.f = tbd.f_k

        push!(result_orb, copy(orb))
        push!(result_r,   r_i_k)
        push!(result_v,   v_i_k)
    end

    (result_orb, result_r, result_v)
end
