# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#    API for the orbit propagators in SatelliteToolbox.jl.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export get_epoch, init_orbit_propagator, propagate!, propagate_to_epoch!, step!

"""
    init_orbit_propagator(T, args...; kwargs...)

Initialize the orbit propagator of type `T`. The arguments `args` and keywords
`kwargs` depends of the propagator type.

"""
init_orbit_propagator

"""
    get_epoch(orbp)

Return the epoch of the propagator `orbp` [JD].

"""
get_epoch

"""
    propagate!(orbp::OrbitPropagator{T}, t::Number) where T
    propagate!(orbp::OrbitPropagator{T}, t::AbstractVector) where T

If `t` is a number, then propagate `orbp` by `t` [s] from the orbit epoch.
Otherwise, if `t` is an array, then propagate the orbit in `orbp` using the time
instants defined in the vector `t` [s].

In both cases, the orbit propagator algorithm is the one related to the
structure `orbp`.

# Returns

* The position vector represented in inertial frame in each time instant [m].
* The velocity vector represented in inertial frame in each time instant [m].

If `t` is an array, then those values will be an array containing the
information related to each epoch in `t`.

"""
propagate!

function propagate!(orbp::OrbitPropagator{T}, t::AbstractVector) where T
    # Output.
    num_elems  = length(t)
    result_orb = Vector{KeplerianElements{T}}(undef, num_elems)
    result_r   = Vector{SVector{3,T}}(undef, num_elems)
    result_v   = Vector{SVector{3,T}}(undef, num_elems)

    @inbounds for k = 1:num_elems
        # Propagate the orbit.
        r_i_k, v_i_k = propagate!(orbp, t[k])
        result_r[k]  = r_i_k
        result_v[k]  = v_i_k
    end

    return result_r, result_v
end

"""
    propagate_to_epoch!(orbp::OrbitPropagator{T}, JD::Number) where T
    propagate_to_epoch!(orbp::OrbitPropagator{T}, JD::AbstractVector) where T

If `t` is a number, then propagate `orbp` until the epoch `JD` [Julian Day].
Otherwise, if `JD` is an array, then propagate the orbit in `orbp` using the
epochs defined in the vector `t` [Julian Day].

In both cases, the orbit propagator algorithm is the one related to the
structure `orbp`.

# Returns

* The position vector represented in inertial frame in each time instant [m].
* The velocity vector represented in inertial frame in each time instant [m].

If `JD` is an array, then those values will be an array containing the
information related to each epoch in `JD`.

"""
propagate_to_epoch!(orbp::OrbitPropagator,
                    JD::Union{Number, AbstractVector}) where T =
    propagate!(orbp, (JD .- get_epoch(orbp))*86400)

"""
    step!(orbp::OrbitPropagator{T}, Δt::Number)

Propagate the orbit in `orbp` by `Δt` [s] using the algorithm of `orbp`. The new
parameters will be written in `orbp`.

# Returns

* The position vector represented in the inertial frame after the step [m].
* The velocity vector represented in the inertial frame after the step [m].

"""
step!

################################################################################
#                              Iterator interface
################################################################################

# There functions allow broadcast when using the orbit propagators.
Base.iterate(orbp::OrbitPropagator) = (orbp, nothing)
Base.iterate(orbp::OrbitPropagator, ::Nothing) = nothing
Base.length(orbp::OrbitPropagator) = 1
Base.eltype(orbp::T) where T<:OrbitPropagator = T
