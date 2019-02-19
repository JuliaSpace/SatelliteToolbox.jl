#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions to compute the satellite ground trace.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export ground_trace

"""
    function ground_trace(orbp::OrbitPropagator{N}, eop_data::Union{Nothing, EOPData_IAU1980, EOPData_IAU2000A} = nothing; ECI = TEME(), ECEF = PEF(), span = 1.0) where N

Compute the ground trace of the object with orbit defined by `orbp`.

By default, it considers that the orbit elements on the propagator are
represented in the True Equator, Mean Equinox (TEME) reference frame and the
ground trace will be computed in the Pseudo-Earth Fixed (PEF) reference frame.
Hence, no EOP data is needed. However, this can be changed by the keywords
presented as follows.

# Keywords

* `eop_data`: EOP data that will be used to convert the ECI reference frame to
              the ECEF reference frame. If `nothing`, then it will not be used
              (see `rECItoECEF`). (**Default** = `nothing`)
* `ECI`: ECI frame in which the orbit elements in `orbp` are represented.
         (**Default** = `TEME()`)
* `ECEF`: ECEF frame that will be used to compute the ground trace.
          (**Default** = `PEF()`)
* `span`: Defines for how much time the ground trace will be computed. The unit
          is the orbit period. (**Default** = 1.0)
* `dt`: Time interval between two samples [s]. (**Default** = 10.0)

# Returns

A vector of tuples with the pairs `(latitude,longitude)` of the ground trace.

"""
function ground_trace(orbp::OrbitPropagator{N};
                      eop_data::Union{Nothing, EOPData_IAU1980, EOPData_IAU2000A} = nothing,
                      ECI = TEME(), ECEF = PEF(), span = 1.0, dt = 10.0) where N

    # Get the orbit period.
    T = period(orbp.orb, :J2)

    # Copy orbit structure so that it is not modified by `propagate`.
    orbp_c = deepcopy(orbp)

    # Compute the points represented in the inertial reference frame.
    (o, r_i, ~) = propagate!(orbp_c, 0.0:dt:T*span)

    # Get the epochs in Julian Day of each instant.
    JD = map(x->x.t, o)

    # Convert from the ECI to the ECEF frame.
    if eop_data == nothing
        r_e = map( (t,v_i)->rECItoECEF(ECI, ECEF, t)*v_i, JD, r_i )
    else
        r_e = map( (t,v_i)->rECItoECEF(ECI, ECEF, t, eop_data)*v_i, JD, r_i )
    end

    # Convert to Geodetic coordinates.
    geod = ECEFtoGeodetic.(r_e)
    return map(x->(x[1],x[2]), geod)
end

