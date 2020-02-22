# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions to convert osculating elements to mean elements.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A; Crawford, P (2008). SGP4 orbit determination. AIAA/AAS
#       Astrodynamics Specialist Conference, Honoulu, HI.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export rv_to_mean_elements_sgp4, rv_to_tle

"""
    rv_to_mean_elements_sgp4(vJD::AbstractVector{T}, vr::AbstractVector{Tv}, vv::AbstractVector{Tv}, W = I; max_it = 10000, sgp4_gc = sgp4_gc_wgs84) where {T,Tv<:AbstractVector}

Compute the mean elements for SGP4 based on the position `vr` and velocity
vectors `vr` represented in TEME reference frame. The epoch of those
measurements [Julian Day] must be in `vJD`.

The matrix `W` defined the weights for the least-square algorithm.

The variable `max_it` defines the maximum allowed number of iterations.

The variable `sgp4_gc` defines which constants should be used when running SGP4.

# Returns

* The epoch of the elements [Julian Day].
* The mean elements for SGP4 algorithm:
    - Mean motion [rad/s];
    - Eccentricity [];
    - Inclination [rad];
    - Right ascension of the ascending node [rad];
    - Argument of perigee [rad];
    - Mean anomaly [rad];
    - BSTAR.
* The covariance matrix of the mean elements estimation.

"""
function rv_to_mean_elements_sgp4(vJD::AbstractVector{T}, vr::AbstractVector{Tv},
                                  vv::AbstractVector{Tv}, W = I; max_it = 10000,
                                 sgp4_gc = sgp4_gc_wgs84) where
    {T,Tv<:AbstractVector}

    # Number of measurements.
    num_meas   = length(vr)

    # Number of states that will be fitted. In this cases they are:
    #
    #   - Mean motion;
    #   - Eccentricity;
    #   - Inclination;
    #   - Right ascension of the ascending node;
    #   - Argument of perigee;
    #   - Mean anomaly;
    #   - BSTAR.
    num_states = 7

    # Initial guess of the mean elements.
    #
    # NOTE: x₁ is the previous estimate and x₂ is the current estimate.
    #
    # TODO: Should we compute the mean motion using a better algorithm?
    o = rv_to_kepler(vr[end], vv[end])
    x₁ = @SVector T[sqrt(m0/o.a^3), o.e, o.i, o.Ω, o.ω, f_to_M(o.e, o.f), 0]
    x₂ = x₁
    P  = SMatrix{num_states,num_states,T}(Matrix(I,num_states,num_states))

    # Loop until the maximum allowed iteration.
    @inbounds for it = 1:max_it
        x₁ = x₂

        # Variables to store the summations to compute the least square fitting
        # algorithm.
        ΣAᵀWA = @SMatrix zeros(num_states, num_states)
        ΣAᵀWb = @SVector zeros(num_states)

        # We will interpolate backwards in time, so that the mean elements
        # in x₂ are related to the newest measurement.
        for k = num_meas:-1:1
            # Obtain the measured ephemerides.
            y = vcat(vr[k],vv[k])

            # Obtain the computed ephemerides considering the current estimate
            # of the mean elements.
            Δt = (vJD[k] - vJD[end])*86400
            r̂, v̂, ~ = _sgp4_si(Δt, sgp4_gc, vJD[k], x₁...)
            ŷ = vcat(r̂,v̂)

            # Compute the error.
            b = y-ŷ

            # Compute the Jacobian.
            A = _sgp4_jacobian(Δt, vJD[end], x₁)

            # Accumulation.
            ΣAᵀWA += A'*W*A
            ΣAᵀWb += A'*W*b
        end

        # Update the estimate.
        P  = pinv(ΣAᵀWA)
        δx = P*ΣAᵀWb
        x₂ = x₁ + δx

        # Make sure that all parameters are inside the allowed intervals.
        if x₂[2] < 0
            x₂ = setindex(x₂, -x₂[2], 2)
        end

        # Obtain the maximum deviation.
        σ_max = maximum( abs.( (x₂-x₁) ) )

        println("PROGRESS: $it / $max_it, σ_max = $σ_max")

        # If the threshold is achieved, then just stop the loop.
        σ_max < 1e-6 && break
    end

    # Return the mean elements for SGP4 and the covariance matrix.
    return vJD[end], x₂, P
end

"""
    rv_to_tle(args...; name::String = "UNDEFINED", sat_num::Int = 9999, classification::Char = 'U', int_designator = "999999", elem_set_number::Int = 0, rev_num, kwargs...)

Convert a set of position and velocity vectors represented in TEME reference
frame to a TLE. The arguments `args` and keywords `kwargs` are the same as those
described in the function `rv_to_mean_elements_sgp4`.

Additionally, the user can specify some parameters of the generated TLE.

This function prints the TLE to `stdout` using the function `print_tle` and also
returns the TLE string.

"""
function rv_to_tle(args...;
                   name::String = "UNDEFINED",
                   sat_num::Int = 9999,
                   classification::Char = 'U',
                   int_designator = "999999",
                   elem_set_number::Int = 0,
                   rev_num::Int = 0,
                   kwargs...)

    # Convert the position and velocity vectors to mean elements.
    JD,x,P = rv_to_mean_elements_sgp4(args...; kwargs...)

    # Compute the data as required by the TLE format.
    dt  = JDtoDate(DateTime, JD)
    dt₀ = JDtoDate(DateTime, DatetoJD(year(dt), 1, 1, 0, 0, 0))

    dt_year    = year(dt)
    epoch_year = dt_year < 1980 ? dt_year - 1900 : dt_year - 2000
    epoch_day  = (dt-dt₀).value/1000/86400 + 1

    # Construct the TLE.
    tle = TLE(name,
              sat_num,
              classification,
              int_designator,
              epoch_year,
              epoch_day,
              JD,
              0.0,
              0.0,
              x[7],
              elem_set_number,
              0,
              x[3]*180/pi,
              x[4]*180/pi,
              x[2],
              x[5]*180/pi,
              x[6]*180/pi,
              x[1]*43200/pi,
              rev_num,
              0)

    tle_str = tle_to_str(tle)

    # Print the TLE to the `stdout`.
    print(tle_str)

    # Return the TLE string.
    return tle_str
end


################################################################################
#                              Private functions
################################################################################

# Compute the SGP4 algorithm considering all variables in SI.
function _sgp4_si(Δt, sgp4_gc, epoch, n_0, e_0, i_0, Ω_0, ω_0, M_0, bstar)
    r,v,sgp4d = sgp4(Δt/60, sgp4_gc, epoch, 60n_0, e_0, i_0, Ω_0, ω_0, M_0, bstar)
    return 1000r, 1000v, sgp4d
end

function _sgp4_jacobian(Δt, JD₀, x₀::AbstractVector{T}; ϵ = 1e-10, sgp4_gc = sgp4_gc_wgs84) where T
    num_states = 7
    dim_obs    = 6
    M          = Matrix{T}(undef, dim_obs, num_states)

    # Auxiliary variables.
    x₁ = copy(x₀)
    x₂ = copy(x₀)

    @inbounds for i = 1:dim_obs, j = 1:num_states
        x₁     = setindex(x₁, x₁[j] - ϵ, j)
        x₂     = setindex(x₂, x₂[j] + ϵ, j)

        r,v,~  = _sgp4_si(Δt, sgp4_gc, JD₀, x₁...)
        f_x₁   = vcat(r,v)
        r,v,~  = _sgp4_si(Δt, sgp4_gc, JD₀, x₂...)
        f_x₂   = vcat(r,v)

        M[i,j] = (f_x₂[i]-f_x₁[i])/(2ϵ)

        x₁     = setindex(x₁, x₀[j], j)
        x₂     = setindex(x₂, x₀[j], j)
    end

    return M
end
