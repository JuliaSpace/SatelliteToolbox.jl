#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    This function computes the sunlight, penumbra, and umbra times of the
#    satellite for one orbit every day in a year.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Longo, C. R. O., Rickman, S. L (1995). Method for the Calculation of
#       Spacecraft Umbra and Penumbra Shadow Terminator Points. NASA Technical
#       Paper 3547.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export eclipse_time_summary

"""
    eclipse_time_summary(JD₀::Number, a::Number, e::Number, i::Number, RAAN::Number, w::Number, Δd::Integer, relative::Bool = false, Δt₀::AbstractFloat = -1.0)

Compute the eclipse time of an orbit with semi-major axis `a` [m], eccentricity
`e`, inclination `i` [rad], initial right ascension of the ascending node `RAAN`
[rad], and initial argument of perigee `w` [rad]. The orbit epoch, which is also
the day in which the analysis will begin, is `JD₀` [Julian Day]. The analysis
will be performed for each day during `Δd` days.

This function will compute the eclipse time of one orbit per day.

If the argument `relative` is `true`, then the computed times will be relative
to the nodal period [%]. Otherwise, they will be computed in seconds. By
default, `relative = false`.

The argument `Δt₀` can be used to select the time step in which the orbit will
be propagated. Notice that this algorithm performs a numerical search to find
the beginning of each section (sunlight, penumbra, and umbra) with millisecond
precision. Thus, selecting a high number for `Δt₀` will make the analysis
faster, but the accuracy is lost if a region time span is smalled than `Δt₀`. If
this parameter is omitted or if it is negative, then the time step will be
selected automatically to match a mean anomaly step of 5°.

All the analysis is performed using a J2 orbit propagator.

# Returns

The following table:

        day | Sunlight Time | Penumbra Time | Umbra Time
       -----+---------------+---------------+------------

"""
function eclipse_time_summary(JD₀::Number, a::Number, e::Number, i::Number,
                              RAAN::Number, w::Number, Δd::Integer,
                              relative::Bool = false, Δt₀::AbstractFloat = -1.0)

    # Get the period of the orbit.
    T = period(a, e, i, :J2)

    # If the user did not specify the integration step, then compute based on
    # a mean anomaly step.
    if Δt₀ <= 0
        # Step of the orbit propagation (mean anomaly) [rad].
        Δm = 0.5*π/180

        # Time step related to the desired mean anomaly step [s].
        n   = angvel(a, e, i, :J2)
        Δt₀ = Δm/n
    end

    # Vector of the days in which the beta angle will be computed.
    days = collect(0:1:Δd-1)

    # Preallocate the output variables.
    p_time = zeros(Δd)  # Penumbra time.
    u_time = zeros(Δd)  # Umbra time.
    s_time = zeros(Δd)  # Sunlight time.

    # Configure the orbit propagator.
    orbp = init_orbit_propagator(Val(:J2), JD₀, a, e, i, RAAN, w, 0)

    # Lambda functions
    # ==========================================================================

    # Function to compute the lightning condition given an instant of the day
    # `t`, the day from orbit epoch `d`, and the Sun vector `s_i`.
    f(t,d,s_i)::Int = begin
        ~, r_i, ~ = propagate!(orbp, 86400d + t)
        return satellite_lighting_condition(r_i, s_i)
    end

    # Return `true` if the lightning condition given an instant of the day `t`,
    # the day from orbit epoch `d`, and the Sun vector `s_i` is equal `state`.
    fb(t,d,s_i,state)::Bool = f(t,d,s_i) == state

    # Accumulate the time step `Δt` according to the state `state`.
    accum(Δt, state, ind, s, p, u)::Nothing = begin
        @inbounds if state == SAT_LIGHTING_SUNLIGHT
            s[ind] += Δt
        elseif state == SAT_LIGHTING_PENUMBRA
            p[ind] += Δt
        elseif state == SAT_LIGHTING_UMBRA
            u[ind] += Δt
        end

        return nothing
    end

    # Loop
    # ==========================================================================

    @inbounds for d in days
        # Get the sun position represented in the Inertial coordinate frame.
        s_i = sun_position_i(JD₀+d)

        # Initial state.
        state = f(0,d,s_i)

        # Compute the eclipse time during one orbit.
        Δt = Δt₀
        k  = Δt

        while true
            new_state = f(k,d,s_i)

            # Check if the state has changed.
            if new_state != state
                # Refine to find the edge.
                k₀ = k - Δt
                k₁ = k
                kc = find_crossing(fb, k₀, k₁, true, false, d, s_i, state)

                # Times to be added in the previous and current states.
                Δ₀ = kc - k₀
                Δ₁ = k₁ - kc

                accum(Δ₀, state,     d+1, s_time, p_time, u_time)
                accum(Δ₁, new_state, d+1, s_time, p_time, u_time)

            # If not, just add the time step to the current state.
            else
                accum(Δt, state, d+1, s_time, p_time, u_time)
            end

            state = new_state

            abs(k - T) < 1e-3 && break

            # Make sure that the last interval will have the exact size so that
            # the end of analysis is the end of the orbit.
            (k + Δt > T) && (Δt = T - k)

            k += Δt
        end
    end

    if (!relative)
        [days s_time p_time u_time]
    else
        total_time = s_time + p_time + u_time
        [days s_time./total_time p_time./total_time u_time./total_time]
    end
end
