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

export satellite_eclipse_time

"""
    function satellite_eclipse_time(JD0::Number, a::Number, e::Number, i::Number, w::Number, RAAN::Number, numDays::Integer, relative::Bool = false)

Compute the eclipse time of a satellite.

# Args

* `JD0`: Initial instant for the analysis [Julian day].
* `a`: Semi-major axis of the orbit [m].
* `e`: Orbit eccentricity.
* `i`: Orbit inclination [rad].
* `w`: Argument of perigee [rad].
* `RAAN`: Right ascension of the ascending node at `JD0` [rad].
* `numDays`: Number of days of the analysis.
* `relative`: Compute the eclipse time relative to the nodal period.

# Returns

The following table:

        day | Sunlight Time | Penumbra Time | Umbra Time
       -----+---------------+---------------+------------

"""
function satellite_eclipse_time(JD0::Number,
                                a::Number,
                                e::Number,
                                i::Number,
                                w::Number,
                                RAAN::Number,
                                numDays::Integer,
                                relative::Bool = false)
    # Constants
    deg2rad = pi/180.0
    rad2deg = 180.0/pi
    day2sec = 24.0*60.0*60.0

    # Step of the orbit propagation (mean anomaly) [rad].
    step = 0.1*deg2rad

    # Initialization of variables.
    theta = 0.0                   # Sun angle relative to the inertial
                                  # coordinate frame.

    days = collect(0:1:numDays-1) # Vector of the days in which the beta angle
                                  # will be computed.

    # Mean anomaly.
    M = collect(0:step:2*pi)

    # Penumbra time.
    p_time = zeros(numDays)

    # Umbra time.
    u_time = zeros(numDays)

    # Sunlight time.
    s_time = zeros(numDays)

    # Semi-lactum rectum.
    p = a*(1.0-e^2)

    # Angular velocity of the orbit [rad/s].
    n = angvel(a, e, i, :J2)

    # Step in time
    tstep = step/n

    # Perturbations.
    #
    # RAAN rotation rate [rad/s].
    dOmega = dRAAN(a, e, i, :J2)

    # Perturbation of the argument of perigee [rad/s].
    dw = dArgPer(a, e, i, :J2)

    # Loop.
    for d in days
        # Get the sun position represented in the Inertial coordinate frame.
        s_i = sun_position_i(JD0+d)

        # Compute the new orbit parameters due to perturbations.
        w_d = w + dw*(d*day2sec)
        RAAN_d = RAAN + dOmega*(d*day2sec)

        for m_k in M
            # Get the satellite position vector represented in the Inertial coordinate
            # frame.
            f = M_to_f(e, m_k)

            (r_i, rt_i) = satellite_position_i(a, e, i, RAAN_d, w_d, f)

            # Check the lighting conditions.
            lighting = satellite_lighting_condition(r_i, s_i)

            if (lighting == SAT_LIGHTING_SUNLIGHT)
                s_time[d+1] += tstep
            elseif (lighting == SAT_LIGHTING_UMBRA)
                u_time[d+1] += tstep
            else
                p_time[d+1] += tstep
            end
        end
    end

    if (!relative)
        [days s_time p_time u_time]
    else
        total_time = s_time + p_time + u_time
        [days s_time./total_time p_time./total_time u_time./total_time]
    end
end
