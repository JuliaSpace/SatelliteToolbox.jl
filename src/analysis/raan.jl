#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions related to the analysis of the Right Ascension of the Ascending
#   Node (RAAN).
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export compute_RAAN_lt, sim_RAAN_J2

"""
    function compute_RAAN_lt(JD::Number, asc_node_lt::Number)

Compute the RAAN (0,2π) \\[rad] so that the orbit plane local time is
`asc_node_lt` [hour] at the Julian day `JD`.

"""
function compute_RAAN_lt(JD::Number, asc_node_lt::Number)
    # Get the sun position at noon (UT) represented in the Inertial ref. frame.
    s_i = sun_position_i(JD)

    # Get the desired angle between the Sun and the ascending node [deg].
    alpha = (asc_node_lt-12)*π/12

    # Get the right ascension of the Sun in the Inertial ref. frame. This is the
    # Sun apparent local time.
    SALT = atan(s_i[2],s_i[1])

    # Get the equation of time to compute the Sun mean local time [rad].
    eot = equation_of_time(JD)

    # Compute the Sun mean local time.
    SMLT = SALT + eot

    # Compute the desired RAAN in the interval 0, 2*pi.
    RAAN = mod(SMLT+alpha, 2π)
end

"""
    function sim_RAAN_J2(a::Number, e::Number, i::Number, RAAN_0::Number, numDays::Integer)

Simulate the RAAN of an orbit with semi-major axis `a` [m], eccentricity `e`,
inclination `i` [rad], and initial RAAN `RAAN_0` [rad] considering J2
perturbations. The analysis is performed for `numDays` days.

# Returns

A `numDays` × 2 matrix in which the i-th line is:

    | day | RAAN (0,2π) [rad] |

"""
function sim_RAAN_J2(a::Number,
                     e::Number,
                     i::Number,
                     RAAN_0::Number,
                     numDays::Integer)

    # Initialization of variables.
    days = collect(0:1:numDays-1) # Vector of the days in which the RAAN will be
                                  # simulated.

    # RAAN rotation rate [rad/day].
    dOmega = dRAAN(a, e, i, :J2)*24*3600

    # Simulate the RAAN for each day considering just the J2 perturbations.
    RAAN = mod.(RAAN_0 .+ dOmega.*days, 2π)

    [days RAAN]
end
