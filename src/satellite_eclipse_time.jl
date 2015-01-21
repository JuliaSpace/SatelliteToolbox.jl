# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#     This function computes the sunlight, penumbra, and umbra times of the
#     satellite for one orbit every day in a year.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#    [1] Longo, C. R. O., Rickman, S. L (1995). Method for the Calculation of
#        Spacecraft Umbra and Penumbra Shadow Terminator Points. NASA Technical
#        Paper 3547.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2014-07-28: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#==#
# 
# @brief Compute the eclipse time of a satellite.
#
# @param[in] t0 Launch date [number of days since 01/01/2000].
# @param[in] a Semi-major axis of the satellite orbit [m].
# @param[in] e Orbit eccentricity.
# @param[in] i Orbit inclination [rad].
# @param[in] w Argument of perigee [rad].
# @param[in] RAAN Right ascension of the ascending node at launch date [rad].
# @param[in] numDays Number of days in the analysis.
#
# @return The beta angle computed for each day in degrees.
#
#==#

function satellite_eclipse_time(t0::Integer,
                                a::FloatingPoint,
                                e::FloatingPoint,
                                i::FloatingPoint,
                                w::FloatingPoint,
                                RAAN::FloatingPoint,
                                numDays::Integer)
    # Constants
    const deg2rad = pi/180.0
    const rad2deg = 180.0/pi
    const day2sec = 24.0*60.0*60.0
    
    # Step of the orbit propagation (mean anomaly) [rad].
    const step = 0.1*deg2rad
    
    # Initialization of variables.
    theta = 0.0            # Sun angle relative to the inertial coordinate frame.
    
    days = [0:1:numDays-1] # Vector of the days in which the beta angle will be
                           # computed.

    # Mean anomaly.
    M = [0:step:2*pi]

    # Penumbra time.
    p_time = zeros(numDays)

    # Umbra time.
    u_time = zeros(numDays)

    # Sunlight time.
    s_time = zeros(numDays)
    
    # Semi-lactum rectum.
    p = a*(1.0-e^2)

    # Period of an orbit [rad/s].
    n = n_J2(a, e, i)

    # Step in time
    tstep = step/n

    # Perturbations.
    #
    # RAAN rotation rate [rad/s].
    dOmega = dRAAN_J2(a, e, i)

    # Perturbation of the argument of perigee [rad/s].
    dw = dw_J2(a, e, i)

    # Loop.
    for d in days
        # Get the sun position represented in the Inertial coordinate frame.
        s_i = sun_position_i(int(t0+d), 43200)

        # Compute the new orbit parameters due to perturbations.
        w_d = w + dw*(d*day2sec)
        RAAN_d = RAAN + dOmega*(d*day2sec)

        for m_k in M
            # Get the satellite position vector represented in the Inertial coordinate
            # frame.
            f = satellite_orbit_compute_f(a, e, i, m_k)

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

    [days s_time p_time u_time]
end