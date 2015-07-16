#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
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
#    Compute the sun radiaton received in a satellite surface.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2015-05-20: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version, based on satellite_sun_angle.jl.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

import Rotations: angle2dcm!

export satellite_sun_radiation_earth_pointing

#==#
# 
# @brief Compute the sun radiation on a surface for an Earth pointing mission.
#
# @param[in] t0 Launch date [number of days since 01/01/2000].
# @param[in] a Semi-major axis of the satellite orbit [m].
# @param[in] e Orbit eccentricity.
# @param[in] i Orbit inclination [rad].
# @param[in] w Argument of perigee [rad].
# @param[in] RAAN Right ascension of the ascending node at launch date [rad].
# @param[in] numDays Number of days in the analysis.
# @param[in] fN_k Function f(s_b) that describes the solar panel normal at each
# k-th sampling step. s_b is the sun vector represented in the body coordinate
# frame.
# @param[in] step (OPTIONAL) Mean anomaly step (default = 0.1 deg).
#
# @return A matrix containing the sun radiation [W/m^2] for each position in
# orbit for each day. If the sun angle is larger than 90 deg or if the satellite
# is in eclipse, then NaN is returned in the matrix.
#
# @note The body reference frame is defined as:
#     _ Z axis points towards the center of Earth;
#     _ Y axis points towards the negative direction of orbit normal;
#     _ X axis completes the right-hand reference frame.
# which is common for Earth pointing satellites.
#
#==#

function satellite_sun_radiation_earth_pointing(t0::Integer,
                                                a::FloatingPoint,
                                                e::FloatingPoint,
                                                i::FloatingPoint,
                                                RAAN::FloatingPoint,
                                                w::FloatingPoint,
                                                numDays::Integer,
                                                fN_k::Function,
                                                step::Float64 = 0.1*pi/180.0)
    # Constants
    const deg2rad = pi/180.0
    const rad2deg = 180.0/pi
    const day2sec = 24.0*60.0*60.0

    # Initialization of variables.
    theta = 0.0            # Sun angle relative to the inertial coordinate frame.
    
    days = [0:1:numDays-1] # Vector of the days in which the beta angle will be
                           # computed.

    # Mean anomaly.
    M = [0:step:2*pi]

    # Period of an orbit [rad/s].
    n = n_J2(a, e, i)

    # Step in time
    tstep = step/n

    # Sun angles.
    sun_radiation = zeros(length(M),numDays)
    
    # Perturbations.
    #
    # RAAN rotation rate [rad/s].
    dOmega = dRAAN_J2(a, e, i)

    # Perturbation of the argument of perigee [rad/s].
    dw = dw_J2(a, e, i)
    
    # DCM that rotates the Inertial reference frame to the orbit reference frame.
    Doi = Array(Float64, (3,3))

    # DCM that rotates the orbit reference frame to the body reference frame.
    #
    # In this case, the body reference frame is defined as:
    #     _ Z axis points towards the center of Earth;
    #     _ Y axis points towards the negative direction of orbit normal;
    #     _ X axis completes the right-hand reference frame.
    # which is common for Earth pointing satellites.
    
    Dbo = [ 0.0 1.0  0.0;
            0.0 0.0 -1.0;
           -1.0 0.0  0.0];
    
    # Loop for each day.
    for d in days
        # Get the sun position represented in the Inertial coordinate frame.
        s_i = sun_position_i(int(t0+d), 43200)
        norm_s_i = norm(s_i)
        
        # Compute the new orbit parameters due to perturbations.
        w_d    = w + dw*(d*day2sec)
        RAAN_d = RAAN + dOmega*(d*day2sec)
        
        # Loop through the orbit.
        for k in [1:length(M)]
            # Get the satellite position vector represented in the Inertial
            # coordinate frame.
            f = satellite_orbit_compute_f(a, e, i, M[k])

            (r_i, rt_i) = satellite_position_i(a, e, i, RAAN_d, w_d, f)
            
            # Check the lighting conditions.
            lighting = satellite_lighting_condition(r_i, s_i)

            if (lighting == SAT_LIGHTING_SUNLIGHT)
                # Convert the sun vector from the Inertial coordinate frame to
                # the body coordinate frame.
                angle2dcm!(Doi, RAAN_d, i, w_d+f, "ZXZ")
                s_b = Dbo*Doi*(s_i/norm_s_i)

                # Vector normal to the solar panel.
                N_k = fN_k(s_b)
                sun_angle_k = acos(dot(s_b,N_k))
                
                # If the sun angle is larger than 90 deg, then the surface is
                # not iluminated. Thus, the angle will be defined as NaN.
                if (sun_angle_k > pi/2)
                    sun_radiation[k, d+1] = NaN
                else
                    sun_radiation[k, d+1] =
                        sunRad/(4*pi*norm_s_i^2)*cos(sun_angle_k)
                end
            else
                # If the satellite is in eclipse, then the surface is not
                # iluminated. Thus, the angle will be defined as NaN.
                sun_radiation[k,d+1] = NaN 
            end
        end 
    end

    sun_radiation
end

#==#
# 
# @brief Compute the sun radiation on a surface for an Earth pointing mission.
#
# @param[in] t0 Launch date [number of days since 01/01/2000].
# @param[in] a Semi-major axis of the satellite orbit [m].
# @param[in] e Orbit eccentricity.
# @param[in] i Orbit inclination [rad].
# @param[in] w Argument of perigee [rad].
# @param[in] RAAN Right ascension of the ascending node at launch date [rad].
# @param[in] numDays Number of days in the analysis.
# @param[in] N Vector normal to the surface represented in the body reference
# frame.
# @param[in] step (OPTIONAL) Mean anomaly step (default = 0.1 deg).
#
# @return A matrix containing the sun radiation [W/m^2] for each position in
# orbit for each day. If the sun angle is larger than 90 deg or if the satellite
# is in eclipse, then NaN is returned in the matrix.
#
# @note The body reference frame is defined as:
#     _ Z axis points towards the center of Earth;
#     _ Y axis points towards the negative direction of orbit normal;
#     _ X axis completes the right-hand reference frame.
# which is common for Earth pointing satellites.
#
#==#

function satellite_sun_radiation_earth_pointing(t0::Integer,
                                                a::FloatingPoint,
                                                e::FloatingPoint,
                                                i::FloatingPoint,
                                                RAAN::FloatingPoint,
                                                w::FloatingPoint,
                                                numDays::Integer,
                                                N::Array{Float64,1},
                                                step::Float64 = 0.1*pi/180.0)
    fN_k(x) = N
    satellite_sun_radiation_earth_pointing(t0, a, e, i, RAAN, w, numDays, fN_k, step)
end
