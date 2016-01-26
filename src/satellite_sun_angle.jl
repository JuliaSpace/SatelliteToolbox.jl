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
#    Compute the sun angle on a satellite surface.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2015-03-05: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#
#    Add function to compute the sun angle depending on a user-provided function
#    that describe the surface normal according to the solar vector represented
#    in the body coordinate frame. Thus, if the solar panel rotates about the
#    body Y axis to maximize the sun angle, then the function must be:
#
#        function fN_k(s_b)
#            theta_p = atan2(-s_b[3], s_b[1])
#            [cos(theta_p); 0.0; -sin(theta_p)]
#        end
#
# 2015-03-05: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

import Rotations: angle2dcm!

export satellite_sun_angle_earth_pointing

"""
### function satellite_sun_angle_earth_pointing(t0::Integer, a::Real, e::Real, i::Real, RAAN::Real, w::Real, numDays::Integer, fN_k::Function, step::Float64 = 0.1*pi/180.0)

Compute the Sun angle on a surface for an Earth-pointing mission.

##### Args

* t0: Launch date [number of days since 01/01/2000].
* a: Semi-major axis of the satellite orbit [m].
* e: Orbit eccentricity.
* i: Orbit inclination [rad].
* w: Argument of perigee [rad].
* RAAN: Right ascension of the ascending node at launch date [rad].
* numDays: Number of days in the analysis.
* fN_k: Function **f(s_b)** that describes the solar panel normal at each k-th
sampling step. Notice that **s_b** is the Sun vector represented in the body
coordinate frame.
* step: (OPTIONAL) Mean anomaly step, *default*: 0.1 deg.

##### Returns

* A matrix containing the sun angle for each position in orbit for each day.

**NOTE**: if the sun angle is larger than 90 deg or if the satellite is in
eclipse, then NaN is returned in the matrix.

##### Remarks

The body reference frame is defined as:

* **Z axis** points towards the center of Earth;
* **Y axis** points towards the negative direction of orbit normal;
* **X axis** completes the right-hand reference frame.

"""

function satellite_sun_angle_earth_pointing(t0::Integer,
                                            a::Real,
                                            e::Real,
                                            i::Real,
                                            RAAN::Real,
                                            w::Real,
                                            numDays::Integer,
                                            fN_k::Function,
                                            step::Float64 = 0.1*pi/180.0)
    # Constants
    const deg2rad = pi/180.0
    const rad2deg = 180.0/pi
    const day2sec = 24.0*60.0*60.0

    # Initialization of variables.
    theta = 0.0                   # Sun angle relative to the inertial
                                  # coordinate frame.

    days = collect(0:1:numDays-1) # Vector of the days in which the sun angle
                                  # will be computed.

    # Mean anomaly.
    M = collect(0:step:2*pi)

    # Period of an orbit [rad/s].
    n = n_J2(a, e, i)

    # Step in time
    tstep = step/n

    # Sun angles.
    sun_angles = zeros(length(M),numDays)

    # Perturbations.
    #
    # RAAN rotation rate [rad/s].
    dOmega = dRAAN_J2(a, e, i)

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
        s_i = sun_position_i(Int(t0+d), 43200)
        norm_s_i = norm(s_i)

        # Compute the new orbit parameters due to perturbations.
        RAAN_d = RAAN + dOmega*(d*day2sec)

        # Loop through the orbit.
        for k in 1:length(M)
            # Get the satellite position vector represented in the Inertial
            # coordinate frame.
            f = satellite_orbit_compute_f(a, e, i, M[k])

            (r_i, rt_i) = satellite_position_i(a, e, i, RAAN_d, w, f)

            # Check the lighting conditions.
            lighting = satellite_lighting_condition(r_i, s_i)

            if (lighting == SAT_LIGHTING_SUNLIGHT)
                # Convert the sun vector from the Inertial coordinate frame to
                # the body coordinate frame.
                angle2dcm!(Doi, RAAN_d, i, f, "ZXZ")
                s_b = Dbo*Doi*(s_i/norm_s_i)

                # Vector normal to the solar panel.
                N_k = fN_k(s_b)
                sun_angle_k = acos(dot(s_b,N_k))

                # If the sun angle is larger than 90 deg, then the surface is
                # not iluminated. Thus, the angle will be defined as NaN.
                if (sun_angle_k > pi/2)
                    sun_angles[k, d+1] = NaN
                else
                    sun_angles[k, d+1] = sun_angle_k
                end
            else
                # If the satellite is in eclipse, then the surface is not
                # iluminated. Thus, the angle will be defined as NaN.
                sun_angles[k,d+1] = NaN
            end
        end
    end

    sun_angles
end

"""
### function satellite_sun_angle_earth_pointing(t0::Integer, a::Real, e::Real, i::Real, RAAN::Real, w::Real, numDays::Integer, N::Array{Float64,1}, step::Float64 = 0.1*pi/180.0)

Compute the Sun angle on a surface for an Earth-pointing mission.

##### Args

* t0: Launch date [number of days since 01/01/2000].
* a: Semi-major axis of the satellite orbit [m].
* e: Orbit eccentricity.
* i: Orbit inclination [rad].
* w: Argument of perigee [rad].
* RAAN: Right ascension of the ascending node at launch date [rad].
* numDays: Number of days in the analysis.
* N: Vector normal to the surface represented in the body reference frame.
* step: (OPTIONAL) Mean anomaly step, *default*: 0.1 deg.

##### Returns

* A matrix containing the sun angle for each position in orbit for each day.

**NOTE**: if the sun angle is larger than 90 deg or if the satellite is in
eclipse, then NaN is returned in the matrix.

##### Remarks

The body reference frame is defined as:

* **Z axis** points towards the center of Earth;
* **Y axis** points towards the negative direction of orbit normal;
* **X axis** completes the right-hand reference frame.

"""

function satellite_sun_angle_earth_pointing(t0::Integer,
                                            a::Real,
                                            e::Real,
                                            i::Real,
                                            RAAN::Real,
                                            w::Real,
                                            numDays::Integer,
                                            N::Array{Float64,1},
                                            step::Float64 = 0.1*pi/180.0)
    fN_k(x) = N
    satellite_sun_angle_earth_pointing(t0, a, e, i, RAAN, w, numDays, fN_k, step)
end
