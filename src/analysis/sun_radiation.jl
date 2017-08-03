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
#    Compute the sun radiation received in a satellite surface.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2015-07-16: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    The function to compute the sun radiation can propagate the orbit using the
#    orbit latitude or the mean anomaly. Notice that if the latter is used, then
#    the average value of the output will be the average sun radiation that the
#    surface will receive, because the angular steps have the same time
#    interval.
#
#    Add functions to compute the sun radiation average per day.
#
# 2015-05-20: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version, based on satellite_sun_angle.jl.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

import Rotations: angle2dcm!

export satellite_sun_radiation_earth_pointing
export satellite_sun_radiation_earth_pointing_mean

"""
### function satellite_sun_radiation_earth_pointing(JD0::Real, a::Real, e::Real, i::Real, RAAN::Real, w::Real, numDays::Integer, fN_k::Function, meanAnomaly::Bool = false, step::Float64 = 0.1*pi/180.0)

Compute the sun radiation on a surface for an Earth pointing mission.

##### Args

* JD0: Julian day of the launch date.
* a: Semi-major axis of the satellite orbit [m].
* e: Orbit eccentricity.
* i: Orbit inclination [rad].
* w: Argument of perigee [rad].
* RAAN: Right ascension of the ascending node at launch date [rad].
* numDays: Number of days in the analysis.
* fN_k: Function **f(s_b)** that describes the solar panel normal at each k-th
sampling step. Notice that **s_b** is the Sun vector represented in the body
coordinate frame.
* meanAnomaly: (OPTIONAL) If **true**, compute using angular steps in the mean
anomaly instead of in the orbit latitude, *default*: **false**.
* step: (OPTIONAL) Mean anomaly step, *default*: 0.1 deg.

##### Returns

* A matrix containing the Sun radiation [W/m²] for each position in orbit for
each day.

**NOTE**: if the sun angle is larger than 90 deg or if the satellite is in
eclipse, then NaN is returned in the matrix.

##### Remarks

The body reference frame is defined as:

* **Z axis** points towards the center of Earth;
* **Y axis** points towards the negative direction of orbit normal;
* **X axis** completes the right-hand reference frame.

If the **mean anomaly** is used, then the average value of the output is the
average sun radiation received by the satellite surface, because every angular
steps have a fixed time interval.

If the **mean anomaly** is used, then the angle interval is [0, 2π]
Otherwise, the angle interval is [-π,π].

"""

function satellite_sun_radiation_earth_pointing(JD0::Real,
                                                a::Real,
                                                e::Real,
                                                i::Real,
                                                RAAN::Real,
                                                w::Real,
                                                numDays::Integer,
                                                fN_k::Function,
                                                meanAnomaly::Bool = false,
                                                step::Float64 = 0.1*pi/180.0)
    # Constants
    const deg2rad = pi/180.0
    const rad2deg = 180.0/pi
    const day2sec = 24.0*60.0*60.0

    # Initialization of variables.
    theta = 0.0                   # Sun angle relative to the inertial
                                  # coordinate frame.

    days = collect(0:1:numDays-1) # Vector of the days in which the eclipse time
                                  # will be computed.

    # Angle.
    ang = (!meanAnomaly) ? collect(-pi:step:pi) : collect(0:step:2*pi)

    # Period of an orbit [rad/s].
    n = n_J2(a, e, i)

    # Step in time
    tstep = step/n

    # Sun radiation.
    sun_radiation = zeros(length(ang),numDays)

    # Perturbations.
    #
    # RAAN rotation rate [rad/s].
    dOmega = dRAAN_J2(a, e, i)

    # Perturbation of the argument of perigee [rad/s].
    dw = dw_J2(a, e, i)

    # DCM that rotates the Inertial reference frame to the orbit reference frame.
    Doi = Array{Float64}(3,3)

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
        s_i = sun_position_i(JD0+d)
        norm_s_i = norm(s_i)

        # Compute the new orbit parameters due to perturbations.
        w_d    = w + dw*(d*day2sec)
        RAAN_d = RAAN + dOmega*(d*day2sec)

        # Loop through the orbit.
        for k in 1:length(ang)
            # Get the satellite position vector represented in the Inertial
            # coordinate frame.

            if (!meanAnomaly)
                f = ang[k]-w_d
            else
                f = satellite_orbit_compute_f(a, e, i, ang[k])
            end

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

                # Normalize N_k.
                N_k = N_k/norm(N_k)

                # Compute the sun angle.
                sun_angle_k = acos(dot(s_b,N_k))

                # If the sun angle is larger than 90 deg, then the surface is
                # not illuminated. Thus, the angle will be defined as NaN.
                if (sun_angle_k > pi/2)
                    sun_radiation[k, d+1] = NaN
                else
                    sun_radiation[k, d+1] =
                        sunRad/(4*pi*norm_s_i^2)*cos(sun_angle_k)
                end
            else
                # If the satellite is in eclipse, then the surface is not
                # illuminated. Thus, the angle will be defined as NaN.
                sun_radiation[k,d+1] = NaN
            end
        end
    end

    sun_radiation
end

"""
### function satellite_sun_radiation_earth_pointing_mean(JD0::Real, a::Real, e::Real, i::Real, RAAN::Real, w::Real, numDays::Integer, fN_k::Function, step::Float64 = 0.1*pi/180.0)

Compute the mean sun radiation on a surface for an Earth-pointing mission.

##### Args

* JD0: Julian day of the launch date.
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

* The mean Sun radiation on a surface [W/m²].

##### Remarks

For more details, *see* **satellite_sun_radiation_earth_pointing**.

"""

function satellite_sun_radiation_earth_pointing_mean(JD0::Real,
                                                     a::Real,
                                                     e::Real,
                                                     i::Real,
                                                     RAAN::Real,
                                                     w::Real,
                                                     numDays::Integer,
                                                     fN_k::Function,
                                                     step::Float64 = 0.1*pi/180.0)
    sun_radiation = satellite_sun_radiation_earth_pointing(JD0,
                                                           a,
                                                           e,
                                                           i,
                                                           RAAN,
                                                           w,
                                                           numDays,
                                                           fN_k,
                                                           true,
                                                           step)

    # Remove the NaNs.
    sun_radiation[find(x->(isnan(x)), sun_radiation)] = 0.0

    # Compute the mean.
    squeeze(mean(sun_radiation, 1), 1)
end

"""
### function satellite_sun_radiation_earth_pointing(JD0::Real, a::Real, e::Real, i::Real, RAAN::Real, w::Real, numDays::Integer, N::Array{Float64,1}, meanAnomaly::Bool = false, step::Float64 = 0.1*pi/180.0)

Compute the sun radiation on a surface for an Earth pointing mission.

##### Args

* JD0: Julian day of the launch date.
* a: Semi-major axis of the satellite orbit [m].
* e: Orbit eccentricity.
* i: Orbit inclination [rad].
* w: Argument of perigee [rad].
* RAAN: Right ascension of the ascending node at launch date [rad].
* numDays: Number of days in the analysis.
* N: Vector normal to the surface represented in the body reference frame.
* meanAnomaly: (OPTIONAL) If **true**, compute using angular steps in the mean
anomaly instead of in the orbit latitude, *default*: **false**.
* step: (OPTIONAL) Mean anomaly step, *default*: 0.1 deg.

##### Returns

* A matrix containing the Sun radiation [W/m²] for each position in orbit for
each day.

**NOTE**: if the sun angle is larger than 90 deg or if the satellite is in
eclipse, then NaN is returned in the matrix.

##### Remarks

The body reference frame is defined as:

* **Z axis** points towards the center of Earth;
* **Y axis** points towards the negative direction of orbit normal;
* **X axis** completes the right-hand reference frame.

If the **mean anomaly** is used, then the average value of the output is the
average sun radiation received by the satellite surface, because every angular
steps have a fixed time interval.

If the **mean anomaly** is used, then the angle interval is [0, 2π]
Otherwise, the angle interval is [-π,π].

"""

function satellite_sun_radiation_earth_pointing(JD0::Real,
                                                a::Real,
                                                e::Real,
                                                i::Real,
                                                RAAN::Real,
                                                w::Real,
                                                numDays::Integer,
                                                N::Array{Float64,1},
                                                meanAnomaly::Bool = false,
                                                step::Float64 = 0.1*pi/180.0)
    fN_k(x) = N
    satellite_sun_radiation_earth_pointing(JD0, a, e, i, RAAN, w, numDays, fN_k,
                                           meanAnomaly, step)
end

"""
### function satellite_sun_radiation_earth_pointing_mean(JD0::Real, a::Real, e::Real, i::Real, RAAN::Real, w::Real, numDays::Integer, N::Array{Float64,1}, step::Float64 = 0.1*pi/180.0)

Compute the mean sun radiation on a surface for an Earth-pointing mission.

##### Args

* JD0: Julian day of the launch date.
* a: Semi-major axis of the satellite orbit [m].
* e: Orbit eccentricity.
* i: Orbit inclination [rad].
* w: Argument of perigee [rad].
* RAAN: Right ascension of the ascending node at launch date [rad].
* numDays: Number of days in the analysis.
* N: Vector normal to the surface represented in the body reference frame.
* step: (OPTIONAL) Mean anomaly step, *default*: 0.1 deg.

##### Returns

* The mean Sun radiation on a surface [W/m²].

##### Remarks

For more details, *see* **satellite_sun_radiation_earth_pointing**.

"""

function satellite_sun_radiation_earth_pointing_mean(JD0::Real,
                                                     a::Real,
                                                     e::Real,
                                                     i::Real,
                                                     RAAN::Real,
                                                     w::Real,
                                                     numDays::Integer,
                                                     N::Array{Float64,1},
                                                     step::Float64 = 0.1*pi/180.0)
    sun_radiation = satellite_sun_radiation_earth_pointing(JD0,
                                                           a,
                                                           e,
                                                           i,
                                                           RAAN,
                                                           w,
                                                           numDays,
                                                           N,
                                                           true,
                                                           step)

    # Remove the NaNs.
    sun_radiation[find(x->(isnan(x)), sun_radiation)] = 0.0

    # Compute the mean.
    squeeze(mean(sun_radiation, 1), 1)
end
