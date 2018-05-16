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
#   Compute the beta angle of a satellite.
#
#   The algorithm was based on:
#
#       Mortari, D., Wilkins, M. P., and Bruccoleri, C.  On Sun-Synchronous
#       Orbits and Associated Constellations
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2015-01-20: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export satellite_beta_angle

"""
    function satellite_beta_angle(JD0::Number, a::Number, e::Number, i::Number, RAAN::Number, numDays::Integer)

Compute the beta angle of a satellite.

##### Args

* `JD0`: Initial instant for the analysis [Julian day].
* `a`: Semi-major axis of the orbit [m].
* `e`: Orbit eccentricity.
* `i`: Orbit inclination [rad].
* `RAAN`: Right ascension of the ascending node at `JD0` [rad].
* `numDays`: Number of days of the analysis.

##### Returns

The beta angle computed for each day in degrees.

"""
function satellite_beta_angle(JD0::Number,
                              a::Number,
                              e::Number,
                              i::Number,
                              RAAN::Number,
                              numDays::Integer)
    # Constants
    rad2deg = 180.0/pi

    # Initialization of variables.
    theta = 0.0                   # Sun angle relative to the inertial
                                  # coordinate frame.

    days = collect(0:1:numDays-1) # Vector of the days in which the beta angle
                                  # will be computed.

    N = [0.0; 0.0; 0.0]           # Versor normal to the orbital plane,
                                  # represented in the inertial coordinate
                                  # frame.

    S = [0.0; 0.0; 0.0]           # Versor that points to the Sun, represented
                                  # in the inertial coordinate frame.

    # Output vector.
    beta = Array{Float64}(numDays,1)

    # RAAN rotation rate [rad/day].
    dOmega = dRAAN(a, e, i, :J2)*24.0*3600.0

    # Loop
    for t in days
        # Compute the RAAN at the day d.
        RAAN_d = RAAN + dOmega*t

        # Compute the versor N represented in the Inertial ref. frame.
        Dio = angle2dcm(-i, -RAAN_d, 0.0, :XZX)
        N_i = Dio*SVector{3}(0,0,1)

        # Compute the Sun position at noon (UT) represented in the Inertial ref.
        # frame.
        S_i = sun_position_i(JD0+t)
        S_i = S_i/norm(S_i)

        # Get the angle between N_i and S_i [deg].
        beta[t+1] = 90.0-acos(dot(N_i,S_i))*rad2deg
    end

    [days beta]
end
