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
#    Compute the beta angle of a satellite.
#
#    The algorithm was based on:
#
#          Mortari, D., Wilkins, M. P., and Bruccoleri, C.
#          On Sun-Synchronous Orbits and Associated Constellations
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2015-01-20: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

import Rotations: angle2dcm!

export satellite_beta_angle

"""
### function satellite_beta_angle(t0::Integer, a::Real, e::Real, i::Real, RAAN::Real, numDays::Integer)

Compute the beta angle of a satellite.

##### Args

* t0: Launch date [number of days since 01/01/2000].
* a: Semi-major axis of the satellite orbit [m].
* e: Orbit eccentricity.
* i: Orbit inclination [rad].
* RAAN: Right ascension of the ascending node at launch date [rad].
* numDays: Number of days of the analysis.

##### Returns

* The beta angle computed for each day in degrees.

"""

function satellite_beta_angle(t0::Integer,
                              a::Real,
                              e::Real,
                              i::Real,
                              RAAN::Real,
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
    beta = Array(Float64, (numDays,1))

    # RAAN rotation rate [rad/day].
    dOmega = dRAAN_J2(a, e, i)*24.0*3600.0

    # DCM that rotates the orbit reference frame to the Inertial reference frame.
    Dio = Array(Float64, (3,3))

    # Loop
    for t in days
        # Compute the RAAN at the day d.
        RAAN_d = RAAN + dOmega*t

        # Compute the versor N represented in the Inertial ref. frame.
        angle2dcm!(Dio, -i, -RAAN_d, 0.0, "XZX")
        N_i = Dio*[0;0;1]

        # Compute the Sun position at noon (UT) represented in the Inertial ref.
        # frame.
        S_i = sun_position_i(Int(t+t0), 43200)
        S_i = S_i/norm(S_i)

        # Get the angle between N_i and S_i [deg].
        beta[t+1] = 90.0-acos(dot(N_i,S_i))*rad2deg
    end

    [days beta]
end
