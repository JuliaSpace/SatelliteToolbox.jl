#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export beta_angle

"""
    function beta_angle(JD₀::Number, a::Number, e::Number, i::Number, RAAN::Number, Δt::Integer)

Compute the beta angle of an orbit with semi-major axis `a` [m], excentricity
`e`, inclination `i` [rad], and initial right ascension of the ascending node
`RAAN` [rad]. The orbit epoch, which is also the day in which the analysis will
begin, is `JD₀` [Julian Day]. The analysis will be performed for the number of
days in `Δt`.

Notice that the simulation considers only the perturbation terms up to J2.

# Returns

An array containing the beta angle [rad] for each day in the analysis.

"""
function beta_angle(JD₀::Number, a::Number, e::Number, i::Number, RAAN::Number,
                    Δt::Integer)
    # Constants
    rad2deg = 180.0/pi

    # Initialization of variables.
    theta = 0.0                   # Sun angle relative to the inertial
                                  # coordinate frame.

    days = collect(0:1:Δt-1)      # Vector of the days in which the beta angle
                                  # will be computed.

    N = [0.0; 0.0; 0.0]           # Versor normal to the orbital plane,
                                  # represented in the inertial coordinate
                                  # frame.

    S = [0.0; 0.0; 0.0]           # Versor that points to the Sun, represented
                                  # in the inertial coordinate frame.

    # Output vector.
    beta = Vector{Float64}(undef,Δt)

    # RAAN rotation rate [rad/day].
    dOmega = dRAAN(a, e, i, :J2)*24.0*3600.0

    # Loop
    for t in days
        # Compute the RAAN at the day d.
        RAAN_d = RAAN + dOmega*t

        # Compute the versor N represented in the Inertial ref. frame.
        Dio = angle_to_dcm(-i, -RAAN_d, 0.0, :XZX)
        N_i = Dio*SVector{3}(0,0,1)

        # Compute the Sun position at noon (UT) represented in the Inertial ref.
        # frame.
        S_i = sun_position_i(JD₀+t)
        S_i = S_i/norm(S_i)

        # Get the angle between N_i and S_i [deg].
        beta[t+1] = 90.0-acos(dot(N_i,S_i))*rad2deg
    end

    [days beta]
end
