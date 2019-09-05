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

Compute the beta angle of an orbit with semi-major axis `a` [m], eccentricity
`e`, inclination `i` [rad], and initial right ascension of the ascending node
`RAAN` [rad]. The orbit epoch, which is also the day in which the analysis will
begin, is `JD₀` [Julian Day]. The analysis will be performed for each day during
`Δt` days.

Notice that the simulation considers only the perturbation terms up to J2.

# Returns

An array with two columns. The first one contains the days of the analysis and
the second one contains the beta angle [rad] for each day.

"""
function beta_angle(JD₀::Number, a::Number, e::Number, i::Number, RAAN::Number,
                    Δt::Integer)

    # Vector of the days in which the beta angle will be computed.
    days = collect(0:1:Δt-1)

    # Output vector.
    β = Vector{Float64}(undef,Δt)

    # RAAN rotation rate [rad/day].
    δΩ = 86400dRAAN(a, e, i, :J2)

    # Loop
    @inbounds @views for d in days
        # Compute the RAAN at the day d.
        Ω = RAAN + δΩ*d

        # Compute the versor N represented in the Inertial ref. frame.
        Dio = angle_to_dcm(-i, -Ω, 0, :XZX)
        N_i = Dio[:,3]  # -> Dio*[0;0;1]

        # Compute the Sun position at noon (UT) represented in the Inertial ref.
        # frame.
        S_i = sun_position_i(JD₀+d)
        S_i = S_i/norm(S_i)

        # Compute the beta angle, which is the angle between the Sun vector and
        # the orbit plane.
        β[d+1] = 90 - acosd(dot(N_i,S_i))
    end

    return [days β]
end
