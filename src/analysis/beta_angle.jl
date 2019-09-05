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
    function beta_angle(JD₀::Number, a::Number, e::Number, i::Number, RAAN::Number, Δt::Integer, pert::Symbol = :J2)

Compute the beta angle of an orbit with semi-major axis `a` [m], eccentricity
`e`, inclination `i` [rad], and initial right ascension of the ascending node
`RAAN` [rad]. The orbit epoch, which is also the day in which the analysis will
begin, is `JD₀` [Julian Day]. The analysis will be performed for each day during
`Δt` days.

The argument `pert` can be used to select the perturbation terms that must be
used when propagating the right ascencion of the ascending node. The possible
values are:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.
* `:J4`: Consider the perturbation terms J2, J4, and J2².

If `pert` is omitted, then it defaults to `:J2`.

# Returns

An array with two columns. The first one contains the days of the analysis and
the second one contains the beta angle [rad] for each day.

"""
function beta_angle(JD₀::Number, a::Number, e::Number, i::Number, RAAN::Number,
                    Δt::Integer, pert::Symbol = :J2)

    # Vector of the days in which the beta angle will be computed.
    days = collect(0:1:Δt-1)

    # Output vector.
    β = Vector{Float64}(undef,Δt)

    # RAAN rotation rate [rad/day].
    δΩ = 86400dRAAN(a, e, i, pert)

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
