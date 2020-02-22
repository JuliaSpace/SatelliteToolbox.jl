#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Many auxiliary functions for ground repeating, Sun-synchronous (GRSS) orbit
#   computations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export adjacent_track_distance_grss, adjacent_track_angle_grss
export compute_ss_orbit_by_num_rev_per_day
export list_ss_orbits_by_rep_period, sort_list_ss_orbits_by_height

"""
    adjacent_track_distance_grss(T::Number, i::Number, To::Int, lat::Number)

Compute the distance between adjacent ground tracks [m] at a given latitude
`lat` [rad] for a ground repeating, Sun-synchronous orbit with period `T` [s],
inclination `i` [rad], and orbit cycle `To` [days].

# Remarks

The functions **does not** check if the orbit is a GRSS orbit.

"""
function adjacent_track_distance_grss(T::Number, i::Number, To::Int, lat::Number)
    # Angle between two adjacent traces.
    theta = (T/To)*pi/43200.0

    # Angle of swath.
    beta = asin(sin(theta)*sin(i))

    # Minimum swath.
    beta*R0*cos(lat)
end

"""
    adjacent_track_distance_grss(a::Number, e::Number, i::Number, To::Int, lat::Number)

Compute the distance between adjacent ground tracks [m] at a given latitude
`lat` [rad] for a ground repeating, Sun-synchronous orbit with semi-major axis
`a` [m], eccentricity `e`, inclination `i` [rad], and orbit cycle `To` [days].

# Remarks

The functions *does not* check if the orbit is a GRSS orbit.

"""
function adjacent_track_distance_grss(a::Number, e::Number, i::Number, To::Int,
                                      lat::Number)
    # Orbit period.
    T = period(a, e, i, :J2)

    # Compute the distance between adjacent ground tracks.
    adjacent_track_distance_grss(T, i, To, lat)
end

"""
    adjacent_track_angle_grss(h::Number, T::Number, i::Number, To::Int, lat::Number)

Compute the angle between two adjacent ground tracks [rad] in a given latitude
`lat` [rad] measured from the satellite position for a ground repeating,
Sun-synchronous orbit with altitude in the Equator `h` [m], period `T` [s],
inclination `i` [rad], and orbit cycle `To` [days].

# Remarks

The functions **does not** check if the orbit is a GRSS orbit.

"""
function adjacent_track_angle_grss(h::Number, T::Number, i::Number, To::Int,
                                   lat::Number)
    # Angle between two adjacent traces.
    theta = (T/To)*pi/43200.0

    # Angle of swath.
    beta = asin(sin(theta)*sin(i))

    # Get the Earth radius in the plane perpendicular to the Equatorial plane.
    Rp = R0*cos(lat)

    # Compute the angle between the two ground tracks.
    b = sqrt(Rp^2+(Rp+h)^2-2*Rp*(Rp+h)*cos(beta/2.0))
    asin(Rp/b*sin(beta/2.0))
end

"""
    adjacent_track_angle_grss(h::Number, a::Number, e::Number, i::Number, To::Int, lat::Number)

Compute the angle between two adjacent ground tracks [rad] in a given latitude
`lat` [rad] measured from the satellite position for a ground repeating,
Sun-synchronous orbit with altitude in the Equator `h` [m], semi-major axis `a`
[m], eccentricity `e`, inclination `i` [rad], and orbit cycle `To` [days].

# Remarks

The functions *does not* check if the orbit is a GRSS orbit.

"""
function adjacent_track_angle_grss(h::Number, a::Number, e::Number, i::Number,
                                   To::Int, lat::Number)
    # Period (J2).
    T = period(a, e, i, :J2)

    # Compute the minimum half FOV.
    adjacent_track_angle_grss(h, T, i, To, lat)
end

"""
    compute_ss_orbit_by_num_rev_per_day(numRevPD::Number, e::Number)

Compute the Sun-synchronous orbit given the number of revolutions per day
`numRevPD` and the eccentricity `e`.

# Returns

* The semi-major axis [m].
* The inclination [rad].
* The residues of the two functions.
* A boolean variable that indicates if the numerical algorithm converged.

"""
function compute_ss_orbit_by_num_rev_per_day(numRevPD::Number, e::Number)
    compute_ss_orbit_by_ang_vel(numRevPD*2*pi/86400.0, e)
end

"""
    list_ss_orbits_by_rep_period(minRep::Int, maxRep::Int, minAlt::Number=-1.0, maxAlt::Number=-1.0, e::Number=0.0)

Compute a list of repeating Sun-synchronous orbits.

# Args

* `minRep`: Minimum repetition time of the orbit [days].
* `maxRep`: Maximum repetition time of the orbit [days].
* `minAlt`: Minimum altitude of the orbits on the list [m].
* `maxAlt`: Minimum altitude of the orbits on the list [m].
* `e`: Eccentricity.

# Returns

A matrix containing the orbits found with the following format:

    Semi-major axis [m] | Altitude [m] | Inclination [rad] | Period [s] | Int | Num | Den
    --------------------|--------------|-------------------|------------|-----|-----|-----

in which the period is Int + Num/Den.

# Remarks

If `minAlt` or `maxAlt` is < 0.0, then the altitude will not be checked when a
orbit is added to the list.

"""
function list_ss_orbits_by_rep_period(minRep::Int,         maxRep::Int,
                                      minAlt::Number=-1.0, maxAlt::Number=-1.0,
                                      e::Number=0.0)
    # Check if the arguments are valid.
    minRep <= 0       && throw(ArgumentError("The minimum repetition time must be greater than 0."))
    maxRep <= 0       && throw(ArgumentError("The maximum repetition time must be greater than 0."))
    minRep > maxRep   && throw(ArgumentError("The minimum repetition time must be smaller or equal to the maximum repetition time."))
    !( 0. <= e < 1. ) && throw(ArgumentError("The eccentricity must be within the interval 0 <= e < 1."))

    # Matrix to store the available orbits.
    ss_orbits = Matrix{Float64}(undef, 0, 7)

    # Integer part of the number of orbits per day.
    intNumOrb = (13, 14, 15, 16, 17)

    # Check if the altitude interval must be verified.
    check_altitude = (minAlt > 0) && (maxAlt > 0)
    check_altitude && (minAlt > maxAlt) && throw(ArgumentError("The minimum altitude must be lower than the maximum altitude."))

    # Loop for the possible repetition times.
    for den = minRep:maxRep
        for num = 0:den-1

            # Check if the fraction num/den is irreducible.
            if gcd(num, den) == 1

                # Loop through the integer parts.
                for ino in intNumOrb
                    addOrbit = false

                    # Number of revolutions per day.
                    numRevPD = ino + num/den

                    (a, i, f1r, f2r, converged) =
                        compute_ss_orbit_by_num_rev_per_day(numRevPD, e)

                    # Check if the orbit is valid.
                    !(@check_orbit(a,e)) && continue

                    # Check if the altitude interval must be verified.
                    if check_altitude
                        if (minAlt < a-R0) && (a-R0 < maxAlt) && (converged)
                            addOrbit = true
                        end
                    else
                        addOrbit = converged
                    end

                    # Check if the orbit must be added to the list.
                    if addOrbit
                        # Compute the period of the orbit considering
                        # perturbations (J2).
                        T = period(a, e, i, :J2)

                        ss_orbits = vcat(ss_orbits, [a a-R0 i T ino num den])
                    end
                end
            end
        end
    end

    # Return the list of orbits.
    ss_orbits
end


"""
    sort_list_ss_orbits_by_height(ss_orbits::Matrix)

Sort the list of Sun-synchronous orbits `ss_orbits` (see
`list_ss_orbits_by_rep_period`) by height and return a new matrix.

"""
sort_list_ss_orbits_by_height(ss_orbits::Matrix) =
    sortslices(ss_orbits, dims=1, by=x->x[1])
