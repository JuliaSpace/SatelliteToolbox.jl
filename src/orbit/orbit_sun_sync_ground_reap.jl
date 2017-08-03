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
#    Many auxiliary functions for ground repeating, Sun-synchronous (GRSS) orbit
#    computations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2015-07-15: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export adjacent_track_distance_grss, adjacent_track_angle_grss
export compute_ss_orbit_by_num_rev_per_day
export list_ss_orbits_by_rep_period, sort_list_ss_orbits_by_height

"""
### function adjacent_track_distance_grss(T::Real, i::Real, To::Real, lat::Real)

Compute the distance between adjacent ground tracks in at a given latitude for a
ground repeating, Sun-synchronous orbit.

##### Args

* T: Orbit period [s].
* i: Inclination [rad].
* To: Orbit cycle [days].
* lat: Latitude [rad].

##### Returns

* The distance between adjacent ground tracks at the given latitude [m].

##### Remarks

The functions *does not* check if the orbit is a GRSS orbit.

"""

function adjacent_track_distance_grss(T::Real, i::Real, To::Real, lat::Real)
    # Angle between two adjacent traces.
    theta = (T/To)*pi/43200.0

    # Angle of swath.
    beta = asin(sin(theta)*sin(i))

    # Minimum swath.
    beta*R0*cos(lat)
end


"""
### function adjacent_track_distance_grss(a::Real, e::Real, i::Real, To::Real, lat::Real)

Compute the distance between adjacent ground tracks in at a given latitude for a
ground repeating, Sun-synchronous orbit.

##### Args

* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].
* To: Orbit cycle [days].
* lat: Latitude [rad].

##### Returns

* The distance between adjacent ground tracks at the given latitude [m].

##### Remarks

The functions *does not* check if the orbit is a GRSS orbit.

"""

function adjacent_track_distance_grss(a::Real, e::Real, i::Real, To::Real,
                                      lat::Real)
    # Orbit period.
    T = t_J2(a,e,i)

    # Compute the distance between adjacent ground tracks.
    adjacent_track_distance_grss(T, i, To, lat)
end

"""
### function adjacent_track_angle_grss(h::Real, T::Real, i::Real, To::Real, lat::Real)

Compute the angle between two adjacent ground tracks in a given latitude
measured on the satellite position for a ground repeating, Sun-synchronous
orbit.

##### Args

* h: Orbit altitude in the Equator [m].
* T: Orbit period [s].
* i: Inclination [rad].
* To: Orbit cycle [days].
* lat: Latitude.

##### Returns

* The angle between adjacent ground tracks at the given latitude [rad].

##### Remarks

The functions *does not* check if the orbit is a GRSS orbit.

"""

function adjacent_track_angle_grss(h::Real, T::Real, i::Real, To::Real,
                                   lat::Real)
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
### function adjacent_track_angle_grss(h::Real, a::Real, e::Real, i::Real, To::Integer, lat::Real)

Compute the angle between two adjacent ground tracks in a given latitude
measured on the satellite position for a ground repeating, Sun-synchronous
orbit.

##### Args

* h: Orbit altitude in the Equator [m].
* a: Semi-major axis [m].
* e: Eccentricity.
* i: Inclination [rad].
* To: Orbit cycle [days].
* lat: Latitude.

##### Returns

* The angle between adjacent ground tracks at the given latitude [rad].

##### Remarks

The functions *does not* check if the orbit is a GRSS orbit.

"""

function adjacent_track_angle_grss(h::Real, a::Real, e::Real, i::Real,
                                   To::Integer, lat::Real)
    # Period (J2).
    T = t_J2(a,e,i)

    # Compute the minimum half FOV.
    adjacent_track_angle_grss(h, T, i, To, lat)
end

"""
### function compute_ss_orbit_by_num_rev_per_day(numRevPD::Real, e::Real)

Compute the Sun-synchronous orbit given the number of revolutions per day and
the eccentricity.

##### Args

* numRevPD: Number of revolutions per day.
* e: Eccentricity.

##### Returns

* The semi-major axis [m].
* The inclination [rad].
* The residues of the two functions.
* A boolean variable that indicates if the numerical algorithm converged.

"""

function compute_ss_orbit_by_num_rev_per_day(numRevPD::Real, e::Real)
    compute_ss_orbit_by_ang_vel(numRevPD*2*pi/86400.0, e)
end

"""
### function list_ss_orbits_by_rep_period(minRep::Int, maxRep::Int, minAlt::Real=-1.0, maxAlt::Real=-1.0, e::Real=0.0)

Compute a list of repeating Sun-synchronous orbits.

##### Args

* minRep: Minimum repetition time of the orbit [days].
* maxRep: Maximum repetition time of the orbit [days].
* minAlt: Minimum altitude of the orbits on the list [m].
* maxAlt: Minimum altitude of the orbits on the list [m].
* e: Eccentricity.

##### Returns

* A matrix containing the orbits found.

##### Remarks


1) If minAlt or maxAlt is < 0.0, then the altitude will not be checked when a
orbit is added to the list.
2) The output matrix has the following format:

        Semi-major axis [m] | Altitude [m] | Period [s] | Int | Num | Den
        --------------------|--------------|------------|-----|-----|----

in which the period is Int + Num/Den.

"""

function list_ss_orbits_by_rep_period(minRep::Int,       maxRep::Int,
                                      minAlt::Real=-1.0, maxAlt::Real=-1.0,
                                      e::Real=0.0)
    numRevPD = 0
    # Check if the arguments are valid.
    if (minRep <= 0)
        throw(ArgumentError("The minimum repetition time must be greater than 0."))
    end

    if (maxRep <= 0)
        throw(ArgumentError("The maximum repetition time must be greater than 0."))
    end

    if (minRep > maxRep)
        throw(ArgumentError("The minimum repetition time must be smaller or equal to the maximum repetition time."))
    end

    if !( 0. <= e < 1. )
        throw(ArgumentError("The eccentricity must be within the interval 0 <= e < 1."))
    end

    # Matrix to store the available orbits.
    ss_orbits = Array{Float64}(0, 7)

    # Integer part of the number of orbits per day.
    intNumOrb = [13., 14., 15., 16., 17.]

    # Loop for the possible repetition times.
    for den = minRep:maxRep
        for num = 0:den-1
            # Check if the fraction num/den is irreducible.
            if ( gcd(num, den) == 1.0 )
                # Loop through the integer parts.
                for ino in intNumOrb
                    addOrbit = false

                    # Number of revolutions per day.
                    numRevPD = ino+Float64(num)/Float64(den)

                    (a, i, f1r, f2r, converged) =
                        compute_ss_orbit_by_num_rev_per_day(numRevPD, e)

                    # Check if the orbit is valid.
                    (!is_orbit_valid(a, e)) && continue

                    # Check if the altitude interval must be verified.
                    if (minAlt > 0) && (maxAlt > 0)
                        if (minAlt < a-R0) && (a-R0 < maxAlt) && (converged)
                            addOrbit = true
                        end
                    else
                        addOrbit = converged
                    end

                    # Check if the orbit must be added to the list.
                    if (addOrbit)
                        # Compute the period of the orbit considering
                        # perturbations (J2).
                        period = t_J2(a, e, i)

                        ss_orbits =
                            vcat(ss_orbits, [a a-R0 i period ino num den])
                    end
                end
            end
        end
    end

    # Return the list of orbits.
    ss_orbits
end


"""
### sort_list_ss_orbits_by_height(ss_orbits::Array{Float64, 2})

Sort the list of Sun-synchronous orbits by height.

##### Args

* ss_orbits: List of Sun-synchronous orbits (*see*
list_ss_orbits_by_rep_period).

##### Returns

* A matrix containing a list of Sun-synchronous orbits sorted by height.

"""

sort_list_ss_orbits_by_height(ss_orbits::Array{Float64, 2}) =
    sortrows(ss_orbits, by=x->x[1])
