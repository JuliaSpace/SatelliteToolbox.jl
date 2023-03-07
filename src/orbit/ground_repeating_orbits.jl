# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Auxiliary functions for ground repeating orbit computations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Remarks
# ==============================================================================
#
#   A ground repeating orbit is any orbit that the number of revolutions per day
#   is a rational number. Hence, this type of orbit repeats its ground trace
#   after a finite number of days.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export ground_repeating_orbit_adjacent_track_angle
export ground_repeating_orbit_adjacent_track_distance

"""
    ground_repeating_orbit_adjacent_track_angle(h::T1, orbit_period::T2, i::T3, orbit_cycle::Integer; kwargs...) where {T1 <: Number, T2 <: Number, T3 <: Number}

Compute the adjacent track distance [m] at Equator in a ground repeating orbit.
The orbit is described by its altitude at Equator `h` [m], orbital period
`orbit_period` [s], inclination `i` [rad], and orbit cycle `orbit_cyle` [day].

!!! warning
    The code does not check if the altitude `h` is correct for the orbit with
    the period `orbit_period`.

!!! note
    Internaly, this function uses the precision obtained by promoting `T1`,
    `T2`, and `T3` to a float-pointing number.

# Keywords

- `R0::Number`: Earth's radius at Equator. (**Default** = `R0`)

# Extended help

## Remarks

A ground repeating orbit is any orbit that the number of revolutions per day is
a rational number. Hence, this type of orbit repeats its ground trace after a
finite number of days.

The information `orbit_period` and `orbit_cyle` is redundant. However, they are
necessary to improve the algorithm precision. Otherwise, the `orbit_cycle` must
be obtained by converting the floating-point number `orbit_period` to a rational
number, leading to numerical problems.
"""
function ground_repeating_orbit_adjacent_track_angle(
    h::T1,
    orbit_period::T2,
    i::T3,
    orbit_cycle::Integer;
    R0::Number = R0
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    T = float(promote_type(T1, T2, T3))

    # Angle between two adjacent tracks at Equator measured from the Earth's
    # center.
    θ = T(orbit_period) / T(orbit_cycle) * T(π) / 43200

    # Angle between two adjacent tracks along the satellite's path measured from
    # the Earth's center. The projection of this angle in the Earth's surface is
    # the distance between the two adjacent tracks.
    β = asin(sin(θ) * sin(T(i)))

    # Compute the angle between the two ground tracks measured from the
    # satellite. α is an auxiliary angle and γ is the angle we are looking for.
    sin_βo2, cos_βo2 = sincos(β / 2)
    α = √(T(R0)^2 + (T(R0) + T(h))^2 - 2T(R0) * (T(R0) + T(h)) * cos_βo2)
    γ = asin(T(R0) / α * sin_βo2)

    return γ
end

"""
    ground_repeating_orbit_adjacent_track_distance(orbit_period::T1, i::T2, orbit_cycle::Integer; kwargs...) where {T1 <: Number, T2 <: Number}

Compute the adjacent track distance [m] at Equator in a ground repeating orbit.
The orbit is described by its orbital period `orbit_period` [s], inclination `i`
[rad], and orbit cycle `orbit_cyle` [day].

!!! note
    Internaly, this function uses the precision obtained by promoting `T1` and
    `T2` to a float-pointing number.

# Keywords

- `R0::Number`: Earth's radius at Equator. (**Default** = `R0`)

# Extended help

## Remarks

A ground repeating orbit is any orbit that the number of revolutions per day is
a rational number. Hence, this type of orbit repeats its ground trace after a
finite number of days.

The information `orbit_period` and `orbit_cyle` is redundant. However, they are
necessary to improve the algorithm precision. Otherwise, the `orbit_cycle` must
be obtained by converting the floating-point number `orbit_period` to a rational
number, leading to numerical problems.
"""
function ground_repeating_orbit_adjacent_track_distance(
    orbit_period::T1,
    i::T2,
    orbit_cycle::Integer;
    R0::Number = R0
) where {T1 <: Number, T2 <: Number}
    T = float(promote_type(T1, T2))

    # Angle between two adjacent tracks at Equator measured from the Earth's
    # center.
    θ = T(orbit_period) / T(orbit_cycle) * T(π) / 43200

    # Angle between two adjacent tracks along the satellite's path measured from
    # the Earth's center. The projection of this angle in the Earth's surface is
    # the distance between the two adjacent tracks.
    β = asin(sin(θ) * sin(T(i)))

    # Distance between two adjacent tracks on the Earth's surface.
    d = T(R0) * β

    return d
end
