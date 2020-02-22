# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Types and structures of TLE.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export TLE

"""
    TLE

This structure contains the same elements of the TLE with the same units.

# Fields

* `name`: Name of the satellite.

## First line

* `sat_num`: Satellite number.
* `classification`: Classification ('U', 'C', or 'S').
* `int_designator`: International designator.
* `epoch_year`: Epoch year (two digits).
* `epoch_day`: Epoch day (day + fraction of the day).
* `epoch`: The epoch represented in Julian Day.
* `dn_o2`: 1st time derivative of mean motion / 2 [rev/day²].
* `ddn_o6`: 2nd time derivative of mean motion / 6 [rev/day³].
* `bstar`: B* drag term.
* `elem_set_number`: Element set number.
* `checksum_l1`: Checksum of the line 1 (modulo 10).

## Second line

* `i`: Inclination [deg].
* `Ω`: Right ascension of the ascending node [deg].
* `e`: Eccentricity.
* `ω`: Argument of perigee [deg].
* `M`: Mean anomaly [deg].
* `n`: Mean motion [rev/day].
* `rev_num`: Revolution number at epoch [rev].
* `checksum_l2`: Checksum of the line 2 (modulo 10).

"""
@with_kw_noshow struct TLE
    name::String

    # First line
    # ==========
    sat_num::Int
    classification::Char
    int_designator::String
    epoch_year::Int
    epoch_day::Float64
    epoch::Float64
    dn_o2::Float64
    ddn_o6::Float64
    bstar::Float64
    elem_set_number::Int
    checksum_l1::Int

    # Second line
    # ===========

    i::Float64
    Ω::Float64
    e::Float64
    ω::Float64
    M::Float64
    n::Float64
    rev_num::Int
    checksum_l2
end
