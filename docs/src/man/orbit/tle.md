TLE
===

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

The TLE, or Two Line Elements set, is a data format that contains information
about the orbit at a specific epoch of an Earth-orbiting object. The information
is split into two lines with 70 characters each. In the following, it is
presented an example of a TLE describing the orbit of the Brazilian satellite
SCD-1 at 25 December 2018:

```
SCD 1                   
1 22490U 93009B   18359.76217587  .00000186  00000-0  84512-6 0  9998
2 22490  24.9694 116.1709 0043211  90.3968  62.0083 14.44539396366163
```

For more information, see
[https://en.wikipedia.org/wiki/Two-line\_element\_set](https://en.wikipedia.org/wiki/Two-line_element_set).

The TLE contains all the necessary information to propagate the orbit using, for
example, the [SGP4](@ref) orbit propagator. Hence, this package contains a set
of functions that helps to load the TLE information to be used in the [Orbit
propagators](@ref).

If the TLEs are stored in one file, then they can be loaded using the function:

```julia
function read_tle(tle_filename::String, verify_checksum::Bool = true)
```

where the `tle_filename` is the file path. Each TLE line has a checksum to
verify the correctness of the data. By default, if the checksum is wrong, then
this function will throw an error. The checksum verification can be avoided by
setting `verify_checksum` to `false`.

This function returns an array of TLEs. Each TLE is an instance of the structure
[`TLE`](@ref).

```julia-repl
julia> tles = read_tle("tles")
2-element Vector{TLE}:
 TLE: SCD 1 (Epoch = 2018-12-25T18:17:31.995)
 TLE: SCD 2 (Epoch = 2018-12-26T05:50:03.289)

julia> tles[1]
                             TLE
    ==========================================================
                            Name: SCD 1
                Satellite number: 22490
        International designator: 93009B
                    Epoch (Year): 18
                     Epoch (Day): 359.76217587
              Epoch (Julian Day): 2458478.26218
              Element set number: 999
                    Eccentricity:   0.00432110 deg
                     Inclination:  24.96940000 deg
                            RAAN: 116.17090000 deg
             Argument of perigee:  90.39680000 deg
                    Mean anomaly:  62.00830000 deg
                 Mean motion (n):  14.44539396 revs/day
               Revolution number: 36616

                              B*: 0.000001 1/[er]

                        1   d
                       ---.--- n: 0.000002 rev/day²
                        2  dt

                        1   d²
                       ---.--- n: 0.000000 rev/day³
                        6  dt²
    ==========================================================

julia> tles[2]
                             TLE
    ==========================================================
                            Name: SCD 2
                Satellite number: 25504
        International designator: 98060A
                    Epoch (Year): 18
                     Epoch (Day): 360.24309362
              Epoch (Julian Day): 2458478.74309
              Element set number: 999
                    Eccentricity:   0.00174310 deg
                     Inclination:  24.99670000 deg
                            RAAN: 319.86640000 deg
             Argument of perigee: 121.38100000 deg
                    Mean anomaly:   9.79390000 deg
                 Mean motion (n):  14.44057429 revs/day
               Revolution number: 6554

                              B*: 0.000011 1/[er]

                        1   d
                       ---.--- n: 0.000002 rev/day²
                        2  dt

                        1   d²
                       ---.--- n: 0.000000 rev/day³
                        6  dt²
    ==========================================================
```

If the TLE is stored in a string, then it can be read using the following
functions:

```julia
function read_tle_from_string(tles::String, verify_checksum::Bool = true)
function read_tle_from_string(tle_l1::String, tle_l2::String, verify_checksum::Bool = false)
```

In the first case, a list of TLEs can be passed in `tles`. In the second case,
the first line and second line of a TLE can be passed in `tle_l1` and `tle_l2`,
respectively. In both cases an array of TLEs is returned.  The argument
`verify_checksum` has the same function as described previously.

```jldoctest
julia> tles = """
       SCD 1
       1 22490U 93009B   18359.76217587  .00000186  00000-0  84512-6 0  9998
       2 22490  24.9694 116.1709 0043211  90.3968  62.0083 14.44539396366163
       SCD 2
       1 25504U 98060A   18360.24309362  .00000218  00000-0  10518-4 0  9996
       2 25504  24.9967 319.8664 0017431 121.3810   9.7939 14.44057429 65541
       """;

julia> read_tle_from_string(tles)
2-element Vector{TLE}:
 TLE: SCD 1 (Epoch = 2018-12-25T18:17:31.995)
 TLE: SCD 2 (Epoch = 2018-12-26T05:50:03.289)

julia> tle_l1 = "1 22490U 93009B   18359.76217587  .00000186  00000-0  84512-6 0  9998"
"1 22490U 93009B   18359.76217587  .00000186  00000-0  84512-6 0  9998"

julia> tle_l2 = "2 22490  24.9694 116.1709 0043211  90.3968  62.0083 14.44539396366163"
"2 22490  24.9694 116.1709 0043211  90.3968  62.0083 14.44539396366163"

julia> read_tle_from_string(tle_l1, tle_l2)
1-element Vector{TLE}:
 TLE: UNDEFINED (Epoch = 2018-12-25T18:17:31.995)
```

It is also available two special types of strings, `tle"` and `tlenc"`, that
automatically loads a set of TLEs. In the first case the checksum is verified
whereas in the second case it is not.

```julia-repl
julia> tles = tle"""
       SCD 1
       1 22490U 93009B   18359.76217587  .00000186  00000-0  84512-6 0  9998
       2 22490  24.9694 116.1709 0043211  90.3968  62.0083 14.44539396366163
       SCD 2
       1 25504U 98060A   18360.24309362  .00000218  00000-0  10518-4 0  9996
       2 25504  24.9967 319.8664 0017431 121.3810   9.7939 14.44057429 65541
       """
2-element Vector{TLE}:
 TLE: SCD 1 (Epoch = 2018-12-25T18:17:31.995)
 TLE: SCD 2 (Epoch = 2018-12-26T05:50:03.289)

julia> tles = tle"""
       SCD 1
       1 22490U 93009B   18359.76217587  .00000186  00000-0  84512-6 0  9998
       2 22490  24.9694 116.1709 0043211  90.3968  62.0083 14.44539396366164
       """
ERROR: LoadError: The TLE file is not valid (error in line 3): Expected checksum: 3, line checksum: 4.
...

julia> tles = tlenc"""
       SCD 1
       1 22490U 93009B   18359.76217587  .00000186  00000-0  84512-6 0  9998
       2 22490  24.9694 116.1709 0043211  90.3968  62.0083 14.44539396366164
       """
1-element Array{TLE,1}:
 TLE: SCD 1 (Epoch = 2018-12-25T18:17:31.995)
```
