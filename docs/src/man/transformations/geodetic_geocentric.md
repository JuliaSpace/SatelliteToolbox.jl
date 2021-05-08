Geodetic and Geocentric
=======================

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

There are three functions that can help to convert between geodetic and
geocentric representations. Notice that currently all Geodetic representations
are based on the WGS84 reference ellipsoid.

## ECEF to Geodetic

It is possible to convert a position vector represented in an Earth-Centered,
Earth-Fixed frame (ECEF) `r_e` to the Geodetic latitude, longitude, and altitude
by the following function:

```julia
function ecef_to_geodetic(r_e::AbstractVector)
```

which returns a tuple with:

* The Geocentric latitude [rad];
* The longitude [rad]; and
* The altitude above the reference ellipsoid [m].

```jldoctest
julia> ecef_to_geodetic([R0;0;0])
(0.0, 0.0, 0.0)

julia> ecef_to_geodetic([0;R0;0])
(0.0, 1.5707963267948966, 0.0)

julia> ecef_to_geodetic([0;0;R0])
(1.5707963267948966, 0.0, 21384.685754820704)
```

# Geodetic to ECEF

The Geodetic latitude `lat` \[rad], longitude `lon` \[rad], and altitude `h`
\[m] can be converted to a vector represented in an ECEF reference frame by the
following function:

```julia
function geodetic_to_ecef(lat::Number, lon::Number, h::Number)
```

in which a 3x1 vector will be returned.

```jldoctest
julia> geodetic_to_ecef(0,0,0)
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
 6.378137e6
 0.0
 0.0

julia> geodetic_to_ecef(deg2rad(-22),deg2rad(-45),0)
3-element StaticArrays.SVector{3, Float64} with indices SOneTo(3):
  4.1835869067109847e6
 -4.1835869067109837e6
 -2.3744128953028163e6
```

# Geodetic to Geocentric

Given a Geodetic latitude `ϕ_gd` \[rad] and altitude above the reference
ellipsoid `h` \[m], one can obtain the Geocentric coordinates (Geocentric
latitude and position from the center of Earth) using the following function:

```julia
function GeodetictoGeocentric(ϕ_gd::Number, h::Number)
```

in which a tuple with two values will be returned:

* The Geocentric latitude \[rad]; and
* The distance from the center of Earth \[m].

!!! note

    The longitude is the same in both Geodetic and Geocentric representations.

```jldoctest
julia> GeodetictoGeocentric(deg2rad(-22), 0)
(-0.38164509973650357, 6.375157677217675e6)

julia> GeodetictoGeocentric(0,0)
(0.0, 6.378137e6)
```
