Orbit propagators
=================

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

Currently, there are three orbit propagators available: **Two Body**, **J2**,
and **SGP4**.  All coded in Julia (no external libraries required).

## Two Body

This algorithm assumes a Keplerian orbit, *i.e.* considers that the Earth is
spherical with the gravitational force computed by Newton's laws.

## J2

This algorithm considers the perturbation terms up to `J2`. The implementation
available here was adapted from [1].

## SGP4

The SGP4 algorithm here was based on [2,3]. It contains the deep space support
that is automatically selected based on the input orbit. Hence, technically, it
is the SPG4/SDP4 algorithm, which will be called just SGP4 here.

## Initialization

All the propagators need to be initialized first using the API function
`init_orbit_proapgator`. The information can be passed in three different ways:

```julia
function init_orbit_proapgator(T, epoch::Number, n_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, M_0::Number)
```

where:

* `T` is the type of the orbit propagator (`Val{:twobody}` for **Two Body**,
  `Val{:J2}` for **J2**, and `Val{:sgp4}` for **SGP4**).
* `epoch`: Initial orbit epoch \[Julian Day].
* `n_0`: Initial angular velocity \[rad/s].
* `e_0`: Initial eccentricity.
* `i_0`: Initial inclination \[rad].
* `Ω_0`: Initial right ascension of the ascending node \[rad].
* `ω_0`: Initial argument of perigee \[rad].
* `M_0`: Initial mean anomaly \[rad].


```julia
function init_orbit_propagator(T, orb_0::Orbit)
```

where:

* `T` is the type of the orbit propagator (`Val{:twobody}` for **Two Body**,
  `Val{:J2}` for **J2**, and `Val{:sgp4}` for **SGP4**).
* `orb_0`: Initial orbital elements (see `Orbit`).

```julia
function init_orbit_propagator(T, tle::TLE)
```

where:

* `T` is the type of the orbit propagator (`Val{:twobody}` for **Two Body**,
  `Val{:J2}` for **J2**, and `Val{:sgp4}` for **SGP4**).
* `tle`: TLE that will be used to initialize the propagator (see [TLE](@ref)).

There are some optional parameters that depend on the orbit propagator type that
can be used to customize the algorithm. Those options are listed as follows:

**Two Body**

* `μ`: Standard gravitational parameter of the central body \[m³/s²]
  (**Default** = `m0`).

**J2 Orbit Propagator**

* `dn_o2`: First time derivative of mean motion divided by 2 \[rad/s²]
  (**Default** = 0).
* `ddn_o6`: Second time derivative of mean motion divided by 6 \[rad/s³]
  (**Default** = 0).
* `j2_gc`: J2 orbit propagator gravitational constants (see `J2_GravCte`)
  (**Default** = `j2_gc_wgs84`).

!!! warning

    The two first options are not available when the TLE is used because this
    information is provided by the TLE.

**SPG4**

* `bstar`: B\* parameter of the SGP4 (**Default** = 0).
* `sgp4_gc`: Gravitational constants (see `SGP4_GravCte`) (**Default** =
  `sgp4_gc_wgs84`).

!!! warning

    The first option is not available when the TLE is used because this
    information is provided by the TLE.

## Propagation

After the orbit propagator is initialized, the orbit can be propagated by the
API functions `propagate!`, `propagate_to_epoch!`, and `step!`.

The function `propagate!` has two signature. The first one is

```julia
function propagate!(orbp, t::Number) where T
```

in which the orbit will be propagated by `t` \[s] **from the orbit epoch**,
which is defined in the initialization and is never changed. This function
returns a tuple with three values:

* The mean Keplerian elements represented in the inertial reference frame
  encapsulated in an instance of the structure `Orbit` \[SI units].
* The position vector represented in the inertial reference frame \[m].
* The velocity vector represented in the inertial reference frame \[m].

The second signature of `propagate!` is:

```julia
function propagate!(orbp, t::AbstractVector) where T
```

where the orbit will be propagated for every value in the vector `t` \[s], which
is a number of seconds **from the orbit epoch**. In this case, an array of
tuples with be returned with each element equivalent to that described for the
first case.

The function `propagate_to_epoch!` also have two signatures similar to
`propagate!`:

```julia
function propagate_to_epoch!(orbp, JD::Number) where T
function propagate_to_epoch!(orbp, JD::AbstractVector) where T
```

It also returns the same information. However, the input argument `JD` is an
epoch \[Julian Day] to which the orbit will be propagated instead of the number
of seconds from the orbit epoch.

!!! warning

    The conversion from Julian Day to seconds that `propagate_to_epoch!` must
    perform can introduce numerical errors.

The `step!` function has the following signature:

```julia
function step!(orbp, Δt::Number)
```

where the orbit is propagated by `Δt` \[s] from the last propagation instant.
This function returns the same information of the first signature of
`propagate!` method.

In all cases, the structure `orbp` is modified by updating the orbit elements
related to the last propagation instant.

!!! note

    All the algorithms can be used to propagate the orbit forward or backward in
    time.

### Reference systems

The inertial reference system in which the propagated values are represented
depends on the reference system used to represent the input data. For TLE
representation, it is very common to use the TEME (True Equator, Mean Equinox)
frame. For more information, see [1].

## Examples

```jldoctest
julia> orbp = init_orbit_propagator(Val{:twobody}, 0.0, 2*pi/6000, 0.001111, 98.405*pi/180, pi/2, 0.0, 0.0);

julia> o,r,v = propagate!(orbp, collect(0:3:24)*60*60);

julia> r
9-element Array{StaticArrays.SArray{Tuple{3},Float64,1,3},1}:
 [5.30792e-7, 7.12871e6, 3.58939e-6]
 [-9.92441e5, 2.19024e6, -6.71674e6]
 [-6.12601e5, -5.78432e6, -4.14603e6]
 [6.12601e5, -5.78432e6, 4.14603e6]
 [9.92441e5, 2.19024e6, 6.71674e6]
 [-5.37523e-7, 7.12871e6, -3.64086e-6]
 [-9.92441e5, 2.19024e6, -6.71674e6]
 [-6.12601e5, -5.78432e6, -4.14603e6]
 [6.12601e5, -5.78432e6, 4.14603e6]
```

```jldoctest
julia> orbp = init_orbit_propagator(Val{:J2}, Orbit(0.0, 7130982.0, 0.001111, 98.405*pi/180, pi/2, 0.0, 0.0));

julia> o,r,v = propagate!(orbp, collect(0:3:24)*60*60);

julia> r
9-element Array{StaticArrays.SArray{Tuple{3},Float64,1,3},1}:
 [5.30372e-7, 7.12306e6, 3.58655e-6]
 [-9.98335e5, 2.14179e6, -6.72549e6]
 [-5.75909e5, -5.83674e6, -4.06734e6]
 [6.65317e5, -5.69201e6, 4.2545e6]
 [9.62557e5, 2.37418e6, 6.65228e6]
 [-1.10605e5, 7.11845e6, -231186.0]
 [-1.02813e6, 1.90664e6, -6.79145e6]
 [-4.82921e5, -5.97389e6, -3.87579e6]
 [750898.0, -5.53993e6, 4.43709e6]
```

```jldoctest
julia> orbp = init_orbit_propagator(Val{:J2}, Orbit(DatetoJD(1986,6,19,0,0,0), 7130982.0, 0.001111, 98.405*pi/180, pi/2, 0.0, 0.0));

julia> o,r,v = propagate_to_epoch!(orbp, DatetoJD(1986,6,19,0,0,0) .+ collect(0:3:24)/24);

julia> r
9-element Array{StaticArrays.SArray{Tuple{3},Float64,1,3},1}:
 [5.30372e-7, 7.12306e6, 3.58655e-6]
 [-9.98335e5, 2.14179e6, -6.72549e6]
 [-5.75909e5, -5.83674e6, -4.06734e6]
 [6.65317e5, -5.69201e6, 4.2545e6]
 [9.62557e5, 2.37418e6, 6.65228e6]
 [-1.10605e5, 7.11845e6, -231186.0]
 [-1.02813e6, 1.90664e6, -6.79145e6]
 [-4.82921e5, -5.97389e6, -3.87579e6]
 [750898.0, -5.53993e6, 4.43709e6]
```

```jldoctest
julia> tle_scd1 = tle"""
       SCD 1
       1 22490U 93009B   18350.91204528  .00000219  00000-0  10201-4 0  9996
       2 22490  24.9683 170.6788 0043029 357.3326 117.9323 14.44539175364603
       """[1];

julia> orbp = init_orbit_propagator(Val{:sgp4}, tle_scd1);

julia> o,r,v = propagate!(orbp, (0:3:24)*60*60);

julia> r
9-element Array{StaticArrays.SArray{Tuple{3},Float64,1,3},1}:
 [2.1104e6, -6.24894e6, 2.71038e6]
 [-5.59246e6, -3.78133e6, 2.1883e6]
 [-5.98838e6, 3.62748e6, -1.13273e6]
 [1.44056e6, 6.29603e6, -3.00473e6]
 [7.02615e6, 791502.0, -1.06173e6]
 [3.607e6, -5.74328e6, 2.21989e6]
 [-4.43043e6, -4.85364e6, 2.68863e6]
 [-6.67554e6, 2.3722e6, -2.79066e5]
 [-1.93293e5, 6.50127e6, -2.89155e6]

julia> v
9-element Array{StaticArrays.SArray{Tuple{3},Float64,1,3},1}:
 [7129.19, 1784.07, -1358.32]
 [4573.31, -5547.04, 2171.25]
 [-3969.35, -5663.64, 2940.09]
 [-7305.14, 1611.56, -49.363]
 [-1211.78, 6739.97, -2945.93]
 [6417.95, 3175.76, -2122.04]
 [5799.59, -4551.63, 1407.45]
 [-2391.64, -6387.69, 3161.66]
 [-7435.44, 128.809, 866.6]

julia> orbp = init_orbit_propagator(Val{:sgp4}, tle_scd1);

julia> o,r,v = step!(orbp, 3*60*60);

julia> o,r,v = step!(orbp, 3*60*60);

julia> o,r,v = step!(orbp, 3*60*60);

julia> o,r,v = step!(orbp, 3*60*60);

julia> o,r,v = step!(orbp, 3*60*60);

julia> o,r,v = step!(orbp, 3*60*60);

julia> o,r,v = step!(orbp, 3*60*60);

julia> o,r,v = step!(orbp, 3*60*60);

julia> r
3-element StaticArrays.SArray{Tuple{3},Float64,1,3}:
 -193293.35025474802
       6.501272877734011e6
      -2.8915511460724953e6

julia> v
3-element StaticArrays.SArray{Tuple{3},Float64,1,3}:
 -7435.439550407856
   128.80933740840044
   866.5999572489231
```

## Low level access

All propagators can be accessed by low-level functions. This allows the user to
have more control about the algorithm and also to reduce the overhead related to
the API functions. If such optimization is necessary, see the functions inside
the directory `./src/orbit/propagators`.

## References

[1] **Vallado, D. A., McClain, W. D (2013).** *Fundamentals of astrodynamics and
applications*. Hawthorne, CA: Microcosm Press.

[2] **Hoots, F. R., Roehrich, R. L (1980).** *Models for Propagation of NORAD
Elements Set*. Spacetrack Report No. 3.

[3] **Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S (2006).**
*Revisiting Spacetrack Report #3: Rev1*. AIAA.

