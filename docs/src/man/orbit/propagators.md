Orbit propagators
=================

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

Currently, there are four orbit propagators available: **Two Body**, **J2**,
**J4** and **SGP4**.  All coded in Julia (no external libraries required).

## Two Body

This algorithm assumes a Keplerian orbit, *i.e.* considers that the Earth is
spherical with the gravitational force computed by Newton's laws.

## J2

This algorithm considers the perturbation terms up to `J2` and the drag effects.
The implementation available here was adapted from [1].

## J4

This algorithm considers the perturbation terms `J2`, `J2²`, and `J4` and the
drag effects. The implementation available here was adapted from [1].

## SGP4

The SGP4 algorithm here was based on [2,3]. It contains the deep space support
that is automatically selected based on the input orbit. Hence, technically, it
is the SPG4/SDP4 algorithm, which will be called just SGP4 here.

## Initialization

All the propagators need to be initialized first using the API function
`init_orbit_proapgator`. The functions signature for each algorithm can be seen
as follows.

### Initialization of Two body, J2, and J4

The orbit propagators **two body**, **J2**, and **J4** can be initialized using
three different methods.

```julia
function init_orbit_proapgator(T, epoch::Number, f_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, f_0::Number)
```

where:

* `T` is the type of the orbit propagator (`Val{:twobody}` for **Two Body**,
  `Val{:J2}` for **J2**, and `Val{:J4}` for **J4**).
* `epoch`: Initial orbit epoch \[Julian Day].
* `a_0`: Initial semi-major axis \[m].
* `e_0`: Initial eccentricity.
* `i_0`: Initial inclination \[rad].
* `Ω_0`: Initial right ascension of the ascending node \[rad].
* `ω_0`: Initial argument of perigee \[rad].
* `f_0`: Initial true anomaly \[rad].

!!! note

    The inputs are the mean orbital elements.

```julia
function init_orbit_propagator(T, orb_0::Orbit)
```

where:

* `T` is the type of the orbit propagator (`Val{:twobody}` for **Two Body**,
  `Val{:J2}` for **J2**, and `Val{:J4}` for **J4**).
* `orb_0`: Initial orbital elements (see `Orbit`).

```julia
function init_orbit_propagator(T, tle::TLE)
```

where:

* `T` is the type of the orbit propagator (`Val{:twobody}` for **Two Body**,
  `Val{:J2}` for **J2**, and `Val{:J4}` for **J4**).
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

**J4 Orbit Propagator**

* `dn_o2`: First time derivative of mean motion divided by 2 \[rad/s²]
  (**Default** = 0).
* `ddn_o6`: Second time derivative of mean motion divided by 6 \[rad/s³]
  (**Default** = 0).
* `j4_gc`: J4 orbit propagator gravitational constants (see `J4_GravCte`)
  (**Default** = `j4_gc_wgs84`).

!!! warning

    The two first options are not available when the TLE is used because this
    information is provided by the TLE.

### Initialization using the angular velocity

If the orbit is defined in terms of the angular velocity (mean motion) instead
of the semi-major axis, then it is possible to use the function `angvel_to_a` to
convert:

```julia
function angvel_to_a(n::Number, e::Number, i::Number, pert::Symbol = :J2)
```

It computes the semi-major axis that will provide an angular velocity `n`
[rad/s] in an orbit with eccentricity `e` and inclination `i` [rad], using the
perturbation terms specified by the symbol `pert`.

Notice that the angular velocity `n` is related to the nodal period, *i.e.* the
time between two consecutive passages by the ascending node.

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.
* `:J4`: Consider the perturbation terms J2, J4, and J2².

If `pert` is omitted, then it defaults to `:J2`.

### Initialization of SGP4

The SGP4/SDP4 propagator is meant to be used together with a TLE. Hence, the
initialization using user-defined orbital elements is not available through the
API. If this is really required, then the user must access the low-level
function `sgp4_init`.

The API function to initialize the SGP4 using a TLE is:

```julia
function init_orbit_propagator(T, tle::TLE, sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84)
```

where:

* `T` must be `Val{:sgp4}`;
* `tle`: TLE that will be used to initialize the propagator (see [TLE](@ref)).
* `sgp4_gc`: Gravitational constants (see `SGP4_GravCte`) (**Default** =
  `sgp4_gc_wgs84`).

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
julia> orbp = init_orbit_propagator(Val{:twobody}, 0.0, 7130982.0, 0.001111, 98.405*pi/180, pi/2, 0.0, 0.0);

julia> o,r,v = propagate!(orbp, collect(0:3:24)*60*60);

julia> r
9-element Array{StaticArrays.SArray{Tuple{3},Float64,1,3},1}:
 [5.30372e-7, 7.12306e6, 3.58655e-6]
 [-987245.0, 2.2796e6, -6.68158e6]
 [-6.3457e5, -5.6651e6, -4.29471e6]
 [5.77611e5, -5.94385e6, 3.90922e6]
 [1.00749e6, 1.82033e6, 6.8186e6]
 [70133.2, 7.1069e6, 4.74655e5]
 [-9.62529e5, 2.72855e6, -6.5143e6]
 [-6.88667e5, -5.3608e6, -4.66083e6]
 [5.18048e5, -6.1958e6, 3.50609e6]
```

```jldoctest
julia> orbp = init_orbit_propagator(Val{:J2}, Orbit(0.0, 7130982.0, 0.001111, 98.405*pi/180, pi/2, 0.0, 0.0));

julia> o,r,v = propagate!(orbp, collect(0:3:24)*60*60);

julia> r
9-element Array{StaticArrays.SArray{Tuple{3},Float64,1,3},1}:
 [5.30372e-7, 7.12306e6, 3.58655e-6]
 [-9.9635e5, 2.18638e6, -6.71137e6]
 [-587229.0, -5.78237e6, -4.14257e6]
 [6.49425e5, -5.77559e6, 4.143e6]
 [9.72829e5, 2.1969e6, 6.71164e6]
 [-76505.6, 7.12265e6, 503.253]
 [-1.01976e6, 2.17562e6, -6.71109e6]
 [-5.24963e5, -5.78847e6, -4.14214e6]
 [711543.0, -5.76813e6, 4.14344e6]
```

```jldoctest
julia> orbp = init_orbit_propagator(Val{:J2}, Orbit(DatetoJD(1986,6,19,0,0,0), 7130982.0, 0.001111, 98.405*pi/180, pi/2, 0.0, 0.0));

julia> o,r,v = propagate_to_epoch!(orbp, DatetoJD(1986,6,19,0,0,0) .+ collect(0:3:24)/24);

julia> r
9-element Array{StaticArrays.SArray{Tuple{3},Float64,1,3},1}:
 [5.30372e-7, 7.12306e6, 3.58655e-6]
 [-9.9635e5, 2.18638e6, -6.71137e6]
 [-587229.0, -5.78237e6, -4.14257e6]
 [6.49425e5, -5.77559e6, 4.143e6]
 [9.72829e5, 2.1969e6, 6.71164e6]
 [-76505.6, 7.12265e6, 503.253]
 [-1.01976e6, 2.17562e6, -6.71109e6]
 [-5.24963e5, -5.78847e6, -4.14214e6]
 [711543.0, -5.76813e6, 4.14344e6]
```

```jldoctest
julia> orbp = init_orbit_propagator(Val{:J4}, Orbit(DatetoJD(1986,6,19,0,0,0), 7130982.0, 0.001111, 98.405*pi/180, pi/2, 0.0, 0.0));

julia> o,r,v = propagate!(orbp, (0:3:24)*60*60);

julia> r
9-element Array{StaticArrays.SArray{Tuple{3},Float64,1,3},1}:
 [5.30372e-7, 7.12306e6, 3.58655e-6]
 [-996359.0, 2.18621e6, -6.71142e6]
 [-587181.0, -5.78257e6, -4.14229e6]
 [6.49494e5, -5.77529e6, 4.14342e6]
 [972787.0, 2.19756e6, 6.71143e6]
 [-76651.4, 7.12265e6, -351.741]
 [-1.0198e6, 2.17463e6, -6.7114e6]
 [-5.24788e5, -5.78918e6, -4.14117e6]
 [7.11718e5, -5.76732e6, 4.14455e6]
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
 -193293.3502548483
       6.501272877734009e6
      -2.8915511460724827e6

julia> v
3-element StaticArrays.SArray{Tuple{3},Float64,1,3}:
 -7435.439550407853
   128.80933740830324
   866.5999572489661
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

