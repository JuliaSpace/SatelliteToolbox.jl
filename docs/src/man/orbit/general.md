General functions
=================

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

This package contains some functions that helps in analysis of orbits.

## Angular velocity

The angular velocity of an object in orbit when considering a Keplerian orbit
(unperturbed model) is given by:

```math
n = n_0 = \sqrt{ \frac{\mu_0}{a^3} }~,
```

where ``\mu_0`` is the standard gravitational parameter for Earth, and ``a`` is
the semi-major axis.

If the perturbation terms up to ``J_2`` are considered, then the angular
velocity is computed by:

```math
n = n_0 + \frac{3}{4} \cdot \frac{R_0^2 \cdot J_2}{a^2\left(1-e^2\right)^2} \cdot n_0 \cdot \left[\sqrt{1-e^2}\cdot(3\cos^2(i)-1) + (5\cos^2(i) - 1) \right]~,
```

where ``e`` is the eccentricity, ``i`` is the inclination, and ``R_0`` is the
Earth equatorial radius.

In this package, the angular velocity \[rad/s] can be computed by the following
functions:

```julia
function angvel(a::Number, e::Number, i::Number, pert::Symbol = :J2)
function angvel(orb::Orbit, pert::Symbol = :J2)
```

where:

* `a` is the semi-major axis \[m];
* `e` is the eccentricity;
* `i` is the inclination \[rad];
* `pert` selects the perturbation terms it should be used, it can be `:J0` or
  `:J2`[^1]; and
* `orb` is an instance of [`Orbit`](@ref).

```jldoctest
julia> angvel(7130982.0, 0.001111, deg2rad(98.405))
0.0010471974485046116

julia> angvel(7130982.0, 0.001111, deg2rad(98.405), :J0)
0.0010484431282179
```

## Time-derivative of the argument of perigee

The time-derivative of the argument of perigee ``\dot{\omega}`` when considering
perturbation terms up to ``J_2`` is:

```math
\dot{\omega} = \frac{3}{4} \cdot \frac{R_0^2 \cdot J_2}{a^2\left(1-e^2\right)^2} \cdot n_0 \cdot (5\cos^2(i) - 1)
```

where ``R_0`` is the Earth equatorial radius, ``a`` is the semi-major axis,
``e`` is the eccentricity, ``i`` is the inclination, and ``n_0`` is the
unperturbed orbital angular velocity.

In the unperturbed model (Keplerian orbit), the time-derivative of the argument
of perigee is always 0.

In this package, the time-derivative of the argument of perigee \[rad/s] can be
computed by the following functions:

```julia
function dArgPer(a::Number, e::Number, i::Number, pert::Symbol = :J2)
function dArgPer(orb::Orbit, pert::Symbol = :J2)
```

where:

* `a` is the semi-major axis \[m];
* `e` is the eccentricity;
* `i` is the inclination \[rad];
* `pert` selects the perturbation terms it should be used, it can be `:J0` or
  `:J2`[^1]; and
* `orb` is an instance of [`Orbit`](@ref).

```jldoctest
julia> dArgPer(7130982, 0.001111, deg2rad(98.405))
-6.082892348533058e-7

julia> dArgPer(7130982, 0.001111, deg2rad(63.435))
-2.433253158726004e-12

julia> dArgPer(7130982, 0.001111, deg2rad(98.405), :J0)
0.0
```

## Time-derivative of the RAAN

The time-derivative of the RAAN (right-ascension of the ascending node)
``\dot{\Omega}`` when considering perturbation terms up to ``J_2`` is:

```math
\dot{\Omega} = -\frac{3}{2} \cdot \frac{R_0^2 \cdot J_2}{a^2\left(1-e^2\right)^2} \cdot n_0 \cdot \cos(i)
```

where ``R_0`` is the Earth equatorial radius, ``a`` is the semi-major axis,
``e`` is the eccentricity, ``i`` is the inclination, and ``n_0`` is the
unperturbed orbital angular velocity.

In the unperturbed model (Keplerian orbit), the time-derivative of the RAAN is
always 0.

In this package, the time-derivative of the RAAN \[rad/s] can be computed by the
following functions:

```julia
function dRAAN(a::Number, e::Number, i::Number, pert::Symbol = :J2)
function dRAAN(orb::Orbit, pert::Symbol = :J2)
```

where:

* `a` is the semi-major axis \[m];
* `e` is the eccentricity;
* `i` is the inclination \[rad];
* `pert` selects the perturbation terms it should be used, it can be `:J0` or
  `:J2`[^1]; and
* `orb` is an instance of [`Orbit`](@ref).

```jldoctest
julia> dRAAN(7130982, 0.001111, deg2rad(98.405))
1.9909533223838115e-7

julia> dRAAN(7130982, 0.001111, deg2rad(98.405), :J0)
0.0
```

## Period

The orbital period of an object in orbit is given by:

```math
T = \frac{2\pi}{n}
```

where ``n`` is the angular velocity as described in [Angular velocity](@ref).

In this package, the orbital period \[s] can be computed by the following
functions:

```julia
function period(a::Number, e::Number, i::Number, pert::Symbol = :J2)
function period(orb::Orbit, pert::Symbol = :J2)
```

where:

* `a` is the semi-major axis \[m];
* `e` is the eccentricity;
* `i` is the inclination \[rad];
* `pert` selects the perturbation terms it should be used, it can be `:J0` or
  `:J2`[^1]; and
* `orb` is an instance of [`Orbit`](@ref).

```jldoctest
julia> period(7130982, 0.001111, deg2rad(98.405))/60
100.00000980636328

julia> period(7130982, 0.001111, deg2rad(98.405), :J0)/60
99.88119746433748
```

---

[^1]: If `pert` is `:J0`, then it will be consider a Keplerian, unperturbed orbit to compute the values. Otherwise, if `pert` is `:J2`, then it will be consider the perturbation terms up to ``J_2`` to compute the values. It `pert` is omitted, then it defaults to `:J2`.

