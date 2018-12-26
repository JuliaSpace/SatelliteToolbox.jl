Anomalies
=========

There are three types of anomalies[^1] that can be used to describe the position
of the satellite in the orbit plane with respect to the argument of perigee:

* The mean anomaly (`M`);
* The eccentric anomaly (`E`); and
* The true anomaly (`f`).

This package contains the following functions that can be used to convert one to
another:

```julia
function M_to_E(e::Number, M::Number, tol::Number = 1e-10)
function M_to_f(e::Number, M::Number, tol::Number = 1e-10)
function E_to_f(e::Number, E::Number)
function E_to_M(e::Number, E::Number)
function f_to_E(e::Number,f::Number)
function f_to_E(orb::Orbit)
function f_to_M(e::Number, f::Number)
function f_to_M(orb::Orbit)
```

where:

* `M` is the mean anomaly \[rad];
* `E` is the eccentric anomaly \[rad];
* `f` is the true anomaly \[rad];
* `e` is the eccentricity;
* `orb` is an instance of the structure [`Orbit`](@ref);
* `tol` is used to select the tolerance for the cases in which the conversion is
  performed by a numerical method, such as the Newton-Raphson algorithm.

All the returned values are in \[rad].

```jldoctest
julia> M_to_E(0.04, pi/4)
0.8144932819286269

julia> M_to_f(0.04, pi/4)
0.8440031124631683

julia> f_to_M(0.04, pi/4)
0.7300148523821107

julia> M_to_f(0, 0.343)
0.3430000000000001

julia> M_to_f(0.04, 0.343)
0.37122803399203647
```

[^1]: In astronomy, anomaly is an angle.
