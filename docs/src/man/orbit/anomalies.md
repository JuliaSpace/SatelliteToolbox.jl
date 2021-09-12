Anomalies
=========

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

There are three types of anomalies[^1] that can be used to describe the position
of the satellite in the orbit plane with respect to the argument of perigee:

* The mean anomaly (`M`);
* The eccentric anomaly (`E`); and
* The true anomaly (`f`).

This package contains the following functions that can be used to convert one to
another:

```julia
function M_to_E(e::Number, M::Number; kwargs...)
function M_to_f(e::Number, M::Number; kwargs...)
function E_to_f(e::Number, E::Number)
function E_to_M(e::Number, E::Number)
function f_to_E(e::Number,f::Number)
function f_to_M(e::Number, f::Number)
```

where:

* `M` is the mean anomaly [rad];
* `E` is the eccentric anomaly [rad];
* `f` is the true anomaly [rad];
* `e` is the eccentricity.

All the returned values are in [rad].

The functions `M_to_E` and `M_to_f` uses the Newton-Raphson algorithm to solve
the Kepler's equation. In this case, the following keywords are available to
configure it:

- `tol::Union{Nothing, Number}`: Tolerance to accept the solution from
    Newton-Raphson algorithm. If `tol` is `nothing`, then it will be
    `eps(T)`, where `T` is a floating-point type obtained from the promotion of
    `T1` and `T2`. (**Default** = `nothing`)
- `max_iterations::Number`: Maximum number of iterations allowed for the
    Newton-Raphson algorithm. If it is lower than 1, then it is set to 10.
    (**Default** = 10)

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
