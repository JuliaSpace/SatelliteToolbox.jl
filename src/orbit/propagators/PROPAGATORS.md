SatelliteToolbox.jl orbit propagators
=====================================

This document describes the design that a propagator must have inside
SatelliteToolbox.jl.

## Initialization

Each propagator must have an initialization function that receives the input
elements and returns a propagator structure. This entity contains all the
initialized variables related to the propagator.

The type of the input elements varies according to the propagator.

No restriction is imposed on the output structure, which can also contain
variables to reduce the computational burden. However, it **must not** have any
variable of the SatelliteToolbox.jl API such as `Orbit`.

This function shall be named as `<propagator identifier>_init`.

```julia
j2d = j2_init(j2_gc, epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, dn_o2, ddn_o6)
```

## Propagation

Each propagator must have a propagation function that receives the propagator
structure and a time `t`. It must return the position and velocity after `t` [T]
from the initialized epoch. The unit of `t` varies according to the propagator.

The propagation must be performed ideally in the same reference frame in which
the input elements were represented. If this is not the case, the documentation
must clearly state it.

The output vectors must be in the same reference frame in which the input
elements were represented. Hence, if the propagator transforms the inputs to
another frame, then it must convert them back.

For the sake of simplification, it is advised, although not imposed, that the
propagator structure is updated with the latest orbital elements using any
representation.

The propagator should return the vectors preferably in SI. If this is not the
case, the documentation must state clearly.

This function shall be named as `<propagator identified>!`.

```julia
r_i, v_i = j2!(j2d, 10)
```

## API

SatelliteToolbox.jl has a propagator API to improve the usability. This API has
the following requirements, which must be met by all propagators that uses it.

Every propagator must have a structure derived from `OrbitPropagator` with the
following requirement:

```julia
struct OrbitPropagator<Propagator name>{Tepoch, T} <: OrbitPropagator{Tepoch, T}
  <Any field required by the propagator>
end
```

where `Tepoch` is the type used to represent the epoch of the input elements,
whereas `T` is the type used for the internal variables.

### Initialization

The initialization is performed by the function:

```julia
init_orbit_propagator(T, args...; kwargs...)
```

where `T = Val{<Orbit propagator symbol>}`, and it must return and object of
type `OrbitPropagator<Propagator name>`. The arguments and keywords depends on
the propagator and must be documented in the docstring. The propagator must
record the epoch during the initialization, which must be kept constant
during the entire object existence. It also needs to record the instant of the
last propagation.

### Epoch

Each propagator must return the initial element epoch in Julian Day by the
function:

```julia
get_epoch(orbp)
```

Notice that this value must never change during the existence of the object
`orbp`.

### Propagation

The following functions must be overloaded by each propagator.

```julia
propagate!(orbp, t)
```

Propagate the orbit of propagator `orbp` by `t` [s] from the epoch. This
function must return the propagated position and velocity represented in the
same reference frame used in the initialization. All units must be in SI.

```julia
step!(orbp, dt)
```

Propagate the orbit of propagator `orbp` by `dt` [s] from the instant of the
last propagation. This function must return the propagated position and velocity
represented in the same reference frame used in the initialization. All units
must be in SI.

We also have the function `propagate_to_epoch!`, but the default implementation
should work for all propagators.

### Mean elements (optional)

The function

```julia
get_mean_elements(orbp)
```

should return the mean elements using the structure `KeplerianElements` related
to the latest propagation. Notice that this is an **optional** feature. If the
propagator does not implement it, then it will return `nothing`.
