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

This function shall be named as `<propagator identified>!`.

```julia
r_i, v_i = j2!(j2d, 10)
```
