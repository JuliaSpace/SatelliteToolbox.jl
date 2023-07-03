ISS Observation Instants
========================

In this tutorial, we will compute the instants in which we can theoretically observe the ISS
from our location. Actually, we will compute the moments an observer on the Earth's surface
has a direct line of sight to the ISS.

## Theory

We can verify if we have a line of sight to an object by computing the elevation ``\lambda``
of its position vector represented in a local reference frame. Let's use the NED
(North-East-Down) reference: its X axis points toward the North, its Y axis points toward
the East, and its Z axis points toward the Earth's center, as shown in the following figure.

![Elevation in NED Reference Frame](../assets/iss_observation/ned.pdf)

If the angle ``\lambda`` is greater than zero, we can theoretically observe the object we
are analyzing because it is above the horizon at the desired location. We can compute this
angle using:

```math
\begin{equation*}
  \lambda = atan\left(-\frac{R_{NED,z}}{\sqrt{R_{NED,x}^2 + R_{NED,y}^2}}\right)~,
\end{equation*}
```

where ``R_{NED,i}`` it the ``i``-axis component of the object position vector ``\vec{R}``
represented in the NED reference frame.

## Algorithm

We need to perform the following tasks to check when ISS enters our field of view:

1. Obtain the ISS position during the desired period;
2. Convert the position vector to the NED reference frame;
3. Obtain the elevation for each instant; and
4. Check when the elevation is greater than zero.

## Code

Before starting, let's load all the packages in the **SatelliteToolbox.jl** ecosystem:

```julia-repl
julia> using SatelliteToolbox
```

We first need to obtain the mean elements of the ISS to propagate its orbit. We can do this
by using the Celestrak TLE fetcher. The following code creates the fetcher and downloads the
latest available ISS TLE:

```julia-repl
julia> tles = fetch_tles(f; satellite_name = "ISS (ZARYA)")
[ Info: Fetch TLEs from Celestrak using satellite name: "ISS (ZARYA)" ...
1-element Vector{TLE}:
 TLE: ISS (ZARYA) (Epoch = 2023-07-02T16:14:23.672)

julia> iss_tle = tles[1]
TLE:
                      Name : ISS (ZARYA)
          Satellite number : 25544
  International designator : 98067A
        Epoch (Year / Day) : 23 / 183.67666287 (2023-07-02T16:14:23.672)
        Element set number : 999
              Eccentricity :   0.00045810
               Inclination :  51.64230000 deg
                      RAAN : 251.73890000 deg
       Argument of perigee :  98.25980000 deg
              Mean anomaly :  37.01280000 deg
           Mean motion (n) :  15.50610282 revs / day
         Revolution number : 40419
                        B* :   0.00021501 1 / er
                     ṅ / 2 :   0.00012101 rev / day²
                     n̈ / 6 :            0 rev / day³
```

The mean elements we obtained are encoded in a TLE by NORAD. Hence, we must use the
SGP4/SDP4 algorithm to propagate its orbit. We can initialize the propagator as follows:

```julia-repl
julia> orbp = Propagators.init(Val(:SGP4), iss_tle)
OrbitPropagatorSgp4{Float64, Float64}:
   Propagator name : SGP4 Orbit Propagator
  Propagator epoch : 2023-07-02T16:14:23.672
  Last propagation : 2023-07-02T16:14:23.672
```

Let's say we want to check when the ISS will stay under our field of view within one day
from this TLE epoch. The following code propagates the TLE for one day using a step of one
second:

```julia-repl
julia> ret = Propagators.propagate!.(orbp, 0:1:86400)
86401-element Vector{Tuple{StaticArraysCore.SVector{3, Float64}, StaticArraysCore.SVector{3, Float64}}}:
 ([4.3241138662552e6, 3.658818213831297e6, 3.737410008328138e6], [-1520.7633082033235, 6181.330955217756, -4275.4053760819315])
 ([4.322590349691385e6, 3.6649972144279014e6, 3.733132227568888e6], [-1526.2732556957171, 6176.663987437403, -4280.179483156526])
 ([4.321061324147931e6, 3.67117154413503e6, 3.7288496754281903e6], [-1531.7812633324513, 6171.989143201469, -4284.948124245782])
 ⋮
 ([-4.523980159202273e6, -3.810874556321396e6, -3.3522994708746294e6], [1497.4675552368813, -5876.71613186301, 4670.189355793778])
 ([-4.522479829177213e6, -3.8167488378395094e6, -3.347627153102706e6], [1503.2022238337156, -5871.8808804082, 4674.448700122841])
```

Each element in the returned array is a tuple with two vectors. The first is the satellite
position [m], and the second is the velocity [m / s]. We only need to check its position.
Thus, let's obtain this information for each instant:

```julia-repl
julia> vr_teme = first.(ret)
86401-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [4.3241138662552e6, 3.658818213831297e6, 3.737410008328138e6]
 [4.322590349691385e6, 3.6649972144279014e6, 3.733132227568888e6]
 [4.321061324147931e6, 3.67117154413503e6, 3.7288496754281903e6]
 ⋮
 [-4.523980159202273e6, -3.810874556321396e6, -3.3522994708746294e6]
 [-4.522479829177213e6, -3.8167488378395094e6, -3.347627153102706e6]
```

This step concludes the first step of the algorithm.

The information obtained by the SGP4 is represented in the True-Equator, Mean-Equinox (TEME)
reference frame, which is a quasi-inertial frame. Now, we need to convert it to a frame
fixed on Earth since we need to compute the elevation angle in a specific position at
Earth's surface. We can do this using the function `r_eci_to_ecef`. It returns a matrix that
rotates an Earth-centered inertial (ECI) frame to an Earth-centered, Earth-fixed (ECEF)
frame.  For our simple example, we will use the PEF (pseudo-Earth fixed) frame as the ECEF.
All the vectors returned by the propagator can be converted as follows:

```julia-repl
julia> vr_pef = r_eci_to_ecef.(TEME(), PEF(), Propagators.epoch(orbp) .+ (collect(0:1:86400) ./ 86400)) .* vr_teme
86401-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [-3.1517745650685215e6, -4.706509167226944e6, 3.737410008328138e6]
 [-3.148954811670437e6, -4.711801726204846e6, 3.733132227568888e6]
 [-3.146131833297915e6, -4.717088716679321e6, 3.7288496754281903e6]
 ⋮
 [3.3857219006283223e6, 4.850365818830511e6, -3.3522994708746294e6]
 [3.383108804842527e6, 4.855406297217666e6, -3.347627153102706e6]
```

!!! note
    The code `Propagators.epoch(orbp) .+ (collect(0:1:86400) ./ 86400)` obtains the Julian
    Day [UTC] of each propagation instant.

We can now convert the ECEF vectors to the NED frame using the function `ecef_to_ned`.
However, we must input the geodetic location we are analyzing. Here, we will use the
location of the city of São José dos Campos, SP, Brazil:

- **Latitude**: 23.1791 S.
- **Longitude**: 45.8872 W.
- **Altitude**: 593 m.

```julia-repl
julia> vr_ned = ecef_to_ned.(vr_pef, -23.1791 |> deg2rad, -45.8872 |> deg2rad, 593; translate = true)
86401-element Vector{StaticArraysCore.SVector{3, Float64}}:
 [3.8867951998111717e6, -5.538957086708884e6, 6.7568967896751845e6]
 [3.885130946017213e6, -5.540616594591433e6, 6.749915536859117e6]
 [3.8834616159196747e6, -5.542269910950413e6, 6.742934017679935e6]
 ⋮
 [-3.540243127040135e6, 5.8070592063407935e6, 6.090776675161325e6]
 [-3.5380883263308997e6, 5.808691621810904e6, 6.097614606058563e6]
```

This step concludes the second step of the algorithm.

!!! note
    We must set the keyword `translate` to `true` because we do not want only to rotate the
    reference frame. We also wish the frame origin for the returned vectors to be translated
    from the Earth's center to the city location.

The elevation in the third step can be easily computed using the presented formula:

```julia-repl
julia> vλ = map(v -> atand(-v[3], sqrt(v[1]^2 + v[2]^2)), vr_ned)
86401-element Vector{Float64}:
 -44.95878130917832
 -44.92746141286671
 -44.89614037059145
   ⋮
 -41.8461848156166
 -41.876994453587045
```

Finally, we need to find the indices related to elevations greater than zero:

```julia-repl
julia> ids = findall(>=(0), vλ)
3102-element Vector{Int64}:
  1409
  1410
  1411
     ⋮
 85164
 85165
```

Those are the instants that we will have a direct line of sight to the ISS. We can discover
the related time using:

```julia-repl
julia> using Dates

julia> getindex(Propagators.epoch(orbp) .+ (collect(0:1:86400) ./ 86400) .|> julian2datetime, ids)
3102-element Vector{DateTime}:
 2023-07-02T16:37:51.672
 2023-07-02T16:37:52.672
 2023-07-02T16:37:53.672
 ⋮
 2023-07-03T15:53:46.672
 2023-07-03T15:53:47.672
```
