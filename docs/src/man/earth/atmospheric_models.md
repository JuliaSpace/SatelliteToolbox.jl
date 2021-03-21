Earth atmospheric models
========================

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

This package implements natively three atmospheric models:

* Exponential atmospheric model according to [1];
* Jacchia-Roberts 1971;
* [Jacchia-Bowman 2008](http://sol.spacenvironment.net/jb2008/); and
* [NRLMSISE-00](https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php).

## Exponential atmospheric model

This model assumes that the atmospheric density is computed by:

```math
\rho(h) = \rho_0 \cdot exp \left\lbrace - \frac{h - h_0}{H} \right\rbrace~,
```

in which ``\rho_0``, ``h_0``, and ``H`` are parameters obtained from tables.
Reference [1] provides a discretization of those parameters based on the
selected height ``h`` that was obtained after evaluation of some accurate
models.

In this toolbox, the model can be evaluated using the following function:

```julia
function expatmosphere(h::Number)
```

in which `h` is the desired height in meters.

```jldoctest
julia> expatmosphere(700e3)
3.614e-14
```

!!! warning

    Notice that this model does not consider important effects such as the Sun
    activity, the geomagnetic activity, the local time at the desired location,
    etc.  Hence, although this can be used for fast evaluations, the accuracy is
    not good.

## Jacchia-Roberts 1971

This is an analytic atmospheric model based on the Jacchia's 1970 model. It
was published in:

> **Roberts, C. E (1971)**. *An analytic model for upper atmosphere densities
> based upon Jacchia's 1970 models*. **Celestial mechanics**, Vol. 4 (3-4), p.
> 368-377, DOI: 10.1007/BF01231398.

Although it is quite old, this model is still used for some applications, like
computing the estimated reentry time for an object on low Earth orbits.

In this toolbox, the model can be evaluated using the following function:

```julia
function jr1971(JD::Number, glat::Number, glon::Number, h::Number, F10::Number, F10ₐ::Number, Kp::Number)
```

in which:

* `JD`: Julian day.
* `glat`: Geodetic latitude [rad].
* `glon`: Geodetic longitude [rad].
* `h`: Altitude [m].
* `F10`: 10.7-cm solar flux \[10⁻²² W/(M² Hz)].
* `F10ₐ`: 10.7-cm averaged solar flux, 81-day centered on input time.
* `Kp`: Kp geomagnetic index (with a delay of 3 hours).

Unfortunately, it does not support fetching the space indices automatically yet.

This function returns an object of type `JR1971_Output` that contains the
following fields:

* `nN2`: Number density of N₂ [1/m³].
* `nO2`: Number density of O₂ [1/m³].
* `nO`: Number density of O [1/m³].
* `nAr`: Number density of Ar [1/m³].
* `nHe`: Number density of He [1/m³].
* `nH`: Number density of H [1/m³].
* `rho`: Total density [kg/m³].
* `T_exo`: Exospheric temperature [K].
* `Tz`: Temperature at the selected altitude [K].

```jldoctest
julia> jr1971(DatetoJD(2018, 6, 19, 18, 35, 0), deg2rad(-22), deg2rad(-45), 700e3, 79, 73.5, 1.34)
JR1971_Output{Float64}
  nN2: Float64 2.8434980991303828e7
  nO2: Float64 174222.87498004676
  nO: Float64 1.4139107014677634e11
  nAr: Float64 8.972570981074634
  nHe: Float64 8.773885389988534e11
  nH: Float64 5.702781005702269e10
  rho: Float64 9.684902904883958e-15
  T_exo: Float64 832.0244272210394
  Tz: Float64 832.0204436414625
```

## Jacchia-Bowman 2008

This is an empirical thermospheric density model based on the Jacchia theory. It
was published in:

> **Bowman, B. R., Tobiska, W. K., Marcos, F. A., Huang, C. Y., Lin, C. S.,
> Burke, W. J (2008)**. *A new empirical thermospheric density model JB2008
> using new solar and geomagnetic indices.* **In the proeceeding of the AIAA/AAS
> Astrodynamics Specialist Conference**, Honolulu, Hawaii.

For more information, visit
[http://sol.spacenvironment.net/jb2008](http://sol.spacenvironment.net/jb2008).

In this toolbox, the model can be evaluated using the following functions:

```julia
function jb2008(JD::Number, glat::Number, glon::Number, h::Number)
function jb2008(JD::Number, glat::Number, glon::Number, h::Number, F10::Number, F10ₐ::Number, S10::Number, S10ₐ::Number, M10::Number, M10ₐ::Number, Y10::Number, Y10ₐ::Number, DstΔTc::Number)
```

in which:

* `JD`: Julian day.
* `glat`: Geocentric latitude [rad].
* `glon`: Geocentric longitude [rad].
* `h`: Altitude [m].
* `F10`: 10.7-cm solar flux \[10⁻²² W/(M² Hz)] \(Tabular time 1 day earlier).
* `F10ₐ`: 10.7-cm averaged solar flux, 81-day centered on input time (Tabular
          time 1 day earlier).
* `S10`: EUV index (26-34 nm) scaled to F10.7 (Tabular time 1 day earlier).
* `S10ₐ`: EUV 81-day averaged centered index (Tabular time 1 day earlier).
* `M10`: MG2 index scaled to F10.7 (Tabular time 2 days earlier).
* `M10ₐ`: MG2 81-day averaged centered index (Tabular time 2 days earlier).
* `Y10`: Solar X-ray & Lya index scaled to F10.7 (Tabular time 5 days earlier).
* `Y10ₐ`: Solar X-ray & Lya 81-day averaged centered index (Tabular time 5 days
          earlier).
* `DstΔTc`: Temperature variation related to the Dst.

If the parameters related with the space indices are not provided (first
signature), then they will be automatically obtained. This, however, requires
the initialization of the space indices (see [Space indices]).

These functions returns an object of type `JB2008_Output` that contains the
following fields:

* `nN2`: Number density of N₂ [1/m³].
* `nO2`: Number density of O₂ [1/m³].
* `nO`: Number density of O [1/m³].
* `nAr`: Number density of Ar [1/m³].
* `nHe`: Number density of He [1/m³].
* `nH`: Number density of H [1/m³].
* `rho`: Total density [kg/m³].
* `T_exo`: Exospheric temperature [K].
* `Tz`: Temperature at the selected altitude [K].

```jldoctest
julia> jb2008(DatetoJD(2018, 6, 19, 18, 35, 0), deg2rad(-22), deg2rad(-45), 700e3, 79, 73.5, 55.1, 53.8, 78.9, 73.3, 80.2, 71.7, 50)
JB2008_Output{Float64}
  nN2: Float64 2.6541724729332495e7
  nO2: Float64 193981.21643718384
  nO: Float64 7.674571609797285e10
  nAr: Float64 13.375957587876071
  nHe: Float64 4.642052516165976e11
  nH: Float64 4.072455917445681e10
  rho: Float64 5.193318161219548e-15
  T_exo: Float64 819.7144509572893
  Tz: Float64 826.7686603272322
```

## NRLMSISE-00

The NRLMSIS-00 empirical atmosphere model was developed by Mike Picone, Alan
Hedin, and Doug Drob based on the MSISE90 model:

> **Picone, J. M., Hedin, A. E., Drob, D. P., Aikin, A. C (2002)**. *NRLMSISE-00
> empirical model of the atmosphere: Statistical comparisons and scientific
> issues*. **Journal of Geophysical Research: Space Physics**, Vol. 107 (A12),
> p. SIA 15-1 -- SIA 15-16, DOI: 10.1029/2002JA009430.
 
In this toolbox, the model can be evaluated using the following function:

```julia
function nrlmsise00(JD::Number, alt::Number, g_lat::Number, g_long::Number [, f107A::Number, f107::Number, ap::Union{Number,AbstractVector}]; output_si::Bool = true, dversion::Bool = true)
```

in which:

* `JD`: Julian Day [UTC].
* `alt`: Altitude [m].
* `g_lat`: Geodetic latitude [rad].
* `g_long`: Geodetic longitude [rad].
* `f107A`: 81 day average of F10.7 flux (centered on day of year `doy`).
* `f107`: Daily F10.7 flux for previous day.
* `ap`: Magnetic index.

If the keyword `output_si` is `true`, then the output will be in \[m⁻³] for
species number density and \[kg/m⁻³] for the total density. Otherwise, the units
will be \[cm⁻³] and \[g/cm⁻³], respectively.

The keyword `dversion` can be used to select which algorithm will be used to
compute the model. If it is set to `true`, then it will call the `gtd7d`
function that includes the anomalous oxygen in the total density. Otherwise, the
function `gtd7` will be called and the anomalous oxygen will not be added in the
total density.

The parameter `ap` can be a number or a vector. If it is a number, then it must
be the daily magnetic index. If it is a vector, then it must contain 7 elements
as follows:

| Index | Description                                                                   |
|-------|:------------------------------------------------------------------------------|
|     1 | Daily AP.                                                                     |
|     2 | 3 hour AP index for current time.                                             |
|     3 | 3 hour AP index for 3 hours before current time.                              |
|     4 | 3 hour AP index for 6 hours before current time.                              |
|     5 | 3 hour AP index for 9 hours before current time.                              |
|     6 | Average of eight 3 hour AP indices from 12 to 33 hours prior to current time. |
|     7 | Average of eight 3 hour AP indices from 36 to 57 hours prior to current time. |

If the parameters related with the space indices are not provided (`f107A`,
`f107`, and `ap`), then they will be automatically obtained. This, however,
requires the initialization of the space indices (see [Space indices]).

The function return an object of type `NRLMSISE00_Output` that contains the
following fields:

* `den_N`: Nitrogen number density [U].
* `den_N2`: N₂ number density [U].
* `den_O`: Oxygen number density [U].
* `den_aO`: Anomalous Oxygen number density [U].
* `den_O2`: O₂ number density [U].
* `den_H`: Hydrogen number density [U].
* `den_He`: Helium number density [U].
* `den_Ar`: Argon number density [U].
* `den_Total`: Total mass density \[T/U] \(this value has different meanings for
  routines `gtd7` and `gtd7d`).
* `T_exo`: Exospheric temperature [K].
* `T_alt`: Temperature at the selected altitude [K].
* `flags`: Flags used to compute NRLMSISE-00 model.

Notice that:

* If `flags[:output_m_kg]` is `false`, then [U] is \[cm⁻³] and [T] is \[g/cm⁻³].
* If `flags[:output_m_kg]` is `true`, then [U] is \[m⁻³] and [T] is \[kg/m⁻³].

```jldoctest
julia> nrlmsise00(DatetoJD(2018, 6, 19, 18, 35, 0), 700e3, deg2rad(-22), deg2rad(-45), 73.5, 79, 5.13)
NRLMSISE00_Output{Float64}
  den_N: Float64 5.597834653523333e9
  den_N2: Float64 5.743312510585916e7
  den_O: Float64 1.2705655159941983e11
  den_aO: Float64 2.4185558056141124e9
  den_O2: Float64 340464.9752380828
  den_H: Float64 1.2667781795293002e11
  den_He: Float64 6.248499395447589e11
  den_Ar: Float64 23.18462060029951
  den_Total: Float64 7.930928885098513e-15
  T_exo: Float64 837.4122645268103
  T_alt: Float64 837.4119807046409
  flags: NRLMSISE00_Flags
```

!!! note

    If the user wants more control over the NRLMSISE-00, they can use the
    low-level functions `gtd7` and `gtd7d`, which has the same functionality as
    available in the FORTRAN version of the model. See the documentation of the
    functions for more information.

## References

[1] **Vallado, D. A., McClain, W. D (2013).** *Fundamentals of astrodynamics and
applications*. Hawthorne, CA: Microcosm Press.
