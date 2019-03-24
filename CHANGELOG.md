SatelliteToolbox.jl Changelog
=============================

Version 0.6.0
-------------

- ![BREAKING][badge-breaking] ![Enhancement][badge-enhancement] The API
  initialization functions of all orbit propagators have changed! Please, see
  the documentation for more details. The analytical propagators Two Body, J2,
  and J4 are now initialized using the mean orbit elements or TLEs, whereas the
  analytical orbit propagator SGP4 can now only be initialized with TLEs.
- ![BREAKING][badge-breaking] ![Enhancement][badge-enhancement] Some API
  functions of `[d]legendre` were removed to simplify the code.
- ![BREAKING][badge-breaking] All the orbital elements in the `Orbit` structure
  must be of the same type now. This does not affect the initialization function
  `Orbit`, which automatically handles the conversion.
- ![BREAKING][badge-breaking] ![Enhancement][badge-enhancement] The user can now
  select the maximum order that will be used in gravitational models. However,
  this is a breaking change because selecting the maximum degree as 0 will make
  the methods use only the first term, leading to the same result as selecting
  it as 1. To use all the terms, the maximum degree must be set to -1. Notice
  that this is still the default behavior if the maximum degree is omitted.
- ![Deprecation][badge-deprecation] All deprecations related to the reference
  frame transformations were removed.
- ![Bugfix][badge-bugfix] The `Orbit` structure created when SPG4 were being
  initialized was not converting the semi-major axis from km to m. Notice that
  this was not affecting the propagation.
- ![Bugfix][badge-bugfix] When using the API to initialize the J2 orbit
  propagator, the parameter `n_0` is the mean motion, *i.e.* the mean angular
  velocity between two consecutive passages to the perigee. The algorithm
  obtained from Vallado's book, on the other hand, seems to consider that `n_0`
  is the rate of change of the mean anomaly. Hence, we must subtract the
  time-derivative of the perigee in order to make the algorithms compatible.
- ![Bugfix][badge-bugfix] Treat special cases when converting orbital elements
  to state vector (position and velocity). (Issue [#25][gh-issue-25])
- ![Feature][badge-feature] The package now supports the entire IAU-2006/2010
  theory (CIO approach).
- ![Feature][badge-feature] It is now possible to compute the satellite ground
  trace using the function `ground_trace`.
- ![Feature][badge-feature] J4 orbit propagator. This new propagator considers
  the perturbation terms J2, J2², and J4.
- ![Feature][badge-feature] Many general orbit analysis functions now support
  perturbation terms up to J4.
- ![Feature][badge-feature] The function `angvel_to_a` can be used to convert
  the orbit mean motion into the semi-major axis.
- ![Feature][badge-feature] The satellite state vector can now be represented by
  the structure `SatelliteStateVector`, which can be initialized by the function
  `satsv`.
- ![Feature][badge-feature] The satellite state vector can be convert between
  all the supported reference frames using `svECItoECI`, `svECItoECEF`,
  `svECEFtoECI`, and `svECEFtoECEF`. Notice that this conversion take into
  account the Earth rotation rate when applicable.
- ![Feature][badge-feature] The user can now select the maximum order that will
  be used when computing the Legendre associated functions or they derivatives.
- ![Enhancement][badge-enhancement] The function `kepler_to_rv` can now receive
  an instance of `Orbit` as input.
- ![Enhancement][badge-enhancement] All colors are now handled by
  [Crayons.jl](https://github.com/KristofferC/Crayons.jl).
- ![Enhancement][badge-enhancement] Many conversion of reference frames using
  the API can now be called without the EOP data, leading to the default
  IAU-76/FK5 theory. (Issue [#12][gh-issue-12])
- ![Enhancement][badge-enhancement] Improve printing functions of `Orbit` and
  `TLE`. (Issue [#21][gh-issue-21])

Version 0.5.1
-------------

- ![Bugfix][badge-bugfix] The gravity model computations were neglecting the
  input variable `n_max` and were using always the maximum available degree.
  (Issue [#22][gh-issue-22])
- ![Bugfix][badge-bugfix] `dlegendre` function was accessing unallocated data if
  the maximum degree is 1.
- ![Feature][badge-feature] New function in the orbit propagators API:
  `propagate_to_epoch!`. (Issue [#20][gh-issue-20])

Version 0.5.0
-------------

- ![BREAKING][badge-breaking] ![Enhancement][badge-enhancement] The NRLMSISE-00
  configuration is not handle by a `Dict` anymore. The configuration is now
  based on a new structure called `NRLMSISE00_Flags`. This reduced the execution
  time of the algorithm by 50%.
- ![BREAKING][badge-breaking] ![Enhancement][badge-enhancement] The macro
  `@check_orbit` that verifies whether an orbit is valid now returns a boolean
  instead of throwing an exception. This provided a huge performance gain in the
  functions that use it.
- ![BREAKING][badge-breaking] ![Enhancement][badge-enhancement] The function
  `GeodetictoECEF` now returns an `SVector`. This yielded a performance gain of
  27%.
- ![BREAKING][badge-breaking] The `igrf12` function now always returns an
  `SVector`.
- ![BREAKING][badge-breaking] All the functions related to the Sun
  (`sun_position_i` and `sun_velocity_i`) now returns an `SVector`.
- ![BREAKING][badge-breaking] All deprecated functions that were defined prior
  to this version were removed.
- ![Deprecation][badge-deprecation] ![Enhancement][badge-enhancement] The frame
  transformations API were simplified by removing the variable that specified
  the model. The IAU2000A theory will be initially implemented using the CIO
  approach. Hence, the model that must be used can be inferred only by looking
  the selected frames and the EOP data. The old style functions are now
  deprecated. (Issue #18)
- ![Deprecation][badge-deprecation] The function `satellite_orbit_compute_f` is
  now deprecated in favor of `M_to_f`. (Issue #16)
- ![Bugfix][badge-bugfix] The `Orbit` structure returned by `propagate!`
  function was not being copied before the assignment. Hence, if an array of
  instants was passed, then all the returned values would have the same `Orbit`
  structure.
- ![Bugfix][badge-bugfix] The low level functions related to the J2 orbit
  propagator were not being exported.
- ![Bugfix][badge-bugfix] The conversion from ECEF to Geodetic had an bug due to
  the singularity in the poles.
- ![Bugfix][badge-bugfix] The conversion from Julian Day to `DateTime` by the
  function `JDtoDate` was fixed.
- ![Bugfix][badge-bugfix] The function `JDtoDate` was not allowing the month to
  be 12. Hence, when the month was December, then it was returning 0. (Issue
  #19).
- ![Bugfix][badge-bugfix] The IAU2000A EOP data was not being parsed correctly.
- ![Bugfix][badge-bugfix] The Sun position algorithm was computing a slightly
  wrong value.
- ![Feature][badge-feature] The initial version of the package documentation is
  now available.
- ![Feature][badge-feature] The EOP data is now downloaded by the package
  RemoteFiles.jl. This modification was required because the HTTP.jl package was
  failing to download EOP data for some reason, since the request was returning
  the 403 error code. This is a temporary workaround until the issue #2 is
  fixed.
- ![Feature][badge-feature] A new generic parser to ICGEM files was added.
  Hence, the user can now load any ICGEM file to compute the gravitational
  force. (Issue #17)
- ![Feature][badge-feature] It was added an algorithm to compute the Sun
  velocity vector.
- ![Feature][badge-feature] The exponential atmospheric model was added.
- ![Feature][badge-feature] The Jacchia-Roberts 1971 atmospheric model was
  added.
- ![Feature][badge-feature] The Jacchia-Bowman 2008 atmospheric model was added.
- ![Feature][badge-feature] It was added the support to automatically fetch the
  necessary files to compute many space indices, like F10.7, Kp, Ap, etc.
- ![Enhancement][badge-enhancement] The macro `@evalpoly` is now used every time
  a polynomial must be evaluated, which provides a small performance gain and
  also uses a more stable algorithm.
- ![Enhancement][badge-enhancement] Huge performance increase. In the following,
  it is presented the reduction of the computational time in some functions:
    - `legendre (conv)`: 40.6% faster.
    - `dlegendre (conv)`: 42.6% faster.
    - `rECEFtoECI (ITRF -> GCRF)`: 99.5% faster.
    - `rECItoECEF (GCRF -> ITRF)`: 99.5% faster.
    - `rECEFtoECEF (ITRF -> PEF)`: 99.8% faster.
    - `rECItoECI (MOD -> GCRF)`: 99.8% faster.
    - `GeodetictoECEF`: 71.6% faster.
    - `Two body orbit propagator`: 99.8% faster.
    - `J2 orbit propagator`: 99.8% faster.
    - `SGP4 orbit propagator`: 8.2% faster.
    - `compute_g (full EGM96)`: 18.3% faster.
    - `nrlmsise00`: 76.3% faster.
    - `igrf12 (geocentric)`: 11.2% faster.
    - `igrf12 (geodetic)`: 98.6% faster.
    - `sun_position_i`: 20.2% faster.
    - `angvel`: 59.3% faster.
    - `compute_ss_orbit_by_ang_vel`: 90.8% faster.
    - `list_ss_orbits_by_rep_period`: 79.4% faster.
- ![Enhancement][badge-enhancement] The space indices required by the
  atmospheric models can now be automatically fetched.
- ![Enhancement][badge-enhancement] Improvements in the function documentation.

Version 0.4.0
-------------

- Fix warnings due to new syntax of the package **Interpolations.jl**.
- Fix remaining compatibilities issues with Julia v1.0.
- Changes in orbit propagators:
    * `propagate!` can now receive one single epoch.
    * All functions related to the orbit propagators uses `SVector`. This led to
      a huge performance gain in API function `propagate!`. However, this
      **can break existing code**, because the array returned by `propagate!` is
      now an array of `SVector{3,T}` instead of an array of `Vector{T}`.

Version 0.3.2
-------------

- Fix compatibility issues with Parameters.jl v0.10.
- Improve performance of `propagate!` functions by fixing type-stability issues.
- `ECEFtoGeodetic` now accepts `AbstractVector` as input.

Version 0.3.1
-------------

- Fix bug when using the function `igrf12` for dates after 2020.
- Add option `shown_warns` to suppress warning in IGRF.

Version 0.3.0
-------------

- Full support for `Julia 0.7` and `Julia 1.0`.
    * The support for `Julia <= 0.6` is dropped in this version. Hence,
      SatelliteToolbox.jl **will not** work with those versions anymore. If it
      is necessary to use `Julia <= 0.6`, then you must stick with
      `SatelliteToolbox.jl <= 0.2.0`.

Version 0.2.0
-------------

- New model:
    * NRLMSISE-00 Atmosphere Model.

- Earth gravity model:
    * New algorithm to compute the gravity model.
    * Add EGM96, JGM2, and JGM3 to the package.
    * Remove the partial EGM2008 from the package.
    * Increase the performance of the function `compute_U`.
    * Increase the performance of the function `compute_g`.
    * Add tests.

- Reference frame transformations:
    * Improvements in IAU76/FK5 conversion.
    * Add support for TEME reference frame.
    * Add more tests to increase the code coverage.
    * Add the function `rECItoECI`.
    * Add the function `rECEFtoECEF`.

- Orbit functions:
    * Add more tests to increase the code coverage.
    * Add function `change_oe_frame`.
    * Add function to read TLE from strings.
    * Add macros `@tle_str` and `@tlenc_str` to parse a TLE from a string.

- Orbit propagators:
    * General performance improvements in SGP4.
    * The epoch of all propagators are now assumed to be the Julian Day.
    * Add Deep Space computations to the SGP4.
    * Add the tests available in the literature that cover all SGP4/SDP4.

Version 0.1.1
-------------

- Fix some deprecations warnings in julia v0.7.

- Bug fixes:
    * Fix IGRF function `igrf12` when date is between 1995 and 2000.

- Tests:
    * Add tests related to IGRF functions.
    * Fix minor problem in IAU-76/FK5 tests.
    * Fix tests in julia v0.7.

Version 0.1.0
-------------

- Initial version.
    * This version was based on the old package **SatToolbox.jl v0.3.7** that
      was renamed to **SatelliteToolbox** to be submitted to julia METADATA
      repo.

[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/Deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-enhancement]: https://img.shields.io/badge/Enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/Bugfix-purple.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

[gh-issue-12]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/12
[gh-issue-20]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/20
[gh-issue-21]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/21
[gh-issue-22]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/22
[gh-issue-25]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/25
