SatelliteToolbox.jl Changelog
=============================

Version 0.9.4
-------------

- ![Bugfix][badge-bugfix] The conversion from RV to Keplerian elements had a bug
  for Equatorial and elliptical orbits. (Issue [#72][gh-issue-72])

Version 0.9.3
-------------

- ![Enhancement][badge-enhancement] The conversion between ECEF or geocentric
  variables to geodetic variables can now support different ellipsoids. (PR
  [#61][gh-pr-61])
- ![Bugfix][badge-bugfix] The computation of Legendre associated functions were
  having accuracy problems for angles very close to 0. This problem caused
  incorrect geomagnetic field calculations at regions close the geographic poles.

Version 0.9.2
-------------

- ![Bugfix][badge-bugfix] The TLE file is now closed after parsing in `read_tle`
  function, avoiding possible errors when opening many files. (PR
  [#60][gh-pr-60])

Version 0.9.1
-------------

- ![Bugfix][badge-bugfix] The interpolation type returned by Interpolations.jl
  0.13.3 changed. The code inside SatelliteToolbox.jl was updated to be
  compatible with the new API. (Issue [#58][gh-issue-58])
- ![Info][badge-info] The compat bounds were updated.

Version 0.9.0
-------------

- ![BREAKING][badge-breaking] Update the files from which the IERS data is
  parsed. We now use the CSV. Old code that relies on local copies of IERS data
  needs to be updated.
- ![Deprecation][badge-deprecation] ![Enhancement][badge-enhancement] A type
  instability in IERS parsing has been removed. However, this modification
  required a new function API. The old API is now deprecated.

Version 0.8.1
-------------

- ![Bugfix][badge-bugfix] Avoid error when receiving duplicated information in
  `fluxtable` file. (Issue [#53][gh-issue-53])
- ![Info][badge-info] The compat bounds of PrettyTables.jl were updated.

Version 0.8.0
-------------

- ![BREAKING][badge-breaking] ![Enhancement][badge-enhancement] This version
  contains improvements to the algorithm that converts position and velocity to
  mean elements using SGP4 algorithm. However, those modifications changed the
  returned vectors.
- ![BREAKING][badge-breaking] The definition of `_space_indices_itp_constants`
  changed due to an update in Interpolation.jl
- ![Info][badge-info] The compat bounds were updated.

Version 0.7.3
-------------

- ![Bugfix][badge-bugfix] The source code was updated to be compatible with
  PrettyTables.jl v0.10.
- ![Enhancement][badge-enhancement] A type-instability related to EOP data was
  removed, leading to a performance gain.

Version 0.7.2
-------------

- ![Bugfix][badge-bugfix] The function `DatetoJD` was failing in 32-bit
  platforms when calling with a `DateTime` object. (PR [#43][gh-pr-43])
- ![Bugfix][badge-bugfix] Remove a deprecated function in `eclipse_time`.
- ![Enhancement][badge-enhancement] It is now possible to call the IGRF
  functions `igrf` and `igrfd` passing a pre-allocated matrix for the
  computation of the Legendre associated functions. Hence, in this case, no
  allocation will be performed, leading to a huge performance gain.
- ![Info][badge-info] This version support Julia 1.0 and 1.5. The support for
  Julia 1.3 has been dropped.

Version 0.7.1
-------------

- ![Bugfix][badge-bugfix] Some functions related to space indices were
  deprecated. They were replaced by the correct versions.

Version 0.7.0
-------------

- ![BREAKING][badge-breaking] All the deprecated functions in v0.4 and v0.5 were
  removed.
- ![Feature][badge-feature] The IGRF can be computed passing the latitude and
  longitude in degrees by using the function `igrfd`.
- ![Feature][badge-feature] IGRF v13 is now available using the functions
  `igrf13syn` or `igrf`. The old v12 is still available for compatibility
  reasons in the function `igrf12syn`.
- ![Enhancement][badge-enhancement] Following the Julia documentation, all
  functions that uses `Val` for multiple-dispatch now needs a value
  `Val(:symbol)` instead of a type `Val{:symbol}`. All the old signatures still
  work but have been marked as deprecated.
- ![Deprecation][badge-deprecation] All functions that uses `Val{:symbol}` for
  multiple-dispatch are now marked as deprecated in favor of `Val(:symbol)`.

Version 0.6.5
-------------

- ![Bugfix][badge-bugfix] The FTP address of WDC files was changed. (Issue
  [#38][gh-issue-38])
- ![Bugfix][badge-bugfix] A test related to an error thrown by the TLE submodule
  was fixed. (Issue [#35][gh-issue-35])
- ![Feature][badge-feature] It is now possible to compute the gaps between
  accesses to the ground stations.
- ![Feature][badge-feature] Initial version of the algorithm to convert a
  osculating state vectors (position and velocity) to TLE. **This is still in
  alpha stage and should be used with caution.**
- ![Feature][badge-feature] The function `tle_to_str` can be used to convert an
  object of type `TLE` to a string.
- ![Enhancement][badge-enhancement] Improvements in the documentation of
  functions and macros.
- ![Enhancement][badge-enhancement] The minimum required version of
  [OptionalData](https://github.com/helgee/OptionalData.jl) was updated. (PR
  [#36][gh-pr-36])

Version 0.6.4
-------------

- ![Bugfix][badge-bugfix] The download path of the file `fluxtable.txt` was
  changed. (Issue [#34][gh-issue-34])
- ![Enhancement][badge-enhancement] The submodule SatelliteToolboxSGP4 was
  created, which contains all the low-level functions related to the SGP4 orbit
  propagator. (PR [#33][gh-pr-33])
- ![Enhancement][badge-enhancement] The submodule SatelliteToolboxTLE was
  created, which contains all functions related to TLE handling.
- ![Enhancement][badge-enhancement] The TLE parsing algorithm was improved. If
  the epoch year in TLE is higher than 75, then it will be considered in the
  past.
- ![Info][badge-info] This version support Julia 1.0 and 1.3. The support for
  Julia 1.2 has been dropped.

Version 0.6.3
-------------

- ![Bugfix][badge-bugfix] A bug in eclipse time computation was fixed.
- ![Bugfix][badge-bugfix] The function `ground_station_visible` now accepts
  `SVector`. (Issue [#28][gh-issue-28])
- ![Bugfix][badge-bugfix] The eccentricity was not being printed in TLEs.
- ![Deprecation][badge-deprecation] The function that checks if a satellite is
  inside the visibility circle of a ground station has been renamed to
  `ground_station_visible`. The old one is now marked as deprecated.
- ![Deprecation][badge-deprecation] The function that computes the beta angle
  has been renamed to `beta_angle`. The old one is now marked as deprecated.
- ![Deprecation][badge-deprecation] The function that computes the eclipse times
  has been renamed to `eclipse_time_summary`. The old one is now marked as
  deprecated.
- ![Feature][badge-feature] The orbit propagator API now has the function
  `epoch` that returns the current epoch of a propagator.
- ![Feature][badge-feature] The accesses to ground stations can now be computed
  using the function `ground_station_accesses`.
- ![Feature][badge-feature] It is now possible to select what kind of
  perturbations are desired when computing the beta angle in function
  `beta_angle`.
- ![Feature][badge-feature] Add the simplified dipole model to compute the
  geomagnetic field.
- ![Enhancement][badge-enhancement] The function `propagate_to_epoch!` can now
  be used with `OrbitPropagatorJ4`. (PR [#31][gh-pr-31])
- ![Enhancement][badge-enhancement] The function `JDtoDate` now uses the native
  function `julian2datetime` to convert between Date to Julian date represented
  in `DateTime` format.
- ![Enhancement][badge-enhancement] The eclipse time computation algorithm was
  drastically improved by adding an edge find mechanism, leading to an algorithm
  3x faster.
- ![Info][badge-info] This version support Julia 1.0 and 1.2. The support for
  Julia 1.1 has been dropped.

Version 0.6.2
-------------

- ![Bugfix][badge-bugfix] Colors in printing functions were not working.
- ![Bugfix][badge-bugfix] A bug in the conversion `JDtoDate` when seconds are
  0 was fixed. (Issue [#30][gh-issue-30])

Version 0.6.1
-------------

- ![Bugfix][badge-bugfix] The EGM-08 coefficients in J2 and J4 orbit propagator
  algorithm were slightly wrong.
- ![Feature][badge-feature] It was added an algorithm to compute the Moon
  position: `moon_position_i`.
- ![Enhancement][badge-enhancement] The performance of Sun position algorithm
  was increased by 15%.
- ![Enhancement][badge-enhancement] The J2 and J4 orbit propagator algorithms
  now have more built-in coefficients: EGM96, JGM02, and JGM03.
- ![Enhancement][badge-enhancement] The `Project.toml` file was added, the
  `REQUIRE` file was removed, and the registration process is now handled by
  [Registrator.jl](https://github.com/JuliaComputing/Registrator.jl).

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
[gh-issue-28]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/28
[gh-issue-30]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/30
[gh-issue-34]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/34
[gh-issue-35]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/35
[gh-issue-38]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/38
[gh-issue-53]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/53
[gh-issue-58]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/58
[gh-issue-72]: https://github.com/JuliaSpace/SatelliteToolbox.jl/issues/72

[gh-pr-31]: https://github.com/JuliaSpace/SatelliteToolbox.jl/pull/31
[gh-pr-33]: https://github.com/JuliaSpace/SatelliteToolbox.jl/pull/33
[gh-pr-33]: https://github.com/JuliaSpace/SatelliteToolbox.jl/pull/33
[gh-pr-36]: https://github.com/JuliaSpace/SatelliteToolbox.jl/pull/36
[gh-pr-43]: https://github.com/JuliaSpace/SatelliteToolbox.jl/pull/43
[gh-pr-43]: https://github.com/JuliaSpace/SatelliteToolbox.jl/pull/43
[gh-pr-60]: https://github.com/JuliaSpace/SatelliteToolbox.jl/pull/60
[gh-pr-61]: https://github.com/JuliaSpace/SatelliteToolbox.jl/pull/61
