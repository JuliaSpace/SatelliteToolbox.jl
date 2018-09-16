SatelliteToolbox.jl Changelog
=============================

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
