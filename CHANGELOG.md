SatelliteToolbox.jl Changelog
=============================

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
