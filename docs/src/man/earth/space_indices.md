Space indices
=============

The **SatelliteToolbox.jl** can automatically fetch some space indices that are
used in some computations, notably in the earth atmospheric models. First, it is
necessary to initialize the related files, which is done by the function:

```julia
function init_space_indices(...)
```

When called without arguments, it will download all the supported files, if
necessary. For more information about the many configuration options, please,
see the function documentation.

The supported files are:

|      File       | Download frequency |                                                           Information                                                           |
|-----------------|--------------------|---------------------------------------------------------------------------------------------------------------------------------|
| `DTCFILE.TXT`   | Daily              | This file contains the exospheric temperature variation caused by the Dst index. This is used for the JB2008 atmospheric model. |
| `fluxtable.txt` | Daily              | This file contains the F10.7 flux data in different formats.                                                                    |
| `SOLFSMY.TXT`   | Daily              | This files contains the indices necessary for the JB2008 atmospheric model.                                                     |
| WDC Files       | Once / Daily\*     | This set of files contain the Kp and Ap indices.                                                                                |

\*: The WDC files are separated by year. The file related to the current year is
downloaded on a daily-basis and the files related to the previous years are
downloaded only once.

After the initialization of the files, the space indices can be obtained by the
following function:

```julia
function get_space_index(IND, JD::Number, ...)
```

in which `JD` is the Julian Day in which the index will be computed, and `IND`
is the desired space index as described in the following table.

| `IND`       | Space index                                                                                                                  |
|-------------|------------------------------------------------------------------------------------------------------------------------------|
| `F10()`     | 10.7-cm adjusted solar flux                                                                                                  |
| `F10adj()`  | 10.7-cm adjusted solar flux                                                                                                  |
| `F10obs()`  | 10.7-cm observed solar flux                                                                                                  |
| `F10M()`    | Average of 10.7-cm adjusted solar flux                                                                                       |
| `F10Madj()` | Average of 10.7-cm adjusted solar flux                                                                                       |
| `F10Mobs()` | Average of 10.7-cm observed solar flux                                                                                       |
| `Kp()`      | Kp index (daily mean)                                                                                                        |
| `Kp_vect()` | A vector containing the Kp index for the following hours of the day: 0-3h, 3-6h, 6-9h, 9-12h, 12-15h, 15-18h, 18-20h, 20-23h |
| `Ap()`      | Ap index (daily mean)                                                                                                        |
| `Ap_vect()` | A vector containing the Ap index for the following hours of the day: 0-3h, 3-6h, 6-9h, 9-12h, 12-15h, 15-18h, 18-20h, 20-23h |
| `S10()`     | EUV index (26-34 nm) scaled to F10.7                                                                                         |
| `M10()`     | MG2 index scaled to F10.7                                                                                                    |
| `Y10()`     | Solar X-ray & Lya index scaled to F10.7                                                                                      |
| `S81a()`    | EUV 81-day averaged centered index                                                                                           |
| `M81a()`    | MG2 81-day averaged centered index                                                                                           |
| `Y81a()`    | Solar X-ray & Lya 81-day averaged centered index                                                                             |
| `DstΔTc()`  | Exospheric temperature variation due to `Dst`                                                                                |

!!! note

    The index `DstΔTc()` is interpolated to the selected instant of the Julian
    Day, whereas all the other indices are constants within the seleteced day.

!!! note

    For the indices `F10M()`, `F10Madj()`, and `F10Mobs()`, one can use an
    optional keyword `window::Int` that defines the size in days of the moving
    average window. If it is not specified, then it defaults to 81 days.

```julia-repl
julia> init_space_indices()
[ Info: Downloading file 'DTCFILE.TXT' from 'http://sol.spacenvironment.net/jb2008/indices/DTCFILE.TXT'.
[ Info: Downloading file 'fluxtable.txt' from 'ftp://ftp.geolab.nrcan.gc.ca/data/solar_flux/daily_flux_values/fluxtable.txt'.
[ Info: Downloading file 'SOLFSMY.TXT' from 'http://sol.spacenvironment.net/jb2008/indices/SOLFSMY.TXT'.
[ Info: Downloading file 'kp2017.wdc' from 'ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/wdc/kp2017.wdc'.
[ Info: Downloading file 'kp2015.wdc' from 'ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/wdc/kp2015.wdc'.
[ Info: Downloading file 'kp2016.wdc' from 'ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/wdc/kp2016.wdc'.
[ Info: Downloading file 'kp2018.wdc' from 'ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/wdc/kp2018.wdc'.

julia> get_space_index(F10(), date_to_jd(2018, 6, 19, 18, 35, 00))
79.0

julia> get_space_index(F10M(), date_to_jd(2018, 6, 19, 18, 35, 00))
73.47037037037039

julia> get_space_index(F10M(), date_to_jd(2018, 6, 19, 18, 35, 00); window = 51)
74.60196078431372

julia> get_space_index(Ap(), date_to_jd(2018, 6, 19, 18, 35, 00))
5.125
```
