#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   The NRLMSIS-00 empirical atmosphere model was developed by Mike
#   Picone, Alan Hedin, and Doug Drob based on the MSISE90 model.
#
#   The MSISE90 model describes the neutral temperature and densities in Earth's
#   atmosphere from ground to thermospheric heights. Below 72.5 km the model is
#   primarily based on the MAP Handbook (Labitzke et al., 1985) tabulation of
#   zonal average temperature and pressure by Barnett and Corney, which was also
#   used for the CIRA-86. Below 20 km these data were supplemented with averages
#   from the National Meteorological Center (NMC). In addition, pitot tube,
#   falling sphere, and grenade sounder rocket measurements from 1947 to 1972
#   were taken into consideration. Above 72.5 km MSISE-90 is essentially a
#   revised MSIS-86 model taking into account data derived from space shuttle
#   flights and newer incoherent scatter results. For someone interested only in
#   the thermosphere (above 120 km), the author recommends the MSIS-86 model.
#   MSISE is also not the model of preference for specialized tropospheric work.
#   It is rather for studies that reach across several atmospheric boundaries.
#   (quoted from http://nssdc.gsfc.nasa.gov/space/model/atmos/nrlmsise00.html)
#
#   This Julia version of NRLMSISE-00 was converted from the C version
#   implemented and maintained by Dominik Brodowski <devel@brodo.de> and
#   available at http://www.brodo.de/english/pub/nrlmsise/index.html .
#
#   The source code is available at the following git:
#
#       https://git.linta.de/?p=~brodo/nrlmsise-00.git;a=tree
#
#   The conversion also used information available at the FORTRAN source code
#   available at
#
#       https://ccmc.gsfc.nasa.gov/pub/modelweb/atmospheric/msis/nrlmsise00/
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] https://www.brodo.de/space/nrlmsise/index.html
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export conf_nrlmsise00, gtd7, gtd7d, nrlmsise00

################################################################################
#                                  Constants
################################################################################

include("./nrlmsise00_coefs.jl")

################################################################################
#                                     API
################################################################################

"""
    function nrlmsise00(JD::Number, alt::Number, g_lat::Number, g_long::Number [, f107A::Number, f107::Number, ap::Union{Number,AbstractVector}]; output_si::Bool = true, dversion::Bool = true)

**NRLMSISE-00**

Neutral Atmosphere Empirical Model from the surface to lower exosphere.

This routine computes the NRLMSISE-00 outputs (see `NRLMSISE00_Output`) using
the configurations in the structure `nrlmsise00` (see `NRLMSISE00_Structure`).

Notice that the NRLMSISE-00 will be run using the default flags (see
`NRLMSISE00_DEFAULT_FLAGS`). The user can only change the value of
`flags[:output_m_kg]` using the keyword `output_si` to select whether the output
must be converted to SI units. If more control is needed, then the user must
manually call the function `conf_nrlmsise00` and then call `gtd7` or `gtd7d`
with the desired flags.

If the space indices `f107A`, `f107`, and `ap` are missing, then they will be
obtained from the online databases (see `init_space_indices()`).

# Args

* `JD`: Julian Day [UTC].
* `alt`: Altitude [m].
* `g_lat`: Geodetic latitude [rad].
* `g_long`: Geodetic longitude [rad].
* `f107A`: 81 day average of F10.7 flux (centered on day of year `JD`).
* `f107`: Daily F10.7 flux for previous day.
* `ap`: Magnetic index (daily) if it is a number. If it is an array, then see
        **Remarks**.

# Keywords

* `output_si`: (OPTIONAL) If `true`, then the output units will be [m⁻³] for
               species number density and [kg/m⁻³] for the total density.
               Otherwise, the units will be [cm⁻³] and [g/cm⁻³], respectively.
* `dversion`: (OPTIONAL) If `true`, run `gtd7d`. Otherwise, run `gtd7` (see
              **Remarks**).

# Returns

An instance of the structure `NRLMSISE00_Output`. The result in variable
`den_Total` depends on the value of `dversion` (see **Remarks**, **Notes on
input variables**).

# Remarks

1. The densities of `O`, `H`, and `N` are set to `0` below `72.5 km`.
2. The exospheric temperature `T_exo` is set to global average for altitudes
   below `120 km`. The `120 km` gradient is left at global average value for
   altitudes below `72.5 km`.
3. Anomalous oxygen is defined as hot atomic oxygen or ionized oxygen that can
   become appreciable at high altitudes (`> 500 km`) for some ranges of inputs,
   thereby affection drag on satellites and debris. We group these species under
   the term **Anomalous Oxygen**, since their individual variations are not
   presently separable with the drag data used to define this model component.

## AP

If `ap` is a `Vector`, then it must be a vector with 7 dimensions as described
below:

| Index | Description                                                                   |
|-------|:------------------------------------------------------------------------------|
|     1 | Daily AP.                                                                     |
|     2 | 3 hour AP index for current time.                                             |
|     3 | 3 hour AP index for 3 hours before current time.                              |
|     4 | 3 hour AP index for 6 hours before current time.                              |
|     5 | 3 hour AP index for 9 hours before current time.                              |
|     6 | Average of eight 3 hour AP indices from 12 to 33 hours prior to current time. |
|     7 | Average of eight 3 hour AP indices from 36 to 57 hours prior to current time. |

## Notes on input variables

`f107` and `f107A` values used to generate the model correspond to the 10.7 cm
radio flux at the actual distance of the Earth from the Sun rather than the
radio flux at 1 AU. The following site provides both classes of values:

    ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/

`f107`, `f107A`, and `ap` effects are neither large nor well established below
80 km and these parameters should be set to 150, 150, and 4 respectively.

If `dversion` is `true`, then the total mass `den_Total` (see
`NRLMSISE00_Output`) is the sum of the mass densities of the species `He`, `O`,
`N₂`, `O₂`, `Ar`, `H`, and `N`, but **does not** include anomalous oxygen.

If `dversion` is `false`, then total mass `den_Total` (see `NRLMSISE00_Output`)
is the effective total mass density for drag and is the sum of the mass
densities of all species in this model **including** the anomalous oxygen.

"""
function nrlmsise00(JD::Number,
                    alt::Number,
                    g_lat::Number,
                    g_long::Number;
                    output_si::Bool = true,
                    dversion::Bool = true)

    # Get the space indices. If the altitude is lower than 80km, set to default
    # according to the instructions in NRLMSISE-00 source code.
    if alt < 80e3
        f107A = 150.0
        f107  = 150.0
        ap    = 4.0
    else
        # TODO: The online version of NRLMSISE-00 seems to use 90 days, whereas
        # the NRLMSISE-00 source code mentions 81 days.
        f107A = get_space_index(Val{:F10Madj}, JD; window = 90)
        f107  = get_space_index(Val{:F10adj},  JD-1)
        ap    = get_space_index(Val{:Ap},      JD)
    end

    # Call the NRLMSISE-00 model.
    nrlmsise00(JD, alt, g_lat, g_long, f107A, f107, ap; output_si = output_si,
               dversion = dversion)
end

function nrlmsise00(JD::Number,
                    alt::Number,
                    g_lat::Number,
                    g_long::Number,
                    f107A::Number,
                    f107::Number,
                    ap::Union{Number,AbstractVector};
                    output_si::Bool = true,
                    dversion::Bool = true)

    # Constant
    dgtr = 1.74533e-2 # Convert degrees to radians.

    # Convert the Julian Day to Date.
    (Y, M, D, h, m, s) = JDtoDate(JD)

    # Get the number of days since the beginning of the year.
    doy = round(Int,DatetoJD(Y, M, D, 0, 0, 0) - DatetoJD(Y, 1, 1, 0, 0, 0)) + 1

    # Get the number of seconds since the beginning of the day.
    Δds = 3600h + 60m + s

    # Get the local apparent solar time [hours].
    #
    # TODO: To be very precise, I think this should also take into consideration
    # the Equation of Time. However, the online version of NRLMSISE-00 does not
    # use this.
    lst = Δds/3600 + g_long*12/π

    # Create the input structure for NRLMSISE-00 converting the arguments.
    #
    # Notice that, for this algorithm, the year **is not** used.
    nrlmsise00d = conf_nrlmsise00(Y,
                                  doy,
                                  Δds,
                                  alt/1000,
                                  g_lat/dgtr,
                                  g_long/dgtr,
                                  lst,
                                  f107A,
                                  f107,
                                  ap,
                                  NRLMSISE00_Flags(output_m_kg = output_si))

    # Call the NRLMSISE-00 model.
    nrlmsise00_out = (dversion) ? gtd7d(nrlmsise00d) : gtd7(nrlmsise00d)

    nrlmsise00_out
end

################################################################################
#                             Auxiliary Functions
################################################################################

"""
    function conf_nrlmsise00(year::Int, doy::Int, sec::Number, alt::Number, g_lat::Number, g_long::Number, lst::Number, f107A::Number, f107::Number, ap::[Number, AbstractVector], flags::NRLMSISE00_Flags = NRLMSISE00_Flags())

Create the structure with the proper configuration to call the NRLMSISE-00
model.

Notice that the input variables have the same units of the original model.

# Args

* `year`: Year (currently ignored).
* `doy`: Day of year.
* `sec`: Seconds in day [UT].
* `alt`: Altitude [km].
* `g_lat`: Geodetic latitude [deg].
* `g_long`: Geodetic longitude [deg].
* `lst`: Local apparent solar time (hours).
* `f107A`: 81 day average of F10.7 flux (centered on day of year `doy`).
* `f107`: Daily F10.7 flux for previous day.
* `ap`: Magnetic index (daily) if it is a number. If it is an array, then see
        **Remarks**.
* `flags`: (OPTIONAL) An instance of the structure `NRLMSISE00_Flags` with the
            configuration flags for NRLMSISE00. If omitted, then the default
            configurations will be used.

# Returns

An instance of the structure `NRLMSISE00_Structure`.

# Remarks

If `ap` is a `Vector`, then it must be a vector with 7 dimensions as described
below:

| Index | Description                                                                   |
|-------|:------------------------------------------------------------------------------|
|     1 | Daily AP.                                                                     |
|     2 | 3 hour AP index for current time.                                             |
|     3 | 3 hour AP index for 3 hours before current time.                              |
|     4 | 3 hour AP index for 6 hours before current time.                              |
|     5 | 3 hour AP index for 9 hours before current time.                              |
|     6 | Average of eight 3 hour AP indices from 12 to 33 hours prior to current time. |
|     7 | Average of eight 3 hour AP indices from 36 to 57 hours prior to current time. |

## Notes on input variables

UT, Local Time, and Longitude are used independently in the model and are not of
equal importance for every situation. For the most physically realistic
calculation these three variables should be consistent (`lst=sec/3600 +
g_long/15`). The Equation of Time departures from the above formula for
apparent local time can be included if available but are of minor importance.

`f107` and `f107A` values used to generate the model correspond to the 10.7 cm
radio flux at the actual distance of the Earth from the Sun rather than the
radio flux at 1 AU. The following site provides both classes of values:

    ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/

`f107`, `f107A`, and `ap` effects are neither large nor well established below
80 km and these parameters should be set to 150, 150, and 4 respectively.

"""
function conf_nrlmsise00(year::Int,
                         doy::Int,
                         sec::Number,
                         alt::Number,
                         g_lat::Number,
                         g_long::Number,
                         lst::Number,
                         f107A::Number,
                         f107::Number,
                         ap::Number,
                         flags::NRLMSISE00_Flags = NRLMSISE00_Flags())

    # Constants
    # =========

    dgtr = 1.74533e-2 # Convert degrees to radians.
    hr   = 0.2618     # Convert hour angle to radians.

    # Compute auxiliary variables
    # ===========================
    df     = f107  - f107A
    dfa    = f107A - 150
    stloc  = sin(1hr*lst)
    ctloc  = cos(1hr*lst)
    s2tloc = sin(2hr*lst)
    c2tloc = cos(2hr*lst)
    s3tloc = sin(3hr*lst)
    c3tloc = cos(3hr*lst)

    # Compute Legendre polynomials.
    #
    # TODO: Not all coefficients are used. Hence, we could gain performance
    # here.
    plg = legendre(Val{:conv}, pi/2-g_lat*dgtr, 8, false)'

    # Latitude variation of gravity
    # =============================
    #
    # None for flags.time_independent = false.
    xlat      = (!flags.time_independent) ? 45.0 : Float64(g_lat)
    gsurf, re = _glatf(xlat)

    # Create and return the structure
    # ===============================

    # This is necessary to avoid type instability as reported here:
    #
    #   https://github.com/mauro3/Parameters.jl/issues/58
    #
    nrlmsise00d::NRLMSISE00_Structure{Float64} =
        NRLMSISE00_Structure{Float64}(year     = year,
                                      doy      = doy,
                                      sec      = sec,
                                      alt      = alt,
                                      g_lat    = g_lat,
                                      g_long   = g_long,
                                      lst      = lst,
                                      f107A    = f107A,
                                      f107     = f107,
                                      ap       = ap,
                                      flags    = flags,
                                      re       = re,
                                      gsurf    = gsurf,
                                      df       = df,
                                      dfa      = dfa,
                                      plg      = plg,
                                      stloc    = stloc,
                                      ctloc    = ctloc,
                                      s2tloc   = s2tloc,
                                      c2tloc   = c2tloc,
                                      s3tloc   = s3tloc,
                                      c3tloc   = c3tloc)

    nrlmsise00d
end

function conf_nrlmsise00(year::Int,
                         doy::Int,
                         sec::Number,
                         alt::Number,
                         g_lat::Number,
                         g_long::Number,
                         lst::Number,
                         f107A::Number,
                         f107::Number,
                         ap::AbstractVector,
                         flags::NRLMSISE00_Flags = NRLMSISE00_Flags())

    # Constants
    # =========

    dgtr = 1.74533e-2 # Convert degrees to radians.
    hr   = 0.2618     # Convert hour angle to radians.

    # Compute the new flags
    # =====================
    flags.use_ap_array = true

    # Compute auxiliary variables
    # ===========================
    df     = f107  - f107A
    dfa    = f107A - 150
    stloc  = sin(1hr*lst)
    ctloc  = cos(1hr*lst)
    s2tloc = sin(2hr*lst)
    c2tloc = cos(2hr*lst)
    s3tloc = sin(3hr*lst)
    c3tloc = cos(3hr*lst)

    # Compute Legendre polynomials.
    #
    # TODO: Not all coefficients are used. Hence, we could gain performance
    # here.
    plg = legendre(Val{:conv}, pi/2-g_lat*dgtr, 8, false)'

    # Latitude variation of gravity
    # =============================
    #
    # None for flags.time_independent = false.
    xlat      = (!flags.time_independent) ? 45.0 : Float64(g_lat)
    gsurf, re = _glatf(xlat)

    # Create and return the structure
    # ===============================

    # This is necessary to avoid type instability as reported here:
    #
    #   https://github.com/mauro3/Parameters.jl/issues/58
    #
    nrlmsise00d::NRLMSISE00_Structure{Float64} =
        NRLMSISE00_Structure{Float64}(year     = year,
                                      doy      = doy,
                                      sec      = sec,
                                      alt      = alt,
                                      g_lat    = g_lat,
                                      g_long   = g_long,
                                      lst      = lst,
                                      f107A    = f107A,
                                      f107     = f107,
                                      ap_array = ap,
                                      flags    = flags,
                                      re       = re,
                                      gsurf    = gsurf,
                                      df       = df,
                                      dfa      = dfa,
                                      plg      = plg,
                                      stloc    = stloc,
                                      ctloc    = ctloc,
                                      s2tloc   = s2tloc,
                                      c2tloc   = c2tloc,
                                      s3tloc   = s3tloc,
                                      c3tloc   = c3tloc)

    nrlmsise00d
end

################################################################################
#                               Public Functions
################################################################################

"""
    function gtd7(nrlmsise00d::NRLMSISE00_Structure{T}) where T<:Number

**NRLMSISE-00**

Neutral Atmosphere Empirical Model from the surface to lower exosphere.

This routine computes the NRLMSISE-00 outputs (see `NRLMSISE00_Output`) using
the configurations in the structure `nrlmsise00` (see `NRLMSISE00_Structure`).

# Args

* `nrlmsise00d`: An instance of `NRLMSISE00_Structure`.

# Returns

An instance of structure `NRLMSISE00_Output` with the outputs.

In this case, the total mass `den_Total` (see `NRLMSISE00_Output`) is the sum of
the mass densities of the species `He`, `O`, `N₂`, `O₂`, `Ar`, `H`, and `N`, but
**does not** include anomalous oxygen.

# Remarks

1. The densities of `O`, `H`, and `N` are set to `0` below `72.5 km`.
2. The exospheric temperature `T_exo` is set to global average for altitudes
   below `120 km`. The `120 km` gradient is left at global average value for
   altitudes below `72.5 km`.
3. Anomalous oxygen is defined as hot atomic oxygen or ionized oxygen that can
   become appreciable at high altitudes (`> 500 km`) for some ranges of inputs,
   thereby affection drag on satellites and debris. We group these species under
   the term **Anomalous Oxygen**, since their individual variations are not
   presently separable with the drag data used to define this model component.

"""
function gtd7(nrlmsise00d::NRLMSISE00_Structure{T}) where T<:Number

    @unpack_NRLMSISE00_Structure nrlmsise00d

    # Constants
    # =========

    mn2  = 4
    zn2  = SVector{4,T}([72.5,55.0,45.0,32.5])
    mn3  = 5
    zn3  = SVector{5,T}([32.5,20.0,15.0,10.0,0.0])
    zmix = T(62.5)

    # Initialization of variables
    # ===========================
    meso_tn2  = zeros(MVector{4,T})
    meso_tn3  = zeros(MVector{5,T})
    meso_tgn2 = zeros(MVector{2,T})
    meso_tgn3 = zeros(MVector{2,T})

    # Latitude variation of gravity
    # =============================

    xmm = pdm[3,5]

    # Thermosphere / Mesosphere (above zn2[1])
    # ========================================

    (alt < zn2[1]) && (nrlmsise00d.alt = zn2[1])

    out_thermo = gts7(nrlmsise00d)
    nrlmsise00d.alt = alt

    # Unpack the values again because `gts7` may have modified `nrlmsise00d`.
    @unpack_NRLMSISE00_Structure nrlmsise00d

    # If we are above `zn2[1]`, then we do not need to compute anything else.
    (alt >= zn2[1]) && return out_thermo

    # Unpack the output values from the thermospheric portion.
    @unpack_NRLMSISE00_Output out_thermo

    # Check if we must convert to SI.
    dm28m = (flags.output_m_kg) ?  dm28*1.0e6 : dm28

    # Lower Mesosphere / Upper Stratosphere (between `zn3[1]` and `zn2[1]`)
    # =====================================================================

    meso_tgn2[1] = meso_tgn1_2

    meso_tn2[1]  = meso_tn1_5

    meso_tn2[2]  =  pma[1,1]*pavgm[1]/
                    ( 1 - flags.all_tn2_var*_glob7s(pma[1,:], nrlmsise00d) )
    meso_tn2[3]  =  pma[2,1]*pavgm[2]/
                    ( 1 - flags.all_tn2_var*_glob7s(pma[2,:], nrlmsise00d) )
    meso_tn2[4]  =  pma[3,1]*pavgm[3]/
                    ( 1 - flags.all_tn2_var*flags.all_tn3_var*_glob7s(pma[3,:], nrlmsise00d) )
    meso_tn3[1]  = meso_tn2[4]

    meso_tgn2[2] = pavgm[9]*pma[10,1]*
                    ( 1 + flags.all_tn2_var*flags.all_tn3_var*_glob7s(pma[10,:], nrlmsise00d))*
                    meso_tn2[4]^2/(pma[3,1]*pavgm[3])^2

    # Lower Stratosphere and Troposphere (below `zn3[1]`)
    # ===================================================

    if alt < zn3[1]
        meso_tgn3[1] = meso_tgn2[2]
		meso_tn3[2]  = pma[4,1]*pavgm[4]/
                       ( 1 - flags.all_tn3_var*_glob7s(pma[4,:], nrlmsise00d) )
		meso_tn3[3]  = pma[5,1]*pavgm[5]/
                       ( 1 - flags.all_tn3_var*_glob7s(pma[5,:], nrlmsise00d) )
		meso_tn3[4]  = pma[6,1]*pavgm[6]/
                       ( 1 - flags.all_tn3_var*_glob7s(pma[6,:], nrlmsise00d) )
		meso_tn3[5]  = pma[7,1]*pavgm[7]/
                       ( 1 - flags.all_tn3_var*_glob7s(pma[7,:], nrlmsise00d) )
		meso_tgn3[2] = pma[8,1]*pavgm[8]*
                       ( 1 + flags.all_tn3_var*_glob7s(pma[8,:], nrlmsise00d))*
                       meso_tn3[5]*meso_tn3[5]/(pma[7,1]*pavgm[7])^2
    end

    # Linear Transition to Full Mixing Below `zn2[1]`
    # ===============================================

    dmc  = (alt > zmix) ? 1 - (zn2[1]-alt)/(zn2[1] - zmix) : T(0)
    dz28 = den_N2

    # N2 Density
    # ==========

    dmr = den_N2 / dm28m - 1

    den_N2, tz = _densm(re,
                        gsurf,
                        alt,
                        dm28m,
                        xmm,
                        T(0),
                        zn3,
                        meso_tn3,
                        meso_tgn3,
                        zn2,
                        meso_tn2,
                        meso_tgn2)

    den_N2 *= 1 + dmr*dmc

    # He Density
    # ==========

    dmr    = den_He / (dz28 * pdm[1,2]) - 1
    den_He = den_N2 * pdm[1,2] * ( 1 + dmr*dmc )

    # O Density
    # =========

    den_O  = T(0)
    den_aO = T(0)

    # O2 Density
    # ==========

    dmr    = den_O2 / (dz28 * pdm[4,2]) - 1
    den_O2 = den_N2 * pdm[4,2] * ( 1 + dmr*dmc )

    # Ar Density
    # ==========

    dmr    = den_Ar / (dz28 * pdm[5,2]) - 1
    den_Ar = den_N2 * pdm[5,2] * ( 1 + dmr*dmc )

    # H Density
    # =========

    den_H = T(0)

    # N Density
    # =========

    den_N = T(0)

    # Total Mass Density
    # ==================

    den_Total = 1.66e-24( 4den_He +
                         16den_O  +
                         28den_N2 +
                         32den_O2 +
                         40den_Ar +
                           den_H  +
                         14den_N)

    # Check if we must convert the units to SI.
    (flags.output_m_kg) && (den_Total /= 1000)

    # Temperature at Selected Altitude
    # ================================

    dd, T_alt = _densm(re,
                       gsurf,
                       alt,
                       T(1),
                       T(0),
                       tz,
                       zn3,
                       meso_tn3,
                       meso_tgn3,
                       zn2,
                       meso_tn2,
                       meso_tgn2)

    # Create output structure and return.
    #
    # This is necessary to avoid type instability as reported here:
    #
    #   https://github.com/mauro3/Parameters.jl/issues/58
    #

    nrlmsise00_out::NRLMSISE00_Output{T} =
        NRLMSISE00_Output{T}(den_N     = den_N,
                             den_N2    = den_N2,
                             den_O     = den_O,
                             den_aO    = den_aO,
                             den_O2    = den_O2,
                             den_H     = den_H,
                             den_He    = den_He,
                             den_Ar    = den_Ar,
                             den_Total = den_Total,
                             T_exo     = T_exo,
                             T_alt     = T_alt,
                             flags     = flags)

    nrlmsise00_out
end

"""
    function gtd7d(nrlmsise00d::NRLMSISE00_Structure{T}) where T<:Number

**NRLMSISE-00**

Neutral Atmosphere Empirical Model from the surface to lower exosphere.

This routine computes the NRLMSISE-00 outputs (see `NRLMSISE00_Output`) using
the configurations in the structure `nrlmsise00` (see `NRLMSISE00_Structure`).

# Args

* `nrlmsise00d`: An instance of `NRLMSISE00_Structure`.

# Returns

An instance of structure `NRLMSISE00_Output` with the outputs.

In this case, the total mass `den_Total` (see `NRLMSISE00_Output`) is the
effective total mass density for drag and is the sum of the mass densities of
all species in this model **including** the anomalous oxygen.

# Remarks

1. The densities of `O`, `H`, and `N` are set to `0` below `72.5 km`.
2. The exospheric temperature `T_exo` is set to global average for altitudes
   below `120 km`. The `120 km` gradient is left at global average value for
   altitudes below `72.5 km`.
3. Anomalous oxygen is defined as hot atomic oxygen or ionized oxygen that can
   become appreciable at high altitudes (`> 500 km`) for some ranges of inputs,
   thereby affection drag on satellites and debris. We group these species under
   the term **Anomalous Oxygen**, since their individual variations are not
   presently separable with the drag data used to define this model component.

"""
function gtd7d(nrlmsise00d::NRLMSISE00_Structure{T}) where T<:Number
    # Call `gt7d` to compute the NRLMSISE-00 outputs.
    out = gtd7(nrlmsise00d)

    # Update the computation of the total mass density.
    den_Total = 1.66e-24( 4out.den_He +
                         16out.den_O  +
                         28out.den_N2 +
                         32out.den_O2 +
                         40out.den_Ar +
                           out.den_H  +
                         14out.den_N  +
                         16out.den_aO)

    # Check if we must convert the units to SI.
    (nrlmsise00d.flags.output_m_kg) && (den_Total /= 1000)

    # Create the new output and return.
    #
    # This is necessary to avoid type instability as reported here:
    #
    #   https://github.com/mauro3/Parameters.jl/issues/58
    #

    nrlmsise00_out::NRLMSISE00_Output{T} =
        NRLMSISE00_Output(out; den_Total = den_Total)

    nrlmsise00_out
end

################################################################################
#                              Private Functions
################################################################################

include("./nrlmsise00_priv.jl")
