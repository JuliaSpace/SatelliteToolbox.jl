#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Conversions related to the orbit elements.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Schwarz, R (2014). Memorandum No. 2: Cartesian State Vectors to
#       Keplerian Orbit Elements. Available at www.rene-schwarz.com
#
#       https://downloads.rene-schwarz.com/dc/category/18
#       (Accessed on 2017-08-09).
#
#   [2] Vallado, D. A., McClain, W. D (2013). Fundamentals of Astrodynamics
#       and Applications. Microcosm Press.
#
#   [3] Kuga, H. K., Carrara, V., Rao, K. R (2005). Introdução à Mecânica
#       Orbital. 2ª ed. Instituto Nacional de Pesquisas Espaciais.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export change_oe_frame

"""
    change_oe_frame(a::Number, e::Number, i::Number, Ω::Number, ω::Number, f::Number, conv_args...)
    change_oe_frame(oe::Orbit, conv_args...)

Change the reference frame of orbit elements. The orbit elements can be
specified by `a`, `e`, `i`, `Ω`, `ω`, and `f`, or the structure `oe` (see
`Orbit`).

The conversion arguments `conv_args` are **the same** arguments that one should
pass to the function `rECItoECI` to convert between the desired frames. For more
information, see the documentation of the function `rECItoECI`.

# Args

* `a`: Semi-major axis [m].
* `e`: Excentricity.
* `i`: Inclination [rad].
* `Ω`: Right-ascension of the ascending node [rad].
* `ω`: Argument of perigee [rad].
* `f`: True anomaly [rad].
* `conv_args...`: Conversion arguments, which are the same arguments that one
                  would pass to the function `rECItoECI` to convert between the
                  desired frames.

* `oe`: An instance of the structure `Orbit` with the orbit elements that will
        be converted [SI units].

# Returns

An instance of the structure `Orbit` with the Keplerian elements [SI units]
converted to the new frame.

# Examples

```julia-repl
julia> eop = get_iers_eop(:IAU1980);

julia> teme_epoch = DatetoJD(2016,6,1,11,0,0);

julia> tod_epoch  = DatetoJD(2016,1,1,0,0,0);

julia> oe_teme    = Orbit(0,
                          7130.982e3,
                          0.001111,
                          98.405*pi/180,
                          227.336*pi/180,
                          90*pi/180,
                          320*pi/180)

                 Orbit
  =====================================
                  t =        0.0
    Semi-major axis =     7130.9820 km
       Eccentricity =        0.001111
        Inclination =       98.4050 ˚
               RAAN =      227.3360 ˚
    Arg. of Perigee =       90.0000 ˚
       True Anomaly =      320.0000 ˚

julia> oe_j2000 = change_oe_frame(oe_teme, TEME(), J2000(), teme_epoch, eop)

                 Orbit
  ======================================
                  t =        0.0
    Semi-major axis =     7130.9820 km
       Eccentricity =        0.001111
        Inclination =       98.3365 ˚
               RAAN =      227.1345 ˚
    Arg. of Perigee =       90.0604 ˚
       True Anomaly =      320.0000 ˚

julia> oe_tod   = change_oe_frame(oe_teme, TEME(), teme_epoch, TOD(), tod_epoch, eop)

                 Orbit
  ======================================
                  t =        0.0
    Semi-major axis =     7130.9820 km
       Eccentricity =        0.001111
        Inclination =       98.4037 ˚
               RAAN =      227.3306 ˚
    Arg. of Perigee =       90.0014 ˚
       True Anomaly =      320.0000 ˚
```
"""
function change_oe_frame(a::Number,
                         e::Number,
                         i::Number,
                         Ω::Number,
                         ω::Number,
                         f::Number,
                         conv_args...)

    # The approach is to transform the orbit elements to Cartesian
    # representation, convert the frame, and then convert back to orbit
    # elements.
    #
    # NOTE: In my benchmarks, the operation with DCMs are faster than with
    # quaternions after the DCM representation was changed to SMatrix.

    r_o, v_o    = kepler_to_rv(a, e, i, Ω, ω, f)
    D_ECIf_ECIo = rECItoECI(DCM, conv_args...)
    r_f         = D_ECIf_ECIo*r_o
    v_f         = D_ECIf_ECIo*v_o

    return rv_to_kepler(r_f, v_f)
end

function change_oe_frame(oe::Orbit, conv_args...)
    k = convert(KeplerianElements, oe)
    return change_oe_frame(k.a, k.e, k.i, k.Ω, k.ω, k.f, conv_args...)
end
