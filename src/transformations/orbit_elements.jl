# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Conversions related to the orbit elements.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export change_oe_frame

"""
    change_oe_frame(a::Number, e::Number, i::Number, Ω::Number, ω::Number, f::Number, conv_args...)
    change_oe_frame(oe::Orbit, conv_args...)

Change the reference frame of orbit elements. The orbit elements can be
specified by `a`, `e`, `i`, `Ω`, `ω`, and `f`, or the structure `oe` (see
[`Orbit`](@ref)). In the latter, the return value type will match the type of `oe`.

The conversion arguments `conv_args` are **the same** arguments that one should
pass to the function `r_eci_to_eci` to convert between the desired frames. For more
information, see the documentation of the function `r_eci_to_eci`.

# Args

* `a`: Semi-major axis [m].
* `e`: Excentricity.
* `i`: Inclination [rad].
* `Ω`: Right-ascension of the ascending node [rad].
* `ω`: Argument of perigee [rad].
* `f`: True anomaly [rad].
* `conv_args...`: Conversion arguments, which are the same arguments that one
                  would pass to the function `r_eci_to_eci` to convert between the
                  desired frames.

* `oe`: An instance of the structure [`Orbit`](@ref) with the orbit elements
        that will be converted [SI units].

# Returns

Using the first signature, this function returns an instance of
[`KeplerianElements`](@ref). If the second signature is used, then the function
return an element with the same type of the input.

# Examples

```julia-repl
julia> eop = get_iers_eop(:IAU1980);

julia> teme_epoch = DatetoJD(2016,6,1,11,0,0);

julia> tod_epoch  = DatetoJD(2016,1,1,0,0,0);

julia> k_teme     = KeplerianElements(teme_epoch,
                                      7130.982e3,
                                      0.001111,
                                      98.405*pi/180,
                                      227.336*pi/180,
                                      90*pi/180,
                                      320*pi/180)
KeplerianElements{Float64}:
           Epoch :    2.45754e6 (2016-06-01T11:00:00)
 Semi-major axis : 7130.98     km
    Eccentricity :    0.001111
     Inclination :   98.405    °
            RAAN :  227.336    °
 Arg. of Perigee :   90.0      °
    True Anomaly :  320.0      °

julia> k_j2000 = change_oe_frame(k_teme, TEME(), J2000(), teme_epoch, eop)
KeplerianElements{Float64}:
           Epoch :    2.45754e6 (2016-06-01T11:00:00)
 Semi-major axis : 7130.98     km
    Eccentricity :    0.001111
     Inclination :   98.3365   °
            RAAN :  227.134    °
 Arg. of Perigee :   90.0604   °
    True Anomaly :  320.0      °

julia> k_tod = change_oe_frame(k_teme, TEME(), teme_epoch, TOD(), tod_epoch, eop)
KeplerianElements{Float64}:
           Epoch :    2.45754e6 (2016-06-01T11:00:00)
 Semi-major axis : 7130.98     km
    Eccentricity :    0.001111
     Inclination :   98.4037   °
            RAAN :  227.331    °
 Arg. of Perigee :   90.0014   °
    True Anomaly :  320.0      °
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

    k = KeplerianElements(0, a, e, i, Ω, ω, f)

    return change_oe_frame(k, conv_args...)
end

function change_oe_frame(oe::T, conv_args...) where T<:Orbit
    # First, we need to convert to state vector.
    sv_o = convert(OrbitStateVector, oe)
    r_o  = sv_o.r
    v_o  = sv_o.v
    a_o  = sv_o.a

    # NOTE: In my benchmarks, the operation with DCMs are faster than with
    # quaternions after the DCM representation was changed to SMatrix.
    D_ECIf_ECIo = r_eci_to_eci(DCM, conv_args...)
    r_f         = D_ECIf_ECIo*r_o
    v_f         = D_ECIf_ECIo*v_o
    a_f         = D_ECIf_ECIo*a_o

    # Create the new state vector that will be converted to an entity with the
    # type `T`.
    sv_f = orbsv(sv_o.t, r_f, v_f, a_f)

    return convert(T, sv_f)
end
