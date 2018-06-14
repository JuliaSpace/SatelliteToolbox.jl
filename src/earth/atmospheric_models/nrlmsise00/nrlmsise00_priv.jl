#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Private functions for NRLMSISE-00 atmosphere model.
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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-06-10: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

"""
    @inline function _ccor(alt::T, r::T, h1::T, zh::T) where T<:Number

Chemistry / Dissociation correction for MSIS models.

# Args

* `alt`: Altitude.
* `r`: Target ratio.
* `h1`: Transition scale length.
* `zh`: Altitude of `1/2 r`.

# Returns

The chemistry / dissociation correction.

"""
@inline function _ccor(alt::T, r::T, h1::T, zh::T) where T<:Number
    e = (alt - zh) / h1

    (e > +70) && return exp(T(0))
    (e < -70) && return exp(r)

    exp( r / (1 + exp(e)) )
end

"""
    @inline function _ccor2(alt::T, r::T, h1::T, zh::T, h2::T) where T<:Number

Chemistry / Dissociation correction for MSIS models.

# Args

* `alt`: Altitude.
* `r`: Target ration.
* `h1`: Transition scale length.
* `zh`: Altitude of `1/2 r`.
* `h2`: Transition scale length 2.

# Returns

The chemistry / dissociation correction.

"""

@inline function _ccor2(alt::T, r::T, h1::T, zh::T, h2::T) where T<:Number

    e1 = (alt - zh) / h1
    e2 = (alt - zh) / h2

    ( (e1 > +70) || (e2 > +70) ) && return exp(T(0))
    ( (e1 < -70) && (e2 < -70) ) && return exp(r)

    exp(r / ( 1 + ( exp(e1) + exp(e2) )/2 ) )
end

"""
    function _densm(re::T, gsurf::T, alt::T, d0::T, xm::T, tz::T, zn3::StaticVector{N3,T}, tn3::AbstractVector{T}, tgn3::AbstractVector{T}, zn2::StaticVector{N2,T}, tn2::AbstractVector{T}, tgn2::AbstractVector{T}) where T<:Number where N2 where N3

Compute the temperature and density profiles for lower atmosphere.

# Returns

* The density.
* The temperature.

"""
function _densm(re::T,
                gsurf::T,
                alt::T,
                d0::T,
                xm::T,
                tz::T,
                zn3::StaticVector{N3,T},
                tn3::AbstractVector{T},
                tgn3::AbstractVector{T},
                zn2::StaticVector{N2,T},
                tn2::AbstractVector{T},
                tgn2::AbstractVector{T}) where T<:Number where N2 where N3

    # Constants
    # =========

    rgas = T(831.4)

    # Initialization of Variables
    # ===========================

    densm_tmp = d0

    if alt > zn2[1]
        (xm == 0) && (densm_tmp = tz)
        return densm_tmp, tz
    end

    #                Stratosphere / Mesosphere Temperature
    # ==========================================================================

    z     = (alt > zn2[N2]) ? alt : zn2[N2]
    z1    = zn2[1]
    z2    = zn2[N2]
    t1    = tn2[1]
    t2    = tn2[N2]
    zg    = _zeta(re, z, z1)
    zgdif = _zeta(re, z2, z1)

    # Set up spline nodes.
    xs2 = zeros(MVector{N2, T})
    ys2 = zeros(MVector{N2, T})

    for k = 1:N2
        xs2[k] = _zeta(re, zn2[k], z1)/zgdif
        ys2[k] = 1 / tn2[k]
    end

    yd1 = -tgn2[1] / (t1*t1) * zgdif
    yd2 = -tgn2[2] / (t2*t2) * zgdif * ( (re+z2)/(re+z1) )^2

    # Calculate spline coefficients.
    y2out = _spline(xs2, ys2, yd1, yd2)
    x     = zg/zgdif
    y     = _splint(xs2, ys2, y2out, x)

    # Temperature at altitude.
    tz = 1 / y

    if xm != 0
        # Calculate stratosphere / mesosphere density.
        glb  = gsurf / (1 + z1/re)^2
        gamm = xm * glb * zgdif / rgas

        # Integrate temperature profile.
        yi   = _splini(xs2, ys2, y2out, x)
        expl = gamm*yi;
        (expl > 50) && (expl = T(50))

        # Density at altitude.
        densm_tmp *= (t1 / tz) * exp(-expl)
    end

    if alt > zn3[1]
        (xm == 0) && (densm_tmp = tz)
        return densm_tmp, tz
    end

    #               Troposphere / stratosphere temperature.
    # ==========================================================================
    z     = alt
    z1    = zn3[1]
    z2    = zn3[N3]
    t1    = tn3[1]
    t2    = tn3[N3]
    zg    = _zeta(re, z, z1)
    zgdif = _zeta(re, z2, z1)

    # Set up spline nodes.
    xs3 = zeros(MVector{N3, T})
    ys3 = zeros(MVector{N3, T})

    for k = 1:N3
        xs3[k] = _zeta(re, zn3[k], z1) / zgdif
        ys3[k] = 1 / tn3[k]
    end

    yd1 = -tgn3[1] / (t1*t1) * zgdif
    yd2 = -tgn3[2] / (t2*t2) * zgdif * ( (re+z2)/(re+z1) )^2

    # Calculate spline coefficients.
    y2out = _spline(xs3, ys3, yd1, yd2)
    x     = zg/zgdif
    y     = _splint(xs3, ys3, y2out, x)

    # Temperature at altitude.
    tz = 1 / y

    if xm != 0
        # Calculate tropospheric / stratosphere density.
        glb  = gsurf / (1 + z1/re)^2
        gamm = xm * glb * zgdif / rgas;

        # Integrate temperature profile.
        yi   = _splini(xs3, ys3, y2out, x)
        expl = gamm*yi;
        (expl > 50) && (expl = T(50))

        # Density at altitude.
        densm_tmp *= (t1 / tz) * exp(-expl)
    end

    (xm == 0) && (densm_tmp = tz)

    densm_tmp, tz
end

"""
    function _densu(re::T, gsurf::T, alt::T, dlb::T, tinf::T, tlb::T, xm::T, alpha::T, zlb::T, s2::T, zn1::StaticVector{N,T}, tn1::AbstractVector{T}, tgn1::AbstractVector{T}) where T<:Number where N

Compute the temperature and density profiles for MSIS models.

This algorithm uses new lower thermo polynomial.

# Returns

* The density.
* The temperature.

"""

function _densu(re::T,
                gsurf::T,
                alt::T,
                dlb::T,
                tinf::T,
                tlb::T,
                xm::T,
                alpha::T,
                zlb::T,
                s2::T,
                zn1::StaticVector{N,T},
                tn1::AbstractVector{T},
                tgn1::AbstractVector{T}) where T<:Number where N

    x         = T(0)
    rgas      = T(831.4)
    z1        = T(0)
    t1        = T(0)
    zgdif     = T(0)
    mn        = 1
    xs        = zeros(MVector{N,T})
    ys        = zeros(MVector{N,T})
    y2out     = zeros(MVector{N,T})

    # Joining altitudes of Bates and spline.
    za = zn1[1]
    z  = (alt > zn1[1]) ? alt : za

    # Geopotential altitude difference from ZLB.
    zg2 = _zeta(re, z, zlb)

    # Bates temperature.
    tt        = tinf - (tinf - tlb) * exp(-s2*zg2)
    ta        = tt
    tz        = tt
    densu_tmp = tz

    if alt < za
        # Compute the temperature below ZA temperature gradient at ZA from Bates
        # profile.
        dta = (tinf - ta) * s2 * ( (re+zlb)/(re+za) )^2

        tgn1[1] = dta
        tn1[1]  = ta
        z       = (alt > zn1[N]) ? alt : zn1[N]
        z1      = zn1[1]
        z2      = zn1[N]
        t1      = tn1[1]
        t2      = tn1[N]

        # Geopotential difference from z1.
        zg    = _zeta(re,  z, z1)
        zgdif = _zeta(re, z2, z1)

        # Set up spline nodes.
        for k = 1:N
            xs[k] = _zeta(re, zn1[k], z1) / zgdif
            ys[k] = 1 / tn1[k]
        end

        # End node derivatives.
        yd1 = -tgn1[1] / (t1*t1) * zgdif
        yd2 = -tgn1[2] / (t2*t2) * zgdif * ( (re+z2)/(re+z1) )^2

        # Compute spline coefficients.
        y2out = _spline(xs, ys, yd1, yd2)
        x     = zg / zgdif
        y     = _splint(xs, ys, y2out, x)

        # Temperature at altitude.
        tz         = 1 / y
        densu_tmp = tz
    end

    (xm == 0) && return densu_tmp, tz

    # Calculate density above za.
    glb   = gsurf / (1 + zlb/re)^2
    gamma = xm * glb / (s2 * rgas * tinf)
    expl  = exp(-s2 * gamma * zg2)

    ( (expl > 50) || (tt <= 0) ) && (expl = T(50))

    # Density at altitude.
    densu_tmp = dlb * (tlb/tt)^(1 + alpha + gamma) * expl

    (alt >= za) && return densu_tmp, tz

    # Compute density below za.
    glb  = gsurf / (1 + z1/re)^2
    gamm = xm * glb * zgdif / rgas

    # Integrate spline temperatures.
    yi   = _splini(xs, ys, y2out, x)
    expl = gamm * yi

    ( (expl > 50) || (tz <= 0) ) && (expl = T(50))

    # Density at altitude.
    densu_tmp *= (t1 / tz)^(1 + alpha) * exp(-expl)

    densu_tmp, tz
end

"""
    @inline function _dnet(dd::T, dm::T, zhm::T, xmm::T, xm::T) where T<:Number

Turbopause correction for MSIS models.

# Args

* `dd`: Diffusive density.
* `dm`: Full mixed density.
* `zhm`: Transition scale length.
* `xmm`: Full mixed molecular weight.
* `xm`: Species molecular weight.

# Returns

The combined density.

"""
@inline function _dnet(dd::T, dm::T, zhm::T, xmm::T, xm::T) where T<:Number

    a  = zhm / (xmm-xm)

    if !( (dm>0) && (dd>0) )
        warn("dnet log error $dm $dd $xm")

        ( (dd == 0) && (dm == 0) ) && (dd = T(1))
        ( dm == 0 ) && return dd
        ( dd == 0 ) && return dm
    end

    ylog = a * log(dm/dd)

    (ylog < -10) && return dd
    (ylog > +10) && return dm

    dd*( 1 + exp(ylog) )^(1/a)
end

@inline function _glatf(lat::T) where T<:Number
    dgtr = T(1.74533e-2)
    c2   = cos(2dgtr*lat)
    gv   = T(980.616)*(1 - T(0.0026373)*c2)
    reff = 2gv / (T(3.085462e-6) + T(2.27e-9)*c2) * T(1e-5)

    (gv, reff)
end

"""
    function _globe7!(p::AbstractVector{T}, nrlmsise00d::NRLMSISE00_Structure{T}) where T<:Number

Compute G(L) function.

Notice that the parameters `apt` and `apdf` of structure `nrlmsise00d` are
modified.

# Args

* `p`: Vector with the coefficients.
* `nrlmsise00d`: NRLMSISE-00 structure (see `NRLMSISE00_Structure`).

# Returns

The temperature (?).

"""
function _globe7!(p::AbstractVector{T},
                  nrlmsise00d::NRLMSISE00_Structure{T}) where T<:Number
    T<:Number

    @unpack_NRLMSISE00_Structure nrlmsise00d

    # Constants
    # =========
    dgtr = T(1.74533e-2)
    dr   = T(1.72142e-2)
    hr   = T(0.2618)
    sr   = T(7.2722e-5)

    # Initialization of variables
    # ===========================
    t    = zeros(T,15)
    tloc = lst

    cd32 = cos( 1dr*(doy-p[32]) )
    cd18 = cos( 2dr*(doy-p[18]) )
    cd14 = cos( 1dr*(doy-p[14]) )
    cd39 = cos( 2dr*(doy-p[39]) )

    # F10.7 Effect
    # ============

    t[1] = p[20]*df*(1 + p[60]*dfa) + p[21]*df^2 + p[22]*dfa + p[30]*dfa^2
    f1   = 1 + (p[48]*dfa + p[20]*df + p[21]*df^2)*flags[:F107_Mean]
    f2   = 1 + (p[50]*dfa + p[20]*df + p[21]*df^2)*flags[:F107_Mean]

    # Time Independent
    # ================

    t[2] =  p[2]*plg[1,3] + p[3]*plg[1,5] + p[23]*plg[1,7] + p[27]*plg[1,2] +
           p[15]*plg[1,3]*dfa*flags[:F107_Mean]

    # Symmetrical Annual
    # ==================

    t[3] = p[19]*cd32

    # Symmetrical Semiannual
    # ======================

    t[4] = (p[16] + p[17]*plg[1,3])*cd18

    # Asymmetrical Annual
    # ===================

    t[5] = f1*( p[10]*plg[1,2] + p[11]*plg[1,4] )*cd14

    # Asymmetrical Semiannual
    # =======================

    t[6] = p[38]*plg[1,2]*cd39

    # Diurnal
    # =======

    if flags[:diurnal]
        t71  = ( p[12]*plg[2,3] )*cd14*flags[:asym_annual]
        t72  = ( p[13]*plg[2,3] )*cd14*flags[:asym_annual]

        t[7] = f2*( (p[4]*plg[2,2] + p[5]*plg[2,4] + p[28]*plg[2,6] + t71 ) *
                   ctloc + ( p[7]*plg[2,2] +
                             p[8]*plg[2,4] +
                            p[29]*plg[2,6] + t72)*stloc )
    end

    # Semidiurnal
    # ===========

    if flags[:semidiurnal]
        t81  = ( p[24]*plg[3,4] + p[36]*plg[3,6])*cd14*flags[:asym_annual]
        t82  = ( p[34]*plg[3,4] + p[37]*plg[3,6])*cd14*flags[:asym_annual]

        t[8] = f2*( (p[6]*plg[3,3] + p[42]*plg[3,5] + t81)*c2tloc +
                    (p[9]*plg[3,3] + p[43]*plg[3,5] + t82)*s2tloc)
    end

    # Terdiurnal
    # ==========

    if flags[:terdiurnal]
        t91 = (p[94]*plg[4,5] + p[47]*plg[4,7])*cd14*flags[:asym_annual]
        t92 = (p[95]*plg[4,5] + p[49]*plg[4,7])*cd14*flags[:asym_annual]

        t[14] = f2 * ( ( p[40]*plg[4,4] + t91 ) * s3tloc +
                       ( p[41]*plg[4,4] + t92 ) * c3tloc )
    end

    # Magnetic activity based on daily AP
    # ===================================

    if flags[:use_ap_array]
        ap = ap_array

        if p[52] != 0
            exp1 = exp( -10800abs(p[52])/( 1 + p[139]*( 45 - abs(g_lat) ) ) )

            (exp1  > 0.99999) && (exp1  = 0.99999)
            (p[25] < 1.0e-4)  && (p[25] = 1.0e-4)

            apt = _sg0(exp1,p,ap)

            t[9] = apt*( p[51] + p[97]*plg[1,3]+p[55]*plg[1,5] +
                        ( p[126]*plg[1,2] + p[127]*plg[1,4] + p[128]*plg[1,6] )*cd14*flags[:asym_annual] +
                        ( p[129]*plg[2,2] + p[130]*plg[2,4] + p[131]*plg[2,6] )*flags[:diurnal]*cos(hr*(tloc-p[132])))
        end
    else
        apd = ap - 4
        p44 = p[44]
        p45 = p[45]

        (p44 < 0) && (p44 = 1e-5)

        apdf = apd + ( p45 - 1 )*( apd + ( exp(-p44 * apd) - 1 )/p44 )

        if flags[:daily_ap]
            t[9] = apdf*( p[33] + p[46]*plg[1,3] + p[35]*plg[1,5] +
                         ( p[101]*plg[1,2] + p[102]*plg[1,4] + p[103]*plg[1,6])*cd14*flags[:asym_annual] +
                         ( p[122]*plg[2,2] + p[123]*plg[2,4] + p[124]*plg[2,6])*flags[:diurnal]*cos(hr*(tloc-p[125])))
        end
    end

    if flags[:all_ut_long_effects] && (g_long > - 1000)
        # Longitudinal
        # ============

        if flags[:longitudinal]
            t[11] = (1 + p[81]*dfa*flags[:F107_Mean]) * (
                     (  p[65]*plg[2,3] +  p[66]*plg[2,5] + p[67]*plg[2,7]  +
                       p[104]*plg[2,2] + p[105]*plg[2,4] + p[106]*plg[2,6] +
                       flags[:asym_annual]*( p[110]*plg[2,2] +
                                             p[111]*plg[2,4] +
                                             p[112]*plg[2,6])*cd14)*cos(dgtr*g_long) +
                     (  p[91]*plg[2,3] +  p[92]*plg[2,5] +  p[93]*plg[2,7] +
                       p[107]*plg[2,2] + p[108]*plg[2,4] + p[109]*plg[2,6] +
                       flags[:asym_annual]*( p[113]*plg[2,2] +
                                             p[114]*plg[2,4] +
                                             p[115]*plg[2,6])*cd14)*sin(dgtr*g_long)
                    )
        end

        # UT and Mixed UT, Longitude
        # ==========================

        if flags[:ut_mixed_ut_long]
            t[12]  = ( 1 +  p[96]*plg[1,2] )*(1 + p[82]*dfa*flags[:F107_Mean] )*
                     ( 1 + p[120]*plg[1,2]*flags[:asym_annual]*cd14)*
                     ( ( p[69]*plg[1,2] + p[70]*plg[1,4] + p[71]*plg[1,6])*cos(sr*(sec-p[72])))

            t[12] += flags[:longitudinal]*
                     ( p[77]*plg[3,4] + p[78]*plg[3,6] + p[79]*plg[3,8])*
                     cos(sr*(sec-p[80]) + 2*dgtr*g_long)*
                     (1 +p[138]*dfa*flags[:F107_Mean])
                 end

        # UT, Longitude Magnetic Activity
        # ===============================

        if flags[:mixed_ap_ut_long]
            if flags[:use_ap_array]
                if p[52] != 0
                    t[13]=apt*flags[:longitudinal]*( 1 + p[133]*plg[1,2] )*
                        (  p[53]*plg[2,3] +  p[99]*plg[2,5] +  p[68]*plg[2,7] )*cos(dgtr*(g_long-p[98])) +
                        apt*flags[:longitudinal]*flags[:asym_annual]*
                        ( p[134]*plg[2,2] + p[135]*plg[2,4] + p[136]*plg[2,6])*cd14*cos(dgtr*(g_long-p[137])) +
                        apt*flags[:ut_mixed_ut_long]*
                        (  p[56]*plg[1,2] + p[57]*plg[1,4]  +  p[58]*plg[1,6])*cos(sr*(sec-p[59]))
                end
            else
                t[13] = apdf*flags[:longitudinal]*(1 + p[121]*plg[1,2] )*
                    (  p[61]*plg[2,3] +  p[62]*plg[2,5] +  p[63]*plg[2,7])*cos(dgtr*(g_long-p[64])) +
                    apdf*flags[:longitudinal]*flags[:asym_annual]*
                    ( p[116]*plg[2,2] + p[117]*plg[2,4] + p[118]*plg[2,6])*cd14*cos(dgtr*(g_long-p[119])) +
                    apdf*flags[:ut_mixed_ut_long]*
                    (  p[84]*plg[1,2] +  p[85]*plg[1,4] +  p[86]*plg[1,6])*cos(sr*(sec-p[76]))
            end
        end
    end

    # Update the NRLMSISE-00 structure.
    @pack nrlmsise00d = apt, apdf

    # Parameters not used: 82, 89, 99, 139-149.
    tinf = p[31] +
           flags[:F107_Mean]*t[1] +
           flags[:time_independent]*t[2] +
           flags[:sym_annual]*t[3] +
           flags[:sym_semiannual]*t[4] +
           flags[:asym_annual]*t[5] +
           flags[:asym_semiannual]*t[6] +
           flags[:diurnal]*t[7] +
           flags[:semidiurnal]*t[8] +
           flags[:daily_ap]*t[9] +
           flags[:all_ut_long_effects]*t[10] +
           flags[:longitudinal]*t[11] +
           flags[:ut_mixed_ut_long]*t[12] +
           flags[:mixed_ap_ut_long]*t[13] +
           flags[:terdiurnal]*t[14]

    tinf
end

"""
    function _glob7s(p::AbstractVector{T}, nrlmsise00d::NRLMSISE00_Structure{T}) where T<:Number

Version of Globe for lower atmosphere (1999-10-26).

# Args

* `p`: Vector with the coefficients.
* `nrlmsise00d`: NRLMSISE-00 structure (see `NRLMSISE00_Structure`).

# Returns

The temperature (?).

"""
function _glob7s(p::AbstractVector{T},
                 nrlmsise00d::NRLMSISE00_Structure{T}) where T<:Number

    @unpack_NRLMSISE00_Structure nrlmsise00d

    # Constants
    # =========

    dgtr = T(1.74533e-2)
    dr   = T(1.72142e-2)
    hr   = T(0.2618)

    # Confirm parameter set.
    (p[100] == 0) && (p[100] = T(2))
    (p[100] != 2) && error("Wrong parameter set for glob7s.")

    # Initialization of variables.
    t = zeros(MVector{14,T})

    cd32 = cos( 1dr*(doy-p[32]) )
    cd18 = cos( 2dr*(doy-p[18]) )
    cd14 = cos( 1dr*(doy-p[14]) )
    cd39 = cos( 2dr*(doy-p[39]) )

    # F10.7
    # =====

    t[1] = p[22]*dfa

    # Time independent
    # ================

    t[2] =  p[2]*plg[1,3] +  p[3]*plg[1,5] + p[23]*plg[1,7] + p[27]*plg[1,2] +
           p[15]*plg[1,4] + p[60]*plg[1,6]

    # Symmetrical Annual
    # ==================

    t[3] = ( p[19] + p[48]*plg[1,3] + p[30]*plg[1,5] )*cd32

    # Symmetrical Semiannual
    # ======================

    t[4] = ( p[16] + p[17]*plg[1,3] + p[31]*plg[1,5] )*cd18

    # Asymmetrical Annual
    # ===================

    t[5] = ( p[10]*plg[1,2] + p[11]*plg[1,4] + p[21]*plg[1,6] )*cd14

    # Asymmetrical Semiannual
    # =======================

    t[6] = p[38]*plg[1,2]*cd39

    # Diurnal
    # =======

    if flags[:diurnal]
        t71  = p[12]*plg[2,3]*cd14*flags[:asym_annual]
        t72  = p[13]*plg[2,3]*cd14*flags[:asym_annual]
        t[7] = ( p[4]*plg[2,2] + p[5]*plg[2,4] + t71 ) * ctloc +
               ( p[7]*plg[2,2] + p[8]*plg[2,4] + t72 ) * stloc
    end

    # Semidiurnal
    # ===========

    if flags[:semidiurnal]
        t81  = (p[24]*plg[3,4]+p[36]*plg[3,6])*cd14*flags[:asym_annual]
        t82  = (p[34]*plg[3,4]+p[37]*plg[3,6])*cd14*flags[:asym_annual]
        t[8] = (p[6]*plg[3,3] + p[42]*plg[3,5] + t81) * c2tloc +
               (p[9]*plg[3,3] + p[43]*plg[3,5] + t82) * s2tloc
    end

    # Terdiurnal
    # ==========

    if flags[:terdiurnal]
        t[14] = p[40] * plg[4,4] * s3tloc + p[41] * plg[4,4] * c3tloc
    end

    # Magnetic Activity
    # =================
    if flags[:use_ap_array]
        t[9] = p[51]*apt + p[97]*plg[1,3] * apt * flags[:time_independent]
    elseif flags[:daily_ap]
        t[9] = apdf * (p[33] + p[46] * plg[1,3] * flags[:time_independent])
    end

    # Longitudinal
    # ============

    if !( !flags[:all_ut_long_effects] ||
          !flags[:longitudinal] ||
          (g_long<=-1000.0) )

        t[11] = ( 1 + plg[1,2]*( p[81]*flags[    :asym_annual]*cos( 1dr*(doy-p[82]) ) +
                                 p[86]*flags[:asym_semiannual]*cos( 2dr*(doy-p[87]) ) )+
                     p[84]*flags[    :sym_annual]*cos( 1dr*(doy-p[85]) ) +
                     p[88]*flags[:sym_semiannual]*cos( 2dr*(doy-p[89]) ) )*
                ( ( p[65]*plg[2,3] + p[66]*plg[2,5] +
                    p[67]*plg[2,7] + p[75]*plg[2,2] +
                    p[76]*plg[2,4] + p[77]*plg[2,6])*cos(dgtr*g_long) +
                  ( p[91]*plg[2,3] + p[92]*plg[2,5] +
                    p[93]*plg[2,7] + p[78]*plg[2,2] +
                    p[79]*plg[2,4] + p[80]*plg[2,6])*sin(dgtr*g_long))
    end

    tinf = flags[:F107_Mean]*t[1] +
           flags[:time_independent]*t[2] +
           flags[:sym_annual]*t[3] +
           flags[:sym_semiannual]*t[4] +
           flags[:asym_annual]*t[5] +
           flags[:asym_semiannual]*t[6] +
           flags[:diurnal]*t[7] +
           flags[:semidiurnal]*t[8] +
           flags[:daily_ap]*t[9] +
           flags[:all_ut_long_effects]*t[10] +
           flags[:longitudinal]*t[11] +
           flags[:ut_mixed_ut_long]*t[12] +
           flags[:mixed_ap_ut_long]*t[13] +
           flags[:terdiurnal]*t[14]

    tinf
end

"""
    function gts7(nrlmsise00d::NRLMSISE00_Structure{T}) where T<:Number

Thermospheric portion of NRLMSISE-00. This function should not be called to
compute NRLMSISE-00. Use `gtd7` or `gtd7d` instead.

# Args

* `nrlmsise00d`: An instance of `NRLMSISE00_Structure`.

# Returns

An instance of structure `NRLMSISE00_Structure` with the outputs.

"""
function gts7(nrlmsise00d::NRLMSISE00_Structure{T}) where T<:Number

    @unpack_NRLMSISE00_Structure nrlmsise00d

    # Constants
    # =========

    dgtr   = T(1.74533E-2)
    dr     = T(1.72142E-2)
    za     = pdl[2,16]
    alpha  = SVector{9,T}(-0.38, 0.0, 0.0, 0.0, 0.17, 0.0, -0.38, 0.0, 0.0)
    altl   = SVector{8,T}(200.0, 300.0, 160.0, 250.0, 240.0, 450.0, 320.0, 450.0)
    zn1    = SVector{5,T}(za, 110.0, 100.0, 90.0, 72.5)

    # Initialization of Variables
    # ===========================

    T_alt     = T(0)
    meso_tn1  = zeros(MVector{5,T})
    meso_tgn1 = zeros(MVector{2,T})

    # Tinf variations not important below `za` or `zn1[1]`
    # ====================================================

    tinf = (alt > zn1[1]) ? ptm[1] * pt[1] * ( 1 + flags[:all_tinf_var]*_globe7!(pt,nrlmsise00d) ) :
                            ptm[1] * pt[1]

    T_exo = tinf

    # Gradient variations not important below `zn1[5]`
    # ================================================

    g0 = (alt > zn1[5]) ? ptm[4] * ps[1] * ( 1 + flags[:all_s_var]*_globe7!(ps,nrlmsise00d) ) :
                          ptm[4] * ps[1]

    tlb = ptm[2] * (1 + flags[:all_tlb_var]*_globe7!(pd[4,:],nrlmsise00d))*pd[4,1]
    s = g0 / (tinf - tlb)

    # Lower thermosphere temperature variations not significant for density
    # above 300 km.

    if alt < 300
        meso_tn1[2]  = ptm[7]*ptl[1,1]/
                       (1-flags[:all_tn1_var]*_glob7s(ptl[1,:], nrlmsise00d))
        meso_tn1[3]  = ptm[3]*ptl[2,1]/
                       (1-flags[:all_tn1_var]*_glob7s(ptl[2,:], nrlmsise00d))
        meso_tn1[4]  = ptm[8]*ptl[3,1]/
                       (1-flags[:all_tn1_var]*_glob7s(ptl[3,:], nrlmsise00d))
        meso_tn1[5]  = ptm[5]*ptl[4,1]/
                       (1-flags[:all_tn1_var]*flags[:all_tn2_var]*_glob7s(ptl[4,:], nrlmsise00d))
        meso_tgn1[2] = ptm[9]*pma[9,1]*
                       (1 + flags[:all_tn1_var]*flags[:all_tn2_var]*_glob7s(pma[9,:], nrlmsise00d))*
                       meso_tn1[5]^2/(ptm[5]*ptl[4,1])^2
    else
        meso_tn1[2]  = ptm[7]*ptl[1,1]
        meso_tn1[3]  = ptm[3]*ptl[2,1]
        meso_tn1[4]  = ptm[8]*ptl[3,1]
        meso_tn1[5]  = ptm[5]*ptl[4,1]
        meso_tgn1[2] = ptm[9]*pma[9,1]*meso_tn1[5]^2/(ptm[5]*ptl[4,1])^2
    end

    # N2 variation factor at Zlb.
    g28 = flags[:all_nlb_var]*_globe7!(pd[3,:], nrlmsise00d)

    # Variation of Turbopause Height
    # ==============================

    zhf = pdl[2,25]*
          (1 + flags[:asym_annual]*pdl[1,25]*sin(dgtr*g_lat)*cos(dr*(doy-pt[14])))
    xmm = pdm[3,5]

    # N2 Density
    # ==========

    # Diffusive density at Zlb.
    db28 = pdm[3,1]*exp(g28)*pd[3,1]

    # Diffusive density at Alt.
    den_N2, T_alt = _densu(re,
                           gsurf,
                           alt,
                           db28,
                           tinf,
                           tlb,
                           T(28),
                           alpha[3],
                           ptm[6],
                           s,
                           zn1,
                           meso_tn1,
                           meso_tgn1)

    # Turbopause.
    zh28  = pdm[3,3]*zhf
    zhm28 = pdm[3,4]*pdl[2,6]
    xmd   = 28 - xmm

    # Mixed density at Zlb.
    b28, tz = _densu(re,
                     gsurf,
                     zh28,
                     db28,
                     tinf,
                     tlb,
                     xmd,
                     alpha[3]-1,
                     ptm[6],
                     s,
                     zn1,
                     meso_tn1,
                     meso_tgn1)

    if flags[:departures_from_eq] && (alt < altl[3])
        # Mixed density at Alt.
        dm28, tz = _densu(re,
                          gsurf,
                          alt,
                          b28,
                          tinf,
                          tlb,
                          xmm,
                          alpha[3],
                          ptm[6],
                          s,
                          zn1,
                          meso_tn1,
                          meso_tgn1)

        # Net density at Alt.
        den_N2 = _dnet(den_N2, dm28, zhm28, xmm, T(28))
    end

    # He Density
    # ==========

    # Density variation factor at Zlb.
    g4   = flags[:all_nlb_var]*_globe7!(pd[1,:], nrlmsise00d)

    # Diffusive density at Zlb.
    db04 = pdm[1,1]*exp(g4)*pd[1,1]

    # Diffusive density at Alt.
    den_He, T_alt = _densu(re,
                           gsurf,
                           alt,
                           db04,
                           tinf,
                           tlb,
                           T(4),
                           alpha[1],
                           ptm[6],
                           s,
                           zn1,
                           meso_tn1,
                           meso_tgn1)

    if flags[:departures_from_eq] && (alt < altl[3])
        # Turbopause.
        zh04 = pdm[1,3]

        # Mixed density at Zlb.
        b04, T_alt = _densu(re,
                            gsurf,
                            zh04,
                            db04,
                            tinf,
                            tlb,
                            4-xmm,
                            alpha[1]-1,
                            ptm[6],
                            s,
                            zn1,
                            meso_tn1,
                            meso_tgn1)

        # Mixed density at Alt.
        dm04, T_alt = _densu(re,
                             gsurf,
                             alt,
                             b04,
                             tinf,
                             tlb,
                             xmm,
                             T(0),
                             ptm[6],
                             s,
                             zn1,
                             meso_tn1,
                             meso_tgn1)

        zhm04 = zhm28

        # Net density at Alt.
        den_He = _dnet(den_He, dm04, zhm04, xmm, T(4))

        # Correction to specified mixing ration at ground.
        rl   = log(b28*pdm[1,2]/b04)
        zc04 = pdm[1,5]*pdl[2,1]
        hc04 = pdm[1,6]*pdl[2,2]

        # Net density corrected at Alt.
        den_He *= _ccor(alt, rl, hc04, zc04)
    end

    # O Density
    # =========

    # Density variation factor at Zlb.
    g16 = flags[:all_nlb_var]*_globe7!(pd[2,:], nrlmsise00d)

    #  Diffusive density at Zlb.
    db16 = pdm[2,1]*exp(g16)*pd[2,1]

    # Diffusive density at Alt.
    den_O, T_alt = _densu(re,
                          gsurf,
                          alt,
                          db16,
                          tinf,
                          tlb,
                          T(16),
                          alpha[2],
                          ptm[6],
                          s,
                          zn1,
                          meso_tn1,
                          meso_tgn1)

    if flags[:departures_from_eq] && (alt <= altl[2])
        # Turbopause.
        zh16 = pdm[2,3]

        # Mixed density at Zlb.
        b16, T_alt = _densu(re,
                            gsurf,
                            zh16,
                            db16,
                            tinf,
                            tlb,
                            16-xmm,
                            alpha[2]-1,
                            ptm[6],
                            s,
                            zn1,
                            meso_tn1,
                            meso_tgn1)

        # Mixed density at Alt.
        dm16, T_alt = _densu(re,
                             gsurf,
                             alt,
                             b16,
                             tinf,
                             tlb,
                             xmm,
                             T(0),
                             ptm[6],
                             s,
                             zn1,
                             meso_tn1,
                             meso_tgn1)

        zhm16 = zhm28;

        # Net density at Alt.
        den_O  = _dnet(den_O, dm16, zhm16, xmm, T(16))
        rl     = pdm[2,2]*pdl[2,17]*(1+flags[:F107_Mean]*pdl[1,24]*dfa)
        hc16   = pdm[2,6]*pdl[2,4]
        zc16   = pdm[2,5]*pdl[2,3]
        hc216  = pdm[2,6]*pdl[2,5]
        den_O *= _ccor2(alt, rl, hc16, zc16, hc216)

        # Chemistry correction.
        hcc16 = pdm[2,8]*pdl[2,14]
        zcc16 = pdm[2,7]*pdl[2,13]
        rc16  = pdm[2,4]*pdl[2,15]

        # Net density corrected at Alt.
        den_O  *= _ccor(alt, rc16, hcc16, zcc16)
    end

    # O2 Density
    # ==========

    # Density variation factor at Zlb.
    g32 = flags[:all_nlb_var]*_globe7!(pd[5,:], nrlmsise00d)

    # Diffusive density at Zlb.
    db32 = pdm[4,1]*exp(g32)*pd[5,1];

    # Diffusive density at Alt.
    den_O2, T_alt = _densu(re,
                           gsurf,
                           alt,
                           db32,
                           tinf,
                           tlb,
                           T(32),
                           alpha[4],
                           ptm[6],
                           s,
                           zn1,
                           meso_tn1,
                           meso_tgn1)

    if flags[:departures_from_eq]
        if alt <= altl[4]
            # Turbopause.
            zh32=pdm[4,3]

            # Mixed density at Zlb.
            b32, T_alt = _densu(re,
                                gsurf,
                                zh32,
                                db32,
                                tinf,
                                tlb,
                                32-xmm,
                                alpha[4]-1,
                                ptm[6],
                                s,
                                zn1,
                                meso_tn1,
                                meso_tgn1)

            # Mixed density at Alt.
            dm32, T_alt = _densu(re,
                                 gsurf,
                                 alt,
                                 b32,
                                 tinf,
                                 tlb,
                                 xmm,
                                 T(0),
                                 ptm[6],
                                 s,
                                 zn1,
                                 meso_tn1,
                                 meso_tgn1)

            zhm32 = zhm28

            # Net density at Alt.
            den_O2 = _dnet(den_O2, dm32, zhm32, xmm, T(32))

            # Correction to specified mixing ratio at ground.
            rl      = log(b28*pdm[4,2]/b32)
            hc32    = pdm[4,6]*pdl[2,8]
            zc32    = pdm[4,5]*pdl[2,7]
            den_O2 *= _ccor(alt, rl, hc32, zc32)
        end

        # Correction for general departure from diffusive equilibrium above Zlb.
        hcc32  = pdm[4,8]*pdl[2,23]
        hcc232 = pdm[4,8]*pdl[1,23]
        zcc32  = pdm[4,7]*pdl[2,22]
        rc32   = pdm[4,4]*pdl[2,24]*(1 + flags[:F107_Mean]*pdl[1,24]*dfa)

        # Net density corrected at Alt.
        den_O2 *= _ccor2(alt, rc32, hcc32, zcc32, hcc232)
    end

    # Ar Density
    # ==========

    # Density variation factor at Zlb.
    g40 = flags[:all_nlb_var]*_globe7!(pd[6,:], nrlmsise00d)

    # Diffusive density at Zlb.
    db40 = pdm[5,1]*exp(g40)*pd[6,1];

    # Diffusive density at Alt.
    den_Ar, T_alt = _densu(re,
                           gsurf,
                           alt,
                           db40,
                           tinf,
                           tlb,
                           T(40),
                           alpha[5],
                           ptm[6],
                           s,
                           zn1,
                           meso_tn1,
                           meso_tgn1)

    if flags[:departures_from_eq] && ( alt <= altl[5] )
        # Turbopause.
        zh40=pdm[5,3]

        # Mixed density at Zlb.
        b40, T_alt = _densu(re,
                            gsurf,
                            zh40,
                            db40,
                            tinf,
                            tlb,
                            40-xmm,
                            alpha[5]-1,
                            ptm[6],
                            s,
                            zn1,
                            meso_tn1,
                            meso_tgn1)

        # Mixed density at Alt.
        dm40, T_alt = _densu(re,
                             gsurf,
                             alt,
                             b40,
                             tinf,
                             tlb,
                             xmm,
                             T(0),
                             ptm[6],
                             s,
                             zn1,
                             meso_tn1,
                             meso_tgn1)

        zhm40= zhm28

        # Net density at Alt.
        den_Ar = _dnet(den_Ar, dm40, zhm40, xmm, T(40))

        # Correction to specified mixing ratio at ground.
        rl   = log(b28*pdm[5,2]/b40)
        hc40 = pdm[5,6]*pdl[2,10]
        zc40 = pdm[5,5]*pdl[2,9]

        # Net density corrected at Alt.
        den_Ar = den_Ar*_ccor(alt, rl, hc40, zc40)
    end

    # H Density
    # =========

    # Density variation factor at Zlb.
    g1 = flags[:all_nlb_var]*_globe7!(pd[7,:], nrlmsise00d)

    # Diffusive density at Zlb.
    db01 = pdm[6,1]*exp(g1)*pd[7,1];

    # Diffusive density at Alt.
    den_H, T_alt = _densu(re,
                          gsurf,
                          alt,
                          db01,
                          tinf,
                          tlb,
                          T(1),
                          alpha[7],
                          ptm[6],
                          s,
                          zn1,
                          meso_tn1,
                          meso_tgn1)

    if flags[:departures_from_eq] && ( alt <= altl[7] )
        # Turbopause.
        zh01=pdm[6,3]

        # Mixed density at Zlb.
        b01, T_alt = _densu(re,
                            gsurf,
                            zh01,
                            db01,
                            tinf,
                            tlb,
                            1-xmm,
                            alpha[7]-1,
                            ptm[6],
                            s,
                            zn1,
                            meso_tn1,
                            meso_tgn1)

        # Mixed density at Alt.
        dm01, T_alt = _densu(re,
                             gsurf,
                             alt,
                             b01,
                             tinf,
                             tlb,
                             xmm,
                             T(0),
                             ptm[6],
                             s,
                             zn1,
                             meso_tn1,
                             meso_tgn1)

        zhm01 = zhm28

        # Net density at Alt.
        den_H = _dnet(den_H, dm01, zhm01, xmm, T(1))

        # Correction to specified mixing ratio at ground.
        rl     = log(b28*pdm[6,2]*abs(pdl[2,18])/b01)
        hc01   = pdm[6,6]*pdl[2,12]
        zc01   = pdm[6,5]*pdl[2,11]
        den_H *= _ccor(alt, rl, hc01, zc01)

        # Chemistry correction.
        hcc01 = pdm[6,8]*pdl[2,20]
        zcc01 = pdm[6,7]*pdl[2,19]
        rc01  = pdm[6,4]*pdl[2,21]

        # Net density corrected at Alt.
        den_H *= _ccor(alt, rc01, hcc01, zcc01)
    end

    # N Density
    # =========

    # Density variation factor at Zlb.
    g14 = flags[:all_nlb_var]*_globe7!(pd[8,:],nrlmsise00d)

    # Diffusive density at Zlb.
    db14 = pdm[7,1]*exp(g14)*pd[8,1]

    # Diffusive density at Alt.
    den_N, T_alt = _densu(re,
                         gsurf,
                         alt,
                         db14,
                         tinf,
                         tlb,
                         T(14),
                         alpha[8],
                         ptm[6],
                         s,
                         zn1,
                         meso_tn1,
                         meso_tgn1)

    if flags[:departures_from_eq] && ( alt <= altl[8] )
        # Turbopause.
        zh14 = pdm[7,3]

        # Mixed density at Zlb.
        b14, T_alt = _densu(re,
                            gsurf,
                            zh14,
                            db14,
                            tinf,
                            tlb,
                            14-xmm,
                            alpha[8]-1,
                            ptm[6],
                            s,
                            zn1,
                            meso_tn1,
                            meso_tgn1)

        #  Mixed density at Alt.
        dm14, T_alt = _densu(re,
                             gsurf,
                             alt,
                             b14,
                             tinf,
                             tlb,
                             xmm,
                             T(0),
                             ptm[6],
                             s,
                             zn1,
                             meso_tn1,
                             meso_tgn1)

        zhm14 = zhm28

        # Net density at Alt.
        den_N = _dnet(den_N, dm14, zhm14, xmm, T(14))

        # Correction to specified mixing ratio at ground.
        rl     = log(b28*pdm[7,2]*abs(pdl[1,3])/b14)
        hc14   = pdm[7,6]*pdl[1,2]
        zc14   = pdm[7,5]*pdl[1,1]
        den_N *= _ccor(alt, rl, hc14, zc14)

        # Chemistry correction.
        hcc14 = pdm[7,8]*pdl[1,5]
        zcc14 = pdm[7,7]*pdl[1,4]
        rc14  = pdm[7,4]*pdl[1,6]

        # Net density corrected at Alt.
        den_N *= _ccor(alt, rc14, hcc14, zcc14)
    end

    # Anomalous O Density
    # ===================

    g16h          = flags[:all_nlb_var]*_globe7!(pd[9,:], nrlmsise00d)
    db16h         = pdm[8,1]*exp(g16h)*pd[9,1]
    tho           = pdm[8,10]*pdl[1,7]
    den_aO, T_alt = _densu(re,
                           gsurf,
                           alt,
                           db16h,
                           tho,
                           tho,
                           T(16),
                           alpha[9],
                           ptm[6],
                           s,
                           zn1,
                           meso_tn1,
                           meso_tgn1)

    zsht = pdm[8,6]
    zmho = pdm[8,5]
    zsho = _scalh(zmho, T(16), tho, gsurf, re)

    den_aO *= exp( -zsht/zsho*( exp( -( alt - zmho )/zsht ) - 1 ) )

    # Total Mass Density
    # ==================

    den_Total = 1.66e-24( 4den_He +
                         16den_O  +
                         28den_N2 +
                         32den_O2 +
                         40den_Ar +
                           den_H  +
                         14den_N)

    # Temperature at Selected Altitude
    # ================================

    ddum, T_alt = _densu(re,
                         gsurf,
                         abs(alt),
                         T(1),
                         tinf,
                         tlb,
                         T(0),
                         T(0),
                         ptm[6],
                         s,
                         zn1,
                         meso_tn1,
                         meso_tgn1)

    # Output
    # ======

    # Check if we should change the unit.
    if flags[:output_m_kg]
        den_He    *= 1e6
        den_O     *= 1e6
        den_N2    *= 1e6
        den_O2    *= 1e6
        den_Ar    *= 1e6
        den_Total *= 1e6/1000
        den_H     *= 1e6
        den_N     *= 1e6
        den_aO    *= 1e6
    end

    # Repack variables that were modified.
    meso_tn1_5  = meso_tn1[5]
    meso_tgn1_2 = meso_tgn1[2]
    @pack nrlmsise00d = meso_tn1_5, meso_tgn1_2, dm28

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

@inline function _scalh(alt::T, xm::T, temp::T, gsurf::T, re::T) where T<:Number
    rgas = T(831.4)
    g    = gsurf / (1 + alt/re)^2
    rgas * temp / (g * xm)
end

"""
    @inline function _splini(xa::StaticVector{N,T}, ya::StaticVector{N,T}, y2a::StaticVector{N,T}, x::T) where T<:Number where N

Compute the integral of the cubic spline function from `xa[1]` to `x`.

# Args

* `xa`: X components of the tabulated function in ascending order.
* `ya`: Y components of the tabulated function evaluated at `xa`.
* `y2a`: Second derivatives.
* `x`: Abscissa endpoint for integration.

# Returns

The integral of cubic spline function from `xa[1]` to `x`.

"""
@inline function _splini(xa::StaticVector{N,T},
                         ya::StaticVector{N,T},
                         y2a::StaticVector{N,T},
                         x::T) where T<:Number where N
    yi  = T(0)
    klo = 1
    khi = 2

    while ( x > xa[klo] ) && ( khi <= N )
        xx = x

        if khi <= (N-1)
            xx = (x < xa[khi]) ? x : xa[khi]
        end

        h = xa[khi] - xa[klo]
        a = (xa[khi] - xx)/h
        b = (xx - xa[klo])/h

        a² = a^2
        b² = b^2
        a⁴ = a^4
        b⁴ = b^4

        coef_a = -(1 + a⁴)/4 + a²/2
        coef_b = b⁴/4 - b²/2

        yi += h * ( (1 - a²) *  ya[klo] / 2 +     b² *  ya[khi] / 2 +
                   ( coef_a  * y2a[klo]     + coef_b * y2a[khi] ) * h^2 / 6)

        klo += 1
        khi += 1
    end

    yi
end


"""
    @inline function _spline(x::StaticVector{N,T}, y::StaticVector{N,T}, yp1::T, ypn::T) where T<:Number where N

Compute the 2nd derivatives of cubic spline interpolation function tabulated by
`x` and `y` given the 2nd derivatives values at `x[1]` (`yp1`) and at `x[N]`
(`ypn`).

This function was adapted from Numerical Recipes.

# Args

* `x`: X components of the tabulated function in ascending order.
* `y`: Y components of the tabulated function evaluated at `x`.
* `yp1`: 2nd derivative value at `x[1]`.
* `ypn`: 2nd derivative value at `x[N]`.

# Returns

The 2nd derivative of cubic spline interpolation function evaluated at `x`.

# Remarks

Values higher than `1e30` in the 2nd derivatives at the borders (`yp1` and
`ypn`) are interpreted as `0`.

"""
@inline function _spline(x::StaticVector{N,T},
                         y::StaticVector{N,T},
                         yp1::T,
                         ypn::T) where T<:Number where N
    u  = zeros(MVector{N,T})
    y2 = zeros(MVector{N,T})

    if (yp1 > 0.99e30)
        y2[1] = 0
        u[1]  = 0
    else
        y2[1] = T(-0.5)
        u[1]  = ( 3/(x[2]-x[1]) )*( (y[2]-y[1])/(x[2]-x[1]) - yp1)
    end

    for i=2:N-1
        sig   = (x[i]-x[i-1])/(x[i+1] - x[i-1])
        p     = sig * y2[i-1] + 2
        y2[i] = (sig - 1) / p
        m_a   = (y[i+1] - y[i]  ) / (x[i+1] - x[i]  )
        m_b   = (y[i]   - y[i-1]) / (x[i]   - x[i-1])
        u[i]  = (6(m_a - m_b)/(x[i+1] - x[i-1]) - sig * u[i-1])/p
    end

    if ypn > 0.99e30
        qn = T(0)
        un = T(0)
    else
        qn = T(0.5)
        un = ( 3/(x[N] - x[N-1]) )*( ypn - (y[N] - y[N-1])/(x[N] - x[N-1]) )
    end

    y2[N] = (un - qn * u[N-1]) / (qn * y2[N-1] + 1)

    for k = N-1:-1:1
        y2[k] = y2[k] * y2[k+1] + u[k]
    end

    y2
end

"""
    @inline function _splint(xa::StaticVector{N,T}, ya::StaticVector{N,T}, y2a::StaticVector{N,T}, x::T) where T<:Number where N

Compute the cubic spline interpolation value at `x`.

This function was adapted from Numerical Recipes.

# Args

* `xa`: X components of the tabulated function in ascending order.
* `ya`: Y components of the tabulated function evaluated at `xa`.
* `y2a`: Second derivatives.
* `x`: Abscissa endpoint for interpolation.

# Returns

The cubic spline interpolation value at `x`.

"""
@inline function _splint(xa::StaticVector{N,T},
                         ya::StaticVector{N,T},
                         y2a::StaticVector{N,T},
                         x::T) where T<:Number where N
    klo = 1
    khi = N

    while (khi - klo) > 1
        k = round(Int,(khi + klo)/2)

        if xa[k] > x
            khi = k
        else
            klo = k
        end
    end

    h = xa[khi] - xa[klo]

    (h == 0) && error("Bad xa input to splint.")

    a  = (xa[khi] -    x   )/h
    b  = (   x    - xa[klo])/h
    yi = a * ya[klo] + b * ya[khi] +
         ( (a^3 - a) * y2a[klo] + (b^3 - b) * y2a[khi]) * h^2/6

    yi
end

################################################################################
#                                  Equations
################################################################################

@inline _zeta(re::Number, zz::Number, zl::Number) = (zz-zl)*(re+zl)/(re+zz)

# 3hr Magnetic Activity Functions
# ===============================

# Eq. A24d.
@inline _g0(a::Number, p::AbstractVector) =
    (a-4 + (p[26] - 1)*(a-4 + ( exp( -abs(p[25])*(a-4) ) -1 ) / abs(p[25]) ) )

# Eq. A24c.
@inline _sumex(ex::Number) = (1 + (1 - ex^19) / (1 - ex) * sqrt(ex) )

# Eq. A24a.
@inline _sg0(ex::Number, p::AbstractVector, ap::AbstractVector) =
    (_g0(ap[2],p) + (_g0(ap[3],p)*ex + _g0(ap[4],p)*ex^2 + _g0(ap[5],p)*ex^3 +
                    (_g0(ap[6],p)*ex^4 + _g0(ap[7],p)*ex^12.0)*
                    (1-ex^8)/(1-ex)))/_sumex(ex)
