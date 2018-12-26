#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Private functions and variables of the Jacchia-Bowman 2008 Atmospheric
#   Model, a product of the Space Environment Technologies.
#
#   The source code is available at the following git:
#
#       http://sol.spacenvironment.net/jb2008/code.html
#
#   For more information about the model, see:
#
#       http://sol.spacenvironment.net/~JB2006/
#       http://sol.spacenvironment.net/~JB2008/
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Bowman, B. R., Tobiska, W. K., Marcos, F. A., Huang, C. Y., Lin, C. S.,
#       Burke, W. J (2008). A new empirical thermospheric density model JB2008
#       using new solar and geomagnetic indices. AIAA/AAS Astrodynamics
#       Specialist Conference, Honolulu, Hawaii.
#
#   [2] Bowman, B. R., Tobiska, W. K., Marcos, F. A., Valladares, C (2007). The
#       JB2006 empirical thermospheric density model. Journal of Atmospheric and
#       Solar-Terrestrial Physics, v. 70, p. 774-793.

#   [3] Jacchia, L. G (1970). New static models of the thermosphere and
#       exosphere with empirical temperature profiles. SAO Special Report #313.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

################################################################################
#                              Private Constants
################################################################################

# Constants to compute the `ΔTc` correction.
const _jb2008_B =
    (-0.457512297e+01, -0.512114909e+01, -0.693003609e+02,  0.203716701e+03,
      0.703316291e+03, -0.194349234e+04,  0.110651308e+04, -0.174378996e+03,
      0.188594601e+04, -0.709371517e+04,  0.922454523e+04, -0.384508073e+04,
     -0.645841789e+01,  0.409703319e+02, -0.482006560e+03,  0.181870931e+04,
     -0.237389204e+04,  0.996703815e+03,  0.361416936e+02,)

const _jb2008_C =
    (-0.155986211e+02, -0.512114909e+01, -0.693003609e+02,  0.203716701e+03,
      0.703316291e+03, -0.194349234e+04,  0.110651308e+04, -0.220835117e+03,
      0.143256989e+04, -0.318481844e+04,  0.328981513e+04, -0.135332119e+04,
      0.199956489e+02, -0.127093998e+02,  0.212825156e+02, -0.275555432e+01,
      0.110234982e+02,  0.148881951e+03, -0.751640284e+03,  0.637876542e+03,
      0.127093998e+02, -0.212825156e+02,  0.275555432e+01)

# F(z) global model values, 1997 - 2006 fit.
const _jb2008_fzm = (0.2689e+00,-0.1176e-01, 0.2782e-01,-0.2782e-01,0.3470e-03)

# G(t) global model values, 1997 - 2006 fit.
const _jb2008_gtm = (-0.3633e+00, 0.8506e-01, 0.2401e+00,-0.1897e+00,
                     -0.2554e+00,-0.1790e-01, 0.5650e-03,-0.6407e-03,
                     -0.3418e-02,-0.1252e-02)

# Coefficients for high altitude density correction.
const _jb2008_cht = (0.22e0, -0.20e-02, 0.115e-02, -0.211e-05)

################################################################################
#                              Private Functions
################################################################################

"""
    function _jb2008_grav(z::R) where R

Compute the gravity [m/s] at altitude `z` [km] according to the model Jacchia
1971 [3].

"""
@inline function _jb2008_grav(z::R) where R
    # Mean Earth radius [km].
    Re = R(6356.776)

    # Gravity at Earth surface [m/s²].
    g₀ = R(9.80665)

    # Gravity at desired altitude [m/s²].
    g₀*(1 + z/Re)^(-2)
end

"""
    function _jb2008_highaltitude(h::Number, F10ₐ::Number)

Compute the high altitude exospheric density correction factor in altitude `h`
[km] and the averaged 10.7-cm solar flux (81-day centered on input time)
\\[10⁻²² W/(M² Hz)].

This function uses the model in Section 6.2 of [2].

"""
function _jb2008_highaltitude(h::Number, F10ₐ::Number)
    # Auxiliary variables.
    C = _jb2008_cht

    # Compute the high-altitude density correction.
    FρH = 1.0
    @inbounds if 1000 <= h <= 1500
        z = (h - 1000)/500

        # In this case, the `F₁₅₀₀` is the density factor at 1500 km.
        F₁₅₀₀ = C[1] + C[2]*F10ₐ + C[3]*1500 + C[4]*1500*F10ₐ

        ∂F₁₅₀₀∂z = 500( C[3] + C[4]*F10ₐ )

        # In [2], there is an error when computing this coefficients (eq. 11).
        # The partial derivative has the `500` factor and the coefficients have
        # other `500` factor that **is not** present in the source-code.
        c0  = +1
        c1  = +0
        c2  = +3F₁₅₀₀ - ∂F₁₅₀₀∂z - 3
        c3  = -2F₁₅₀₀ + ∂F₁₅₀₀∂z + 2
        FρH = @evalpoly(z, c0, c1, c2, c3)
    elseif h > 1500
        FρH = C[1] + C[2]*F10ₐ + C[3]*h + C[4]*h*F10ₐ
    end

    FρH
end

"""
    function _jb2008_int(z₀::Number, z₁::Number, R::Number, Tx::Number, T∞::Number, δf::Function)

Compute the integral of the function `δf` between `z₀` and `z₁` using the
Newton-Cotes 4th degree method. `R` is a number that defines the step size, `Tx`
is the temperature at the inflection point, and `T∞` is the exospheric
temperature.

The signature of the function `δf` is:

    δf(z, Tx, T∞)

and it must be `_jb2008_δf1` or `_jb2008_δf2`.

This function returns a tuple containing the integral and last value of `z` used
in the numerical algorithm.

"""
function _jb2008_int(z₀::Number, z₁::Number, R::Number, Tx::Number, T∞::Number,
                     δf::Function)
    # Weights for the Newton-Cotes Five-Point Quad. Formula.
    WT = (14/45, 64/45, 24/45, 64/45, 14/45)

    # Compute the number of integration steps.
    #
    # This is computed so that `z₂ = z₁*(zr)^n`. Hence, `zr` is the factor that
    # defines the size of each integration interval.
    al = log(z₁/z₀)
    n  = floor(Int64, al/R) + 1
    zr = exp(al/n)

    # Initialize the integration auxiliary variables.
    zi₁ = z₀
    zj  = 0.0
    Mbj = 0.0
    Tlj = 0.0
    int = 0.0 # This variable stores the integral from `z₀` to `z₁`.

    # For each integration step, use the Newton-Cotes 4th degree formula to
    # integrate (Boole's rule).
    @inbounds for i = 1:n
        zi₀   = zi₁                   # The beginning of the i-th integration step.
        zi₁   = zr*zi₁                # The end of the i-th integration step.
        Δz    = (zi₁ - zi₀)/4         # Step for the i-th integration step.
        int_i = WT[1]*δf(zi₀, Tx, T∞) # First term of the Newton-Cotes 4th degree sum.

        # Compute the Newton-Cotes 4th degree sum.
        for j = 2:5
            zj  = zi₀ + (j-1)*Δz

            # Value of the integrand at `zj`.
            int_i += WT[j] * δf(zj, Tx, T∞)
        end
        int += int_i * Δz
    end

    int, zj
end

"""
    function _jb2008_M(z::R) where R

Compute the mean molecular mass at altitude `z` [km] using the empirical profile
in eq. 1 [3].

"""
@inline function _jb2008_M(z::R) where R
    !(90 <= z < 105.1) && @warn "The empirical model for the mean molecular mass is valid only for 90 <= z <= 105 km."

    @evalpoly(z - 100, +28.15204,
                        -0.085586,
                        +1.2840e-4,
                        -1.0056e-5,
                        -1.0210e-5,
                        +1.5044e-6,
                        +9.9826e-8)
end

"""
    function _jb2008_semiannual(doy::Number, h::Number, F10ₐ::Number, S10ₐ::Number, M10ₐ::Number)

Compute the semiannual variation of the density considering the JB2008 model
[1].

# Args

* `doy`: Day of the year + fraction of the day.
* `h`: Height [km].
* `F10ₐ`: Averaged 10.7-cm flux (81-day centered on input-time)
          \\[10⁻²² W/(M² Hz)].
* `S10ₐ`: EUV 81-day averaged centered index.
* `M10ₐ`: MG2 81-day averaged centered index.

# Returns

* Semiannual F(z) heigh function.
* Semiannual G(t) yearly periodic function.
* Semiannual variation of the density `Δsalog₁₀ρ`.

"""
function _jb2008_semiannual(doy::Number, h::Number, F10ₐ::Number, S10ₐ::Number,
                            M10ₐ::Number)
    # Auxiliary variables.
    B = _jb2008_fzm
    C = _jb2008_gtm
    z = h/1000
    ω = 2π*(doy-1)/365  # See eq. 4 [2].


    # Compute the new 81-day centered solar index for F(z) according to eq. 4
    # [1].
    Fsmjₐ = 1F10ₐ - 0.7S10ₐ - 0.04M10ₐ

    # Compute the semiannual F(z) height function according to eq. 5 [1].
    @inbounds Fz = B[1] + B[2]*Fsmjₐ + B[3]*z*Fsmjₐ + B[4]*z^2*Fsmjₐ + B[5]*z*Fsmjₐ^2

    # Compute the new 81-day centered solar index for G(t) according to eq. 6
    # [1].
    Fsmₐ = 1F10ₐ - 0.75S10ₐ - 0.37M10ₐ

    # Compute the semiannual G(t) yearly periodic function according to eq. 7
    # [1].
    sω,  cω  = sincos(1ω)
    s2ω, c2ω = sincos(2ω)

    @inbounds Gt =        C[1] + C[2]*sω + C[3]*cω + C[4]*s2ω + C[5]*c2ω +
                   Fsmₐ*( C[6] + C[7]*sω + C[8]*cω + C[9]*s2ω + C[10]*c2ω )

    (Fz < 1e-6) && (Fz = 1e-6)

    Δsalog₁₀ρ = Fz*Gt

    (Fz, Gt, Δsalog₁₀ρ)
end

"""
    function _jb2008_T(z::R, Tx::R, T∞::R) where R<:Number

Compute the temperature [K] at height `z` [km] given the temperature `Tx` [K] at
the inflection point, and the exospheric temperature `T∞` [K] according to the
theory of the model Jacchia 1971 [3].

The inflection point is considered to by `z = 125 km`.

"""
function _jb2008_T(z::R, Tx::R, T∞::R) where R<:Number
    # Constants
    # =========

    T₁  = 183    # Temperature at the lower bound [K].
    z₁  = 90     # Altitude of the lower bound [km].
    zx  = 125    # Altitude of the inflection point [km].
    Δz₁ = z₁-zx

    # Check the parameters
    # ====================
    (z  < z₁) && error("The altitude must not be lower than $(z₁) km.")
    (T∞ < 0 ) && error("The exospheric temperature must be positive.")

    # Compute the temperature gradient at the inflection point.
    Gx = R(1.9)*( Tx - T₁ )/( zx - z₁ )

    # Compute the temperature at desire altitude
    # ==========================================

    Δz = z - zx

    if z <= zx
        c₁ = Gx
        c₂ = R(0)
        c₃ = R(-5.1)/(R(5.7)*Δz₁^2)*Gx
        c₄ = R( 0.8)/(R(1.9)*Δz₁^3)*Gx
        T  = @evalpoly(Δz, Tx, c₁, c₂, c₃, c₄)
    else
        A  = R(2/π) * ( T∞ - Tx )
        T  = Tx + A*atan( Gx/A * Δz * ( 1 + R(4.5e-6) * Δz^2.5 ) )
    end

    T
end

"""
    function _jb2008_ΔTc(F10::Number, lst::Number, glat::Number, h::Number)

Compute the correction in the `Tc` for Jacchia-Bowman model.

This correction is mention in [2]. However, the equations do not seem to match
those in the source-code. The ones implemented here are exactly the same as in
the source-code.

# Args

* `F10`: F10.7 flux.
* `lst`: Local solar time (0 - 24) \\[hr].
* `glat`: Geocentric latitude [rad].
* `h`: Altitude [km].

# Returns

The correction `ΔTc` [K].

"""
function _jb2008_ΔTc(F10::Number, lst::Number, glat::Number, h::Number)
    # Auxiliary variables according to [2, p.  784].
    B  = _jb2008_B
    C  = _jb2008_C
    F  = (F10 - 100)/100
    θ  = lst/24
    θ² = θ^2
    θ³ = θ^3
    θ⁴ = θ^4
    θ⁵ = θ^5
    ϕ  = cos(glat)

    ΔTc = 0.0

    # Compute the temperature variation given the altitude.
    @inbounds if 120 <= h <= 200
        ΔTc200 = C[17]     + C[18]*θ*ϕ   + C[19]*θ²*ϕ   + C[20]*θ³*ϕ +
                 C[21]*F*ϕ + C[22]*θ*F*ϕ + C[23]*θ²*F*ϕ

        ΔTc200Δz = C[1]      + B[2]*F     + C[3]*θ*F    + C[4]*θ²*F    +
                   C[5]*θ³*F + C[6]*θ⁴*F  + C[7]*θ⁵*F   + C[8]*θ*ϕ     +
                   C[9]*θ²*ϕ + C[10]*θ³*ϕ + C[11]*θ⁴*ϕ  + C[12]*θ⁵*ϕ   +
                   C[13]*ϕ   + C[14]*F*ϕ  + C[15]*θ*F*ϕ + C[16]*θ²*F*ϕ

        zp  = (h - 120)/80
        ΔTc = (3ΔTc200 - ΔTc200Δz)*zp^2 + (ΔTc200Δz - 2ΔTc200)*zp^3

    elseif 200 < h <= 240
        H = (h - 200)/50

        ΔTc = C[ 1]*H      + B[2]*F*H     + C[ 3]*θ*F*H   + C[ 4]*θ²*F*H   +
              C[ 5]*θ³*F*H + C[6]*θ⁴*F*H  + C[ 7]*θ⁵*F*H  + C[ 8]*θ*ϕ*H    +
              C[9]*θ²*ϕ*H  + C[10]*θ³*ϕ*H + C[11]*θ⁴*ϕ*H  + C[12]*θ⁵*ϕ*H   +
              C[13]*ϕ*H    + C[14]*F*ϕ*H  + C[15]*θ*F*ϕ*H + C[16]*θ²*F*ϕ*H +
              C[17]        + C[18]*θ*ϕ    + C[19]*θ²*ϕ    + C[20]*θ³*ϕ     +
              C[21]*F*ϕ    + C[22]*θ*F*ϕ  + C[23]*θ²*F*ϕ

    elseif 240 < h <= 300
        H = 40/50

        aux1 = C[ 1]*H      + B[ 2]*F*H    + C[ 3]*θ*F*H   + C[ 4]*θ²*F*H   +
               C[ 5]*θ³*F*H + C[ 6]*θ⁴*F*H + C[ 7]*θ⁵*F*H  + C[ 8]*θ*ϕ*H    +
               C[ 9]*θ²*ϕ*H + C[10]*θ³*ϕ*H + C[11]*θ⁴*ϕ*H  + C[12]*θ⁵*ϕ*H   +
               C[13]*ϕ*H    + C[14]*F*ϕ*H  + C[15]*θ*F*ϕ*H + C[16]*θ²*F*ϕ*H +
               C[17]        + C[18]*θ*ϕ    + C[19]*θ²*ϕ    + C[20]*θ³*ϕ     +
               C[21]*F*ϕ    + C[22]*θ*F*ϕ  + C[23]*θ²*F*ϕ

        aux2 = C[ 1]        + B[ 2]*F    + C[ 3]*θ*F + C[ 4]*θ²*F + C[ 5]*θ³*F  +
               C[ 6]*θ⁴*F   + C[ 7]*θ⁵*F + C[ 8]*θ*ϕ + C[ 9]*θ²*ϕ + C[10]*θ³*ϕ  +
               C[11]*θ⁴*ϕ   + C[12]*θ⁵*ϕ + C[13]*ϕ   + C[14]*F*ϕ  + C[15]*θ*F*ϕ +
               C[16]*θ²*F*ϕ

        H = 300/100

        ΔTc300 = B[ 1]        + B[ 2]*F      + B[ 3]*θ*F    + B[ 4]*θ²*F   +
                 B[ 5]*θ³*F   + B[ 6]*θ⁴*F   + B[ 7]*θ⁵*F   + B[ 8]*θ*ϕ    +
                 B[ 9]*θ²*ϕ   + B[10]*θ³*ϕ   + B[11]*θ⁴*ϕ   + B[12]*θ⁵*ϕ   +
                 B[13]*H*ϕ    + B[14]*θ*H*ϕ  + B[15]*θ²*H*ϕ + B[16]*θ³*H*ϕ +
                 B[17]*θ⁴*H*ϕ + B[18]*θ⁵*H*ϕ + B[19]*ϕ

        ΔTc300Δz = B[13]*ϕ + B[14]*θ*ϕ + B[15]*θ²*ϕ + B[16]*θ³*ϕ + B[17]*θ⁴*ϕ +
                   B[18]*θ⁵*ϕ

        aux3 = 3ΔTc300 - ΔTc300Δz - 3aux1 - 2aux2
        aux4 = ΔTc300 - aux1 - aux2 - aux3
        zp   = (h - 240)/60
        ΔTc  = aux1 + aux2*zp + aux3*zp^2 + aux4*zp^3

    elseif 300 < h <= 600
         H  = h/100
        ΔTc = B[ 1]        + B[ 2]*F      + B[ 3]*θ*F    + B[ 4]*θ²*F   +
              B[ 5]*θ³*F   + B[ 6]*θ⁴*F   + B[ 7]*θ⁵*F   + B[ 8]*θ*ϕ    +
              B[ 9]*θ²*ϕ   + B[10]*θ³*ϕ   + B[11]*θ⁴*ϕ   + B[12]*θ⁵*ϕ   +
              B[13]*H*ϕ    + B[14]*θ*H*ϕ  + B[15]*θ²*H*ϕ + B[16]*θ³*H*ϕ +
              B[17]*θ⁴*H*ϕ + B[18]*θ⁵*H*ϕ + B[19]*ϕ

    elseif 600 < h <= 800
        zp   = (h - 600)/100
        hp   = 600/100
        aux1 = B[ 1]         + B[ 2]*F       + B[ 3]*θ*F     + B[ 4]*θ²*F    +
               B[ 5]*θ³*F    + B[ 6]*θ⁴*F    + B[ 7]*θ⁵*F    + B[ 8]*θ*ϕ     +
               B[ 9]*θ²*ϕ    + B[10]*θ³*ϕ    + B[11]*θ⁴*ϕ    + B[12]*θ⁵*ϕ    +
               B[13]*hp*ϕ    + B[14]*θ*hp*ϕ  + B[15]*θ²*hp*ϕ + B[16]*θ³*hp*ϕ +
               B[17]*θ⁴*hp*ϕ + B[18]*θ⁵*hp*ϕ + B[19]*ϕ

        aux2 = B[13]*ϕ    + B[14]*θ*ϕ + B[15]*θ²*ϕ + B[16]*θ³*ϕ + B[17]*θ⁴*ϕ +
               B[18]*θ⁵*ϕ

        aux3 = -(3aux1 + 4aux2)/4
        aux4 =  ( aux1 +  aux2)/4
        ΔTc  = aux1 + aux2*zp + aux3*zp^2 + aux4*zp^3
    end

    ΔTc
end

"""
    function _jb2008_δf1(z, Tx, T∞)

Auxiliary function to compute the integrand in `_jb2008_int`.

"""
@inline function _jb2008_δf1(z, Tx, T∞)
    Mb = _jb2008_M(z)
    Tl = _jb2008_T(z, Tx, T∞)
    g  = _jb2008_grav(z)
    Mb * g / Tl
end

"""
    function _jb2008_δf2(z, Tx, T∞)

Auxiliary function to compute the integrand in `_jb2008_int`.

"""
@inline function _jb2008_δf2(z, Tx, T∞)
    Tl = _jb2008_T(z, Tx, T∞)
    g  = _jb2008_grav(z)
    g / Tl
end

