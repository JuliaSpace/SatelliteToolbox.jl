#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   The Jacchia-Bowman 2008 Atmospheric Model, a product of the Space
#   Environment Technologies.
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

export jb2008

"""
    function jb2008(JD::Number, glat::Number, glon::Number, h::Number)
    function jb2008(JD::Number, glat::Number, glon::Number, h::Number, F10::Number, F10ₐ::Number, S10::Number, S10ₐ::Number, M10::Number, M10ₐ::Number, Y10::Number, Y10ₐ::Number, DstΔTc::Number)

Compute the atmospheric density using the Jacchia-Bowman 2008 (JB2008) model.

If the space indices are not provided (first call), then they will be obtained
from the online database. In this case, the function `init_space_indices()` must
be called first and the function will throw an exception if the selected `JD` is
outside of the available data.

This model is a product of the **Space Environment Technologies**, more
information can be seen in the websites:

http://sol.spacenvironment.net/jb2006/

http://sol.spacenvironment.net/jb2008/

# Args

* `JD`: Julian day.
* `glat`: Geocentric latitude [rad].
* `glon`: Geocentric longitude [rad].
* `h`: Altitude [m].

* `F10`: 10.7-cm solar flux \\[10⁻²² W/(M² Hz)] (Tabular time 1 day earlier).
* `F10ₐ`: 10.7-cm averaged solar flux, 81-day centered on input time (Tabular
          time 1 day earlier).
* `S10`: EUV index (26-34 nm) scaled to F10.7 (Tabular time 1 day earlier).
* `S10ₐ`: EUV 81-day averaged centered index (Tabular time 1 day earlier).
* `M10`: MG2 index scaled to F10.7 (Tabular time 2 days earlier).
* `M10ₐ`: MG2 81-day averaged centered index (Tabular time 2 days earlier).
* `Y10`: Solar X-ray & Lya index scaled to F10.7 (Tabular time 5 days earlier).
* `Y10ₐ`: Solar X-ray & Lya 81-day averaged centered index (Tabular time 5 days
          earlier).
* `DstΔTc`: Temperature variation related to the Dst.

# Returns

An instance of the structure `JB2008_Output` with the computed values.

"""
function jb2008(JD::Number, glat::Number, glon::Number, h::Number)
    # Get the data in the desired Julian Day considering the tabular time of the
    # model.
    F10    = get_F10(JD-1)
    F10ₐ   = get_F10ₐ(JD-1)
    S10    = get_S10(JD-1)
    S10ₐ   = get_S10ₐ(JD-1)
    M10    = get_M10(JD-2)
    M10ₐ   = get_M10ₐ(JD-2)
    Y10    = get_Y10(JD-5)
    Y10ₐ   = get_Y10ₐ(JD-5)
    DstΔTc = get_DstΔTc(JD)

    jb2008(JD, glat, glon, h, F10, F10ₐ, S10, S10ₐ, M10, M10ₐ, Y10, Y10ₐ, DstΔTc)
end

function jb2008(JD::Number, glat::Number, glon::Number, h::Number, F10::Number,
                F10ₐ::Number, S10::Number, S10ₐ::Number, M10::Number,
                M10ₐ::Number, Y10::Number, Y10ₐ::Number, DstΔTc::Number)

    # Constants
    # ==========================================================================

    T₁    = 183.0       # Temperature at the lower bound [K].
    z₁    = 90.0        # Altitude of the lower bound [km].
    zx    = 125.0       # Altitude of the inflection point [km].
    Rstar = 8314.32     # Rstar is the universal gas-constant (mks)
                        # [joules/(K.kmol)].
    ρ₁    = 3.46e-6     # Density at `z₁` [kg/m³].
    A     = 6.02257e26  # Avogadro's constant (mks) [molecules/kmol].

    # Assumed sea-level composition.
    Mb₀  = 28.960
    q₀N₂ = 0.78110
    q₀O₂ = 0.20955
    q₀Ar = 9.3400e-3
    q₀He = 1.2890e-5

    # Molecular weights of each specie [kg/kmol].
    MN₂ = 28.0134
    MO₂ = 31.9988
    MO  = 15.9994
    MAr = 39.9480
    MHe =  4.0026
    MH  =  1.00797

    # Thermal diffusion coefficient for each specie.
    α_N₂ =  0.0
    α_O₂ =  0.0
    α_O  =  0.0
    α_Ar =  0.0
    α_He = -0.38
    α_H₂ =  0.0

    # Values used to stablish height step sizes in the integration process
    # between 90km to 105km, 105km to 500km, and above 500km.
    R1 = 0.010
    R2 = 0.025
    R3 = 0.075

    # Preliminaries
    # ==========================================================================

    # Convert the altitude from [m] to [km].
    h /= 1000

    # Compute the Sun position represented in the inertial reference frame
    # (MOD).
    s_i = sun_position_i(JD)

    # Compute the Sun declination [rad].
    δs = atan( s_i[3], sqrt(s_i[1]*s_i[1] + s_i[2]*s_i[2]) )

    # Compute the Sun right ascension [rad].
    Ωs = atan( s_i[2], s_i[1] )

    # Compute the right ascension of the selected location w.r.t. the inertial
    # reference frame.
    Ωp = glon + JDtoGMST(JD)

    # Compute the hour angle at the selected location, which is the angle
    # measured at the XY plane between the right ascension of the selected
    # position and the right ascension of the Sun.
    H = Ωp - Ωs

    # Compute the local solar time.
    lst = mod( (H + π)*12/π, 24 )

    # Algorithm
    # ==========================================================================

    # Eq. 2 [1], Eq. 14 [3]
    # --------------------------------------------------------------------------
    #
    # Nighttime minimum of the global exospheric temperature distribution when
    # the planetary geomagnetic index Kp is zero.
    ΔF10 =  F10 - F10ₐ
    ΔS10 =  S10 - S10ₐ
    ΔM10 = M10 - M10ₐ
    ΔY10 =  Y10 - Y10ₐ

    Wt = (F10ₐ/240)^(1/4)
    (Wt > 1) && (Wt = 1.0)
    Fsₐ = F10ₐ*Wt + S10ₐ*(1-Wt)
    Tc = 392.4 + 3.227Fsₐ + 0.298ΔF10 + 2.259ΔS10 + 0.312ΔM10 + 0.178ΔY10

    # Eq. 15 [3]
    # --------------------------------------------------------------------------
    η = abs(glat - δs)/2
    θ = abs(glat + δs)/2

    # Eq. 16 [3]
    # --------------------------------------------------------------------------
    τ = H - 0.64577182 + 0.10471976sin(H + 0.75049158)

    # Eq. 17 [3]
    # --------------------------------------------------------------------------
    m = 2.5
    n = 3.0
    R = 0.31
    C = cos(η)^m
    S = sin(θ)^m

    # NOTE: The original equation in [3] does not have the `abs` as in the
    # source-code of JB2008.
    Tl = Tc*( 1 + R*( S + (C-S)*abs( cos(τ/2) )^n ) )

    # Compute the correction to `Tc` considering the local solar time and
    # latitude.
    ΔTc = _jb2008_ΔTc(F10, lst, glat, h)

    # Compute the local exospheric temperature with the geomagnetic storm
    # effect.
    T_exo = Tl + DstΔTc
    T∞    = T_exo + ΔTc

    # Eq. 9 [3]
    # --------------------------------------------------------------------------
    #
    # Temperature at the inflection point `z = 125 km`.
    a  =  444.3807
    b  =    0.02385
    c  = -392.8292
    k  =   -0.0021357
    Tx = a + b*T∞ + c*exp(k*T∞)

    # Eq. 5 [3]
    # --------------------------------------------------------------------------
    #
    # From 90 to 105 km, for a given temperature profile `T[k]`, the density
    # ρ is computed by integrating the barometric equation:
    #
    #                  - _ -     _
    #                 |  M´ |    M´ g
    #   d lnρ´ = d ln | --- | - ------ dz ,
    #                 |  T  |    R* T
    #                  -   -
    #
    # which can be rewritten as:
    #
    #         -     -       _
    #        |  ρ´T  |      M´ g
    #   d ln | ----- | = - ------ dz .
    #        |   _   |      R* T
    #        |   M´  |
    #         -     -
    #
    #
    # Here, we will compute the following integral:
    #
    #
    #      / z₂
    #     |     _
    #     |     M´g
    #     |   ------- dz
    #     |      T
    #     |
    #    /   z₁
    #
    # in which z₁ is the minimum value between `h` and 105 km. The integration
    # will be computed by Newton-Cotes 4th degree method.

    z₂ = min(h, 105.0)

    int, z₂ = _jb2008_int(z₁, z₂, R1, (z)->begin
                              Mb = _jb2008_M(z)
                              Tl = _jb2008_T(z, Tx, T∞)
                              g  = _jb2008_grav(z)
                              Mb * g / Tl
                          end::Float64)

    Mb₁ = _jb2008_M(z₁)
    Tl₁ = _jb2008_T(z₁, Tx, T∞)
    Mb₂ = _jb2008_M(z₂)
    Tl₂ = _jb2008_T(z₂, Tx, T∞)

    # `Mbj` and `Tlj` contains, respectively, the mean molecular mass and local
    # temperature at the boundary of the integration interval.
    #
    # The `1000` factor is to convert `Rstar` to the appropriate units.
    ρ = ρ₁*(Mb₂/Mb₁)*(Tl₁/Tl₂)*exp(-1000int/Rstar)

    # Eq. 2 [3]
    # --------------------------------------------------------------------------
    NM = A*ρ
    N  = NM/Mb₂

    # Eq. 3 [3]
    # --------------------------------------------------------------------------
    NMoMb₀  = NM/Mb₀
    log_nN₂ = log(q₀N₂*NMoMb₀)
    log_nAr = log(q₀Ar*NMoMb₀)
    log_nHe = log(q₀He*NMoMb₀)

    # Eq. 4 [3]
    # --------------------------------------------------------------------------
    log_nO₂ = log(NMoMb₀ * (1 + q₀O₂) - N)
    log_nO  = log(2(N - NMoMb₀))
    log_nH  = 0.0

    Tz = 0.0
    if h <= 105
        Tz = Tl₂
        log_nH₂ = log_nHe - 25
    else
        # Eq.6 [3]
        # ----------------------------------------------------------------------
        #
        # From 100 km to 500 km, neglecting the thermal diffusion coefficient,
        # the eq. 6 in [3] can be written as:
        #
        #                          _
        #               1+αᵢ       M´ g
        #   d ln ( n . T    ) = - ------ dz ,
        #                          R* T
        #
        # where n number density for the i-th specie. This equations must be
        # integrated for each specie we are considering.
        #
        # Here, we will compute the following integral:
        #
        #
        #      / z₃
        #     |
        #     |     g
        #     |   ----- dz
        #     |     T
        #     |
        #    /   z₂
        #
        # in which z₃ is the minimum value between `h` and 500 km. The
        # integration will be computed by Newton-Cotes 4th degree method.

        z₃ = min(h, 500.0)

        int_1, z₃ = _jb2008_int(z₂, z₃, R1, (z)->begin
                                    Tl = _jb2008_T(z, Tx, T∞)
                                    g  = _jb2008_grav(z)
                                    g / Tl
                                end::Float64)

        Tl₃ = _jb2008_T(z₃, Tx, T∞)

        # If `h` is lower than 500 km, then keep integrating until 500 km due to
        # the hydrogen. Otherwise, continue the integration until `h`. Hence, we
        # will compute the following integral:
        #
        #      / z₄
        #     |
        #     |     g
        #     |   ----- dz
        #     |     T
        #     |
        #    /   z₃
        #
        # in which z₄ is the maximum value between `h` and 500 km. The
        # integration will be computed by Newton-Cotes 4th degree method.

        z₄ = max(h, 500.0)

        int_2, z₄ = _jb2008_int(z₃, z₄, (h <= 500) ? R2 : R3, (z)->begin
                                    Tl = _jb2008_T(z, Tx, T∞)
                                    g  = _jb2008_grav(z)
                                    g / Tl
                                end::Float64)

        Tl₄ = _jb2008_T(z₄, Tx, T∞)

        if h <= 500
            T500km    = Tl₄
            Tz        = Tl₃
            log_TfoTi = log(Tl₃/Tl₂)
            H_sign    = +1

            #          g
            # goRT = ------
            #         R*T
            goRT      = 1000int_1/Rstar
        else
            T500km    = Tl₃
            Tz        = Tl₄
            log_TfoTi = log(Tl₄/Tl₂)
            H_sign    = -1

            #          g
            # goRT = ------
            #         R*T
            goRT      = 1000(int_1 + int_2)/Rstar
        end

        log_nN₂ += -(1 + α_N₂)*log_TfoTi - goRT*MN₂
        log_nO₂ += -(1 + α_O₂)*log_TfoTi - goRT*MO₂
        log_nO  += -(1 + α_O )*log_TfoTi - goRT*MO
        log_nAr += -(1 + α_Ar)*log_TfoTi - goRT*MAr
        log_nHe += -(1 + α_He)*log_TfoTi - goRT*MHe

        # Eq. 7 [3]
        # ----------------------------------------------------------------------
        #
        # The equation of [3] will be converted from `log10` to `log`.
        # Furthermore, the units will be converted to [1/cm³] to [1/m³].

        log10_nH_500km = @evalpoly(log10(T∞), 73.13, -39.40, 5.5)
        log_nH_500km   = log(10) * ( log10_nH_500km + 6 )
        #                                             ^
        #                                             |
        #                         This factor converts from [1/cm³] to [1/m³]

        # Compute the hydrogen density based on the value at 500km.
        log_nH = log_nH_500km + H_sign *
            ( log(Tl₄/Tl₃) + ( 1000/Rstar ) * int_2 * MH )
    end

    # Eq. 24 [3] - Seasonal-latitudinal variation
    # ----------------------------------------------------------------------
    #
    # TODO: The term related to the year in the source-code of JB 2008 is
    # different from [3]. This must be verified.

    #         | Modified JD  |
    Φ = mod( ( JD - 2400000.5 - 36204 )/365.2422, 1 )

    Δlog₁₀ρ = 0.02 * ( h - 90 ) * sign(glat) * exp( -0.045( h - 90 ) ) *
        sin(glat)^2 * sin( 2π*Φ + 1.72 )

    # Convert from `log10` to `log`.
    Δlogρ = log(10)*Δlog₁₀ρ

    # Eq. 23 [3] - Semiannual variation
    # ----------------------------------------------------------------------

    if h < 2000
        # Compute the year given the selected Julian Day.
        year, month, day, = JDtoDate(JD)

        # Compute the day of the year.
        doy = JD - DatetoJD(year, 1, 1, 0, 0, 0) + 1

        # Use the new semiannual model from [1].
        Fz, Gz, Δsalog₁₀ρ = _jb2008_semiannual(doy, h, F10ₐ, S10ₐ, M10ₐ)

        (Fz < 0) && (Δsalog₁₀ρ = 0.0)

        # Convert from `log10` to `log`.
        Δsalogρ = log(10)*Δsalog₁₀ρ
    else
        Δsalogρ = 0.0
    end

    # Compute the total variation.
    Δρ = Δlogρ + Δsalogρ

    # Apply to the number densities.
    log_nN₂ += Δρ
    log_nO₂ += Δρ
    log_nO  += Δρ
    log_nAr += Δρ
    log_nHe += Δρ
    log_nH  += Δρ

    # Compute the high altitude exospheric density correction factor.
    #
    # The source-code computes this correction after computing the mass
    # density and mean molecular weight. Hence, the correction is applied
    # only to the total density. However, we will apply the correction to
    # each specie. TODO: Verify if this is reasonable.
    log_FρH  = log(_jb2008_highaltitude(h, F10ₐ))
    log_nN₂ += log_FρH
    log_nO₂ += log_FρH
    log_nO  += log_FρH
    log_nAr += log_FρH
    log_nHe += log_FρH
    log_nH  += log_FρH

    # Compute the mass density and mean molecular weight and convert number
    # density logs from natural to common.
    sum_n  = 0.0
    sum_mn = 0.0

    for (log_n, M) in [ (log_nN₂, MN₂), (log_nO₂, MO₂), (log_nO, MO),
                        (log_nAr, MAr), (log_nHe, MHe), (log_nH, MH) ]
        n = exp(log_n)
        sum_n  += n
        sum_mn += n*M
    end

    nN₂ = exp(log_nN₂)
    nO₂ = exp(log_nO₂)
    nO  = exp(log_nO)
    nAr = exp(log_nAr)
    nHe = exp(log_nHe)
    nH  = exp(log_nH)
    ρ   = sum_mn/A

    # Create and return the output structure.
    JB2008_Output{Float64}(nN₂, nO₂, nO, nAr, nHe, nH, ρ, T_exo, Tz)
end

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
    if 1000 <= h <= 1500
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
    function _jb2008_int(z₀::Number, z₁::Number, R::Number, δf::Function)

Compute the integral of the function `δf` between `z₀` and `z₁` using the
Newton-Cotes 4th degree method. `R` is a number that defined the step size.

This function returns a tuple containing the integral and last value of `z` used
in the numerical algorithm.

"""
function _jb2008_int(z₀::Number, z₁::Number, R::Number, δf::Function)
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
    for i = 1:n
        zi₀   = zi₁           # The beginning of the i-th integration step.
        zi₁   = zr*zi₁        # The end of the i-th integration step.
        Δz    = (zi₁ - zi₀)/4 # Step for the i-th integration step.
        int_i = WT[1]*δf(zi₀) # First term of the Newton-Cotes 4th degree sum.

        # Compute the Newton-Cotes 4th degree sum.
        for j = 2:5
            zj  = zi₀ + (j-1)*Δz

            # Value of the integrand at `zj`.
            int_i += WT[j] * δf(zj)
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
    Fz = B[1] + B[2]*Fsmjₐ + B[3]*z*Fsmjₐ + B[4]*z^2*Fsmjₐ + B[5]*z*Fsmjₐ^2

    # Compute the new 81-day centered solar index for G(t) according to eq. 6
    # [1].
    Fsmₐ = 1F10ₐ - 0.75S10ₐ - 0.37M10ₐ

    # Compute the semiannual G(t) yearly periodic function according to eq. 7
    # [1].
    sω,  cω  = sincos(1ω)
    s2ω, c2ω = sincos(2ω)

    Gt =        C[1] + C[2]*sω + C[3]*cω + C[4]*s2ω + C[5]*c2ω +
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
    if 120 <= h <= 200
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
