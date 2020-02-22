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

include("./jb2008_priv.jl")

export jb2008

"""
    jb2008(JD::Number, glat::Number, glon::Number, h::Number)
    jb2008(JD::Number, glat::Number, glon::Number, h::Number, F10::Number, F10ₐ::Number, S10::Number, S10ₐ::Number, M10::Number, M10ₐ::Number, Y10::Number, Y10ₐ::Number, DstΔTc::Number)

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
    F10    = get_space_index(Val{:F10obs},  JD-1)
    F10ₐ   = get_space_index(Val{:F10Mobs}, JD-1)
    S10    = get_space_index(Val{:S10},     JD-1)
    S10ₐ   = get_space_index(Val{:S81a},    JD-1)
    M10    = get_space_index(Val{:M10},     JD-2)
    M10ₐ   = get_space_index(Val{:M81a},    JD-2)
    Y10    = get_space_index(Val{:Y10},     JD-5)
    Y10ₐ   = get_space_index(Val{:Y81a},    JD-5)
    DstΔTc = get_space_index(Val{:DstΔTc},  JD)

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

    int, z₂ = _jb2008_int(z₁, z₂, R1, Tx, T∞, _jb2008_δf1)

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

        int_1, z₃ = _jb2008_int(z₂, z₃, R1, Tx, T∞, _jb2008_δf2)

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

        int_2, z₄ = _jb2008_int(z₃, z₄, (h <= 500) ? R2 : R3, Tx, T∞, _jb2008_δf2)

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

