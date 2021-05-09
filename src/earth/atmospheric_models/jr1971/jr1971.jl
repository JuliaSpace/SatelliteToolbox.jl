# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   The Jacchia-Roberts 1971 Atmospheric Model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Roberts, C. R (1971). An analytic model for upper atmosphere densities
#       based upon Jacchia's 1970 models.
#
#   [2] Jacchia, L. G (1970). New static models of the thermosphere and
#       exosphere with empirical temperature profiles. SAO Special Report #313.
#
#   [3] Long, A. C., Cappellari Jr., J. O., Velez, C. E., Fuchs, A. J (editors)
#       (1989). Goddard Trajectory Determination System (GTDS) Mathematical
#       Theory (Revision 1). FDD/552-89/0001 and CSC/TR-89/6001.
#
#   [4] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [5] GMT-2668 - http://bugs.gmatcentral.org/browse/GMT-2668
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export jr1971

include("./jr1971_priv.jl")

"""
    jr1971(JD::Number, glat::Number, glon::Number, h::Number, F10::Number, F10ₐ::Number, Kp::Number)

Compute the atmospheric density using the Jacchia-Roberts 1971 model.

# Args

* `JD`: Julian day.
* `glat`: Geodetic latitude [rad].
* `glon`: Geodetic longitude [rad].
* `h`: Altitude [m].
* `F10`: 10.7-cm solar flux \\[10⁻²² W/(M² Hz)].
* `F10ₐ`: 10.7-cm averaged solar flux, 81-day centered on input time.
* `Kp`: Kp geomagnetic index (with a delay of 3 hours).

# Returns

An instance of the structure `JR1971_Output` with the computed values.

"""
function jr1971(JD::Number, glat::Number, glon::Number, h::Number, F10::Number,
                F10ₐ::Number, Kp::Number)

    # Constants
    # ==========================================================================

    @unpack_JR1971_CONSTANTS _jr1971_constants

    # Auxiliary variables
    # ==========================================================================

    Ra² = Ra*Ra       # Mean Earth radius squared [km²].

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
    Ωp = glon + jd_to_gmst(JD)

    # Compute the hour angle at the selected location, which is the angle
    # measured at the XY plane between the right ascension of the selected
    # position and the right ascension of the Sun.
    H = Ωp - Ωs

    # Algorithm
    # ==========================================================================

    # Exospheric temperature
    # ==========================================================================

    # Diurnal variation
    # =================

    # Eq. 14 [2], eq. 4-89 [3]
    # --------------------------------------------------------------------------
    #
    # Nighttime minimum of the global exospheric temperature distribution when
    # the planetary geomagnetic index Kp is zero.
    ΔF10 =  F10 - F10ₐ
    Tc   = 379 + 3.24F10 + 1.3ΔF10

    # Eq. 15 [2], eq. 4-91 [3]
    # --------------------------------------------------------------------------
    η = abs(glat - δs)/2
    θ = abs(glat + δs)/2

    # Eq. 16 [2], eq. 4-91 [3]
    # --------------------------------------------------------------------------
    τ = H + ( -37 + 6sin(H + 43*π/180) )*π/180

    # Eq. 17 [2], eq. 4.90 [3]
    # --------------------------------------------------------------------------
    m = 2.2
    n = 3.0
    R = 0.3
    C = cos(η)^m
    S = sin(θ)^m

    Tl = Tc*( 1 + R*( S + (C-S)*cos(τ/2)^n ) )

    # Variations with geomagnetic activity
    # ====================================

    # Eq. 18 or Eq. 20 [2], eq. 4-94 [3]
    # --------------------------------------------------------------------------

    ΔT∞ = (h < 200) ? 14Kp + 0.02exp(Kp) : 28Kp + 0.03exp(Kp)

    # Eq. 4-95 [3]
    # --------------------------------------------------------------------------
    #
    # Compute the local exospheric temperature with the geomagnetic storm
    # effect.
    T∞ = Tl + ΔT∞

    # Temperature at desired altitude
    # ==========================================================================

    # Eq. 4-96 [3]
    # --------------------------------------------------------------------------
    #
    # Compute the temperature at inflection point `zx`.
    #
    # The values at [1, p. 369] are from an old version of Jacchia 1971 model.
    # We will use the new values available at [2].
    a  =  371.6678
    b  =    0.0518806
    c  = -294.3505
    d  =   -0.00216222
    Tx = a + b*T∞ + c*exp(d*T∞)

    # Compute the temperature at desired point.
    Tz = _jr1971_T(h, Tx, T∞)

    # Compute corrections to density, eqs. 4-96 to 4-101 [3]
    # ==========================================================================

    # If the user wants, then compute the corrections of geomegnetic effect,
    # semi-annual variation, and seasonal latitudinal variation. Otherwise, this
    # function will return the standard density [4].

    Δlog₁₀ρ_g  = 0.0
    Δlog₁₀ρ_sa = 0.0
    Δlog₁₀ρ_lt = 0.0
    Δlog₁₀ρ_c  = 0.0

    # Geomagnetic effect, eq. 4-97 [3]
    # --------------------------------------------------------------------------
    (h < 200) && ( Δlog₁₀ρ_g = 0.012Kp + 1.2e-5exp(Kp) )

    # Semi-annual variation, eqs. 4-98 to 4-99 [3]
    # --------------------------------------------------------------------------
    Φ    = (JD - 2436204.5)/365.2422  # Number of days since January 1, 1958.

    τ_sa = Φ + 0.09544( ( 1/2 + 1/2*sin( 2π*Φ + 6.035 ) )^(1.65) - 1/2 )
    f_z  = ( 5.876e-7h^2.331 + 0.06328 )*exp( -0.002868h )
    g_t  = 0.02835 + ( 0.3817 + 0.17829sin( 2π*τ_sa + 4.137 ) )*
    sin( 4π*τ_sa + 4.259 )

    Δlog₁₀ρ_sa = f_z*g_t

    # Seasonal latitudinal variation, eq. 4-100 [3]
    # --------------------------------------------------------------------------
    Δlog₁₀ρ_lt = 0.014( h - 90 )*exp( -0.0013( h - 90 )^2 )*
    sin( 2π*Φ + 1.72 )*sin( glat )*abs( sin(glat) )

    # Total correction, eq. B-10 [4]
    # --------------------------------------------------------------------------
    Δlog₁₀ρ_c = Δlog₁₀ρ_g + Δlog₁₀ρ_lt + Δlog₁₀ρ_sa

    # Density
    # ==========================================================================

    if h == z₁
        ρ = ρ₁ * 10^Δlog₁₀ρ_c

        # Convert to SI and return.
        return JR1971_Output(nN2   = (ρ * μi.N₂) * Av/Mi.N₂ * 1e6,
                             nO2   = (ρ * μi.O₂) * Av/Mi.O₂ * 1e6,
                             nO    = (ρ * μi.O)  * Av/Mi.O  * 1e6,
                             nAr   = (ρ * μi.Ar) * Av/Mi.Ar * 1e6,
                             nHe   = (ρ * μi.He) * Av/Mi.He * 1e6,
                             nH    = (ρ * μi.H)  * Av/Mi.H  * 1e6,
                             rho   = ρ*1e3,
                             T_exo = T∞,
                             Tz    = Tz )

    elseif z₁ < h <= zx

        # First, we need to find the roots of the polynomial:
        #
        #   P(Z) = c₀ + c₁ ⋅ z + c₂ ⋅ z² + c₃ ⋅ z³ + c₄ ⋅ z⁴
        c₀ = ( 35^4*Tx/(Tx - T₁) + Ca[1] )/Ca[5]
        c₁ = Ca[2]/Ca[5]
        c₂ = Ca[3]/Ca[5]
        c₃ = Ca[4]/Ca[5]
        c₄ = Ca[5]/Ca[5]
        r₁, r₂, x, y = _jr1971_roots([c₀, c₁, c₂, c₃, c₄])

        # f and k, [1. p. 371]
        # ----------------------------------------------------------------------
        f = 35^4*Ra²/Ca[5]
        k = -g₀/( Rstar*(Tx - T₁) )

        # U(ν), V(ν), and W(ν) functions, [1, p. 372]
        # ----------------------------------------------------------------------

        # TODO: Create copies of the roots, which are the variables used in the
        # functions U(ν), V(ν), and W(ν). Otherwise, the code is not type
        # stable. This seems related to:
        #
        #   https://github.com/JuliaLang/julia/issues/15276
        #
        cr₁, cr₂, cx, cy = r₁, r₂, x, y

        U(ν) = (ν + Ra)^2*(ν^2 - 2cx*ν + cx^2 + cy^2)*(cr₁ - cr₂)
        V(ν) = (ν^2 - 2cx*ν + cx^2 + cy^2)*(ν - cr₁)*(ν - cr₂)

        # The equation for `W(ν)` in [1] was:
        #
        #   W(ν) = r₁*r₂*Ra²*(ν + Ra) + (x^2 + y^2)*Ra*(Ra*ν + r₁*r₂)
        #
        # However, [3,4] mention that this is not correct, and must be replaced
        # by:
        W(ν) = cr₁*cr₂*Ra*(Ra + ν)*( Ra + ( cx^2 + cy^2 )/ν )

        # X, [1, p. 372]
        # ------------------------------------------------------------------
        X = -2r₁*r₂*Ra*( Ra² + 2x*Ra + x^2 + y^2)

        # The original paper [1] divided into two sections: 90 km to 105 km, and
        # 105 km to 125 km. However, [3] divided into 90 to 100 km, and 100 km
        # to 125 km.

        if h <= z₂

            # Altitudes between 90 km and 100 km
            # ==================================================================

            # S(z) polynomial, [1, p. 371]
            # ------------------------------------------------------------------
            B₀, B₁, B₂, B₃, B₄, B₅ = [ α[k] + β[k]*Tx/(Tx - T₁) for k = 1:6 ]

            S(z) = @evalpoly(z, B₀, B₁, B₂, B₃, B₄, B₅)

            # Auxiliary variables, [1. p. 372]
            # ------------------------------------------------------------------
            p₂ =   S(r₁)/U(r₁)
            p₃ =  -S(r₂)/U(r₂)
            p₅ =  S(-Ra)/V(-Ra)

            # There is a typo in the fourth term in [1] that was corrected in
            # [3].
            p₄ = ( B₀ - r₁*r₂*Ra²*( B₄ + ( 2x + r₁ + r₂ - Ra )*B₅ ) -
                  r₁*r₂*Ra*( x^2 + y^2 )*B₅ + r₁*r₂*( Ra² - ( x^2 + y^2 ) )*p₅ +
                  W(r₁)*p₂ + W(r₂)*p₃ )/X

            p₆ = B₄ + ( 2x + r₁ + r₂ - Ra )*B₅ - p₅ - 2( x + Ra )*p₄ -
                 ( r₂ + Ra )*p₃ - (r₁ + Ra)*p₂

            p₁ = B₅ - 2p₄ - p₃ - p₂

            # F₁ and F₂, [1, p. 372-373]
            # ------------------------------------------------------------------
            log_F₁ = p₁ * log( ( h + Ra )/( z₁ + Ra ) ) +
                     p₂ * log( ( h - r₁ )/( z₁ - r₁ ) ) +
                     p₃ * log( ( h - r₂ )/( z₁ - r₂ ) ) +
                     p₄ * log( (  h^2 - 2x*h  + x^2 + y^2 )/
                               ( z₁^2 - 2x*z₁ + x^2 + y^2 ) )

            # This equation in [4] is wrong, since `f` is multiplying `A₆`. We
            # will use the one in [3].
            F₂ = ( h - z₁ )*( Aa[7] + p₅/( (h + Ra)*(z₁ + Ra) ) ) +
                 p₆/y*atan( y*(h - z₁) / ( y^2 + (h - x)*(z₁ - x) ) )

            # Compute the density, eq. 13 [1]
            # ------------------------------------------------------------------
            Mz = _jr1971_M(h)
            ρ  = ρ₁*Mz*T₁/( M₁*Tz )*exp( k*( log_F₁ + F₂ ) ) * 10^Δlog₁₀ρ_c

            # Convert to SI and return.
            return JR1971_Output(nN2   = (ρ * μi.N₂) * Av/Mi.N₂ * 1e6,
                                 nO2   = (ρ * μi.O₂) * Av/Mi.O₂ * 1e6,
                                 nO    = (ρ * μi.O)  * Av/Mi.O  * 1e6,
                                 nAr   = (ρ * μi.Ar) * Av/Mi.Ar * 1e6,
                                 nHe   = (ρ * μi.He) * Av/Mi.He * 1e6,
                                 nH    = (ρ * μi.H)  * Av/Mi.H  * 1e6,
                                 rho   = ρ*1e3,
                                 T_exo = T∞,
                                 Tz    = Tz )
        else

            # Altitudes between 100 km and 125 km
            # ==================================================================

            # First, we need to compute the temperature and density at 100 km.
            T_₁₀₀ = _jr1971_T(z₂, Tx, T∞)

            # References [3,4] suggest to compute the density using a polynomial
            # fit, so that the computational burden can be reduced:
            #
            ρ_₁₀₀ = @evalpoly(T∞, ζ[1], ζ[2], ζ[3], ζ[4], ζ[5], ζ[6], ζ[7])*M₀

            # However, it turns out that this approach leads to discontinuity at
            # 100 km. This was also seen by GMAT [5].
            #
            # TODO: Check if we need to fix this.

            # Auxiliary variables, [1, p. 374]
            # ------------------------------------------------------------------

            q₂ =  1/U(r₁)
            q₃ = -1/U(r₂)
            q₅ =  1/V(-Ra)
            q₄ = (1 + r₁*r₂*(Ra² - (x^2 + y^2))*q₅ + W(r₁)*q₂ + W(r₂)*q₃)/X
            q₆ = -q₅ - 2(x + Ra)*q₄ - (r₂ + Ra)*q₃ - (r₁ + Ra)*q₂
            q₁ = -2q₄ - q₃ - q₂

            # F₃ and F₄, [1, p. 374]
            # ------------------------------------------------------------------
            log_F₃ = q₁ * log( ( h + Ra )/( z₂ + Ra ) ) +
                     q₂ * log( ( h - r₁ )/( z₂ - r₁ ) ) +
                     q₃ * log( ( h - r₂ )/( z₂ - r₂ ) ) +
                     q₄ * log( ( h^2 - 2x*h + x^2 + y^2 )/
                               ( z₂^2 - 2x*z₂ + x^2 + y^2 ) )

            F₄ = ( q₅*( h - z₂ ) )/( ( h + Ra )*( Ra + z₂ ) ) +
                 q₆/y*atan( y*(h - z₂)/( y^2 + (h - x)*(z₂ - x) ) )

            # Compute the density of each specie [3]
            # ------------------------------------------------------------------

            ρi = zeros(5)

            for i = 1:5
                ρi[i] = ρ_₁₀₀*Mi[i]/M₀*μi[i]*( T_₁₀₀/Tz )^(1 + αi[i])*
                        exp( Mi[i]*k*f*(log_F₃ + F₄ ) )
            end

            # Apply the corrections.
            ρi .*= 10^Δlog₁₀ρ_c

            # Convert to SI and return.
            return JR1971_Output(nN2   = ρi[1] * Av/Mi.N₂ * 1e6,
                                 nO2   = ρi[2] * Av/Mi.O₂ * 1e6,
                                 nO    = ρi[3] * Av/Mi.O  * 1e6,
                                 nAr   = ρi[4] * Av/Mi.Ar * 1e6,
                                 nHe   = ρi[5] * Av/Mi.He * 1e6,
                                 nH    = 0.0,
                                 rho   = sum(ρi)*1e3,
                                 T_exo = T∞,
                                 Tz    = Tz )
        end

    else

        # Altitudes higher than 125 km
        # ======================================================================

        # First, we need to compute the density at 125 km.

        # References [3,4] suggest to compute the density using a polynomial
        # fit, so that the computational burden can be reduced:

        ρi_₁₂₅ = zeros(5)

        for i = 1:5
            aux = @evalpoly(T∞, δij[i][1], δij[i][2], δij[i][3], δij[i][4],
                                δij[i][5], δij[i][6], δij[i][7])

            ρi_₁₂₅[i] = Mi[i]*10^aux/Av
        end

        # However, it turns out that this approach leads to discontinuity at 100
        # km. This was also seen by GMAT [5].
        #
        # TODO: Check if we need to fix this.

        # Compute `l` according to eq. 4-136 [3]
        # ----------------------------------------------------------------------
        l = @evalpoly(T∞, la[1], la[2], la[3], la[4], la[5])

        # Eq. 25' [1]
        # ----------------------------------------------------------------------
        γi = values(Mi).*(g₀*Ra²/( Rstar*l*T∞ )*(T∞ - Tx)/(Tx - T₁)*(zx - z₁)/(Ra + zx))

        # Eq. 25 [1]
        # ----------------------------------------------------------------------
        ρi = zeros(5)

        for i = 1:5
            ρi[i] = ρi_₁₂₅[i]*(Tx/Tz)^(1 + αi[i] + γi[i])* ( (T∞ - Tz)/
                                                             (T∞ - Tx) )^γi[i]
        end

        # Correction of seasonal variations of helium by latitude, eq. 4-101 [3]
        # ----------------------------------------------------------------------

        Δlog₁₀ρ_He = 0.65/(23.439291*π/180)*abs(δs)*
        ( sin( π/4 - glat*δs/( 2abs(δs) ) )^3 - 0.35355)

        ρi[_jr1971_id.He] *= 10^(Δlog₁₀ρ_He)

        # For altitude higher than 500 km, we must account for H.
        # ----------------------------------------------------------------------

        ρH = 0.0

        if h > 500
            # Compute the temperature and the H density at 500 km.
            T_₅₀₀       = _jr1971_T(500.0, Tx, T∞)
            log₁₀_T_₅₀₀ = log10(T_₅₀₀)
            ρH_₅₀₀      = Mi.H/Av*10^( 73.13 -
                                      ( 39.4 - 5.5log₁₀_T_₅₀₀ )*log₁₀_T_₅₀₀ )

            # Compute the H density at desired altitude.
            γH = Mi.H*g₀*Ra²/( Rstar*l*T∞ )*
                   (T∞ - Tx)/(Tx - T₁)*
                   (zx - z₁)/(Ra + zx)

            ρH = ρH_₅₀₀*(T_₅₀₀/Tz)^( 1 + αi.H + γH )*( (T∞ - Tz   )/
                                                       (T∞ - T_₅₀₀) )^γH
        end

        # Apply the corrections.
        ρi .*= 10^Δlog₁₀ρ_c
        ρH  *= 10^Δlog₁₀ρ_c

        # Convert to SI and return.
        return JR1971_Output(nN2   = ρi[1] * Av/Mi.N₂ * 1e6,
                             nO2   = ρi[2] * Av/Mi.O₂ * 1e6,
                             nO    = ρi[3] * Av/Mi.O  * 1e6,
                             nAr   = ρi[4] * Av/Mi.Ar * 1e6,
                             nHe   = ρi[5] * Av/Mi.He * 1e6,
                             nH    = ρH    * Av/Mi.H  * 1e6,
                             rho   = (sum(ρi) + ρH)*1e3,
                             T_exo = T∞,
                             Tz    = Tz )
    end
end
