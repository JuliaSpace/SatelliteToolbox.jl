# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Private functions and variables of the Jacchia-Roberts 1971 Atmospheric
#   Model.
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
#   [3] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [4] Long, A. C., Cappellari Jr., J. O., Velez, C. E., Fuchs, A. J (editors)
#       (1989). Goddard Trajectory Determination System (GTDS) Mathematical
#       Theory (Revision 1). FDD/552-89/0001 and CSC/TR-89/6001.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# NamedTuple that will be use when defining constants related to the species.
JR1971_NT{T} = NamedTuple{ (:N₂, :O₂, :O, :Ar, :He, :H), NTuple{6,T} }

"""
Structure with the constants for the Jacchia-Roberts 1971 Atmospheric Model.

"""
@with_kw struct JR1971_CONSTANTS{T}
    # General constants
    # ==========================================================================

    Rstar::T = T(8.31432)     # Rstar is the universal gas-constant [joules/(K.mol)].
    Av::T    = T(6.022045e23) # Avogadro's constant (mks) [molecules/mol].
    Ra::T    = T(6356.766)    # Mean Earth radius [km].
    g₀::T    = T(9.80665)     # Gravity at sea level [m/s²].

    # Values at the sea level
    # ==========================================================================

    M₀::T = 28.960            # Mean molecular mass at sea level [g/mol].

    # Values at lower boundary
    # ==========================================================================

    T₁::T = T(183)            # Temperature at the lower bound [K].
    z₁::T = T(90)             # Altitude of the lower bound [km].
    M₁::T = T(28.82678)       # Mean molecular mass at `z₁` [g/mol].
    ρ₁::T = T(3.46e-9)        # Density at `z₁` [g/cm³].

    # Values at the second division
    # ==========================================================================

    z₂::T = T(100)

    # Values at the inflection point
    # ==========================================================================

    zx::T = T(125)           # Altitude of the inflection point [km].

    # Molecular mass [g/mol]
    # ==========================================================================

    Mi::JR1971_NT{T} = ( N₂ = 28.0134,
                         O₂ = 31.9988,
                         O  = 15.9994,
                         Ar = 39.9480,
                         He =  4.0026,
                         H  = 1.00797 )

    # Thermal diffusion coefficient
    # ==========================================================================

    αi::JR1971_NT{T} = ( N₂ =  0,
                         O₂ =  0,
                         O  =  0,
                         Ar =  0,
                         He = -0.38,
                         H  =  0 )

    # Constituent density * M₀ / ρ_₁₀₀ / Av
    # ==========================================================================

    μi::JR1971_NT{T} = ( N₂ = 0.78110,
                         O₂ = 0.161778,
                         O  = 0.095544,
                         Ar = 0.0093432,
                         He = 0.61471e-5,
                         H  = 0 )

    # Coefficients for series expansion
    # ==========================================================================

    Aa::NTuple{7,T} = ( -435093.363387,
                          28275.5646391,
                           -765.33466108,
                             11.043387545,
                             -0.08958790995,
                              0.00038737586,
                             -0.000000697444 )

    Ca::NTuple{5,T} = ( -89284375.0,
                          3542400.0,
                           -52687.5,
                              340.5,
                               -0.8 )

    la::NTuple{5,T} = (  0.1031445e+5,
                         0.2341230e+1,
                         0.1579202e-2,
                        -0.1252487e-5,
                         0.2462708e-9 )

    α::NTuple{6,T} = ( 3144902516.672729,
                       -123774885.4832917,
                          1816141.096520398,
                           -11403.31079489267,
                               24.36498612105595,
                                0.008957502869707995 )

    β::NTuple{6,T} = ( -52864482.17910969,
                          -16632.50847336828,
                              -1.308252378125,
                               0,
                               0,
                               0)

    ζ::NTuple{7,T} = (  0.1985549e-10,
                       -0.1833490e-14,
                        0.1711735e-17,
                       -0.1021474e-20,
                        0.3727894e-24,
                       -0.7734110e-28,
                        0.7026942e-32 )

    δij::JR1971_NT{NTuple{7,T}} =
        (
          N₂ = (  0.1093155e2   ,
                  0.1186783e-2  ,    # (1/K)
                 -0.1677341e-5  ,    # (1/K)²
                  0.1420228e-8  ,    # (1/K)³
                 -0.7139785e-12 ,    # (1/K)⁴
                  0.1969715e-15 ,    # (1/K)⁵
                 -0.2296182e-19 ),   # (1/K)⁶

          O₂ = (  0.9924237e1   ,
                  0.1600311e-2  ,    # (1/K)
                 -0.2274761e-5  ,    # (1/K)²
                  0.1938454e-8  ,    # (1/K)³
                 -0.9782183e-12 ,    # (1/K)⁴
                  0.2698450e-15 ,    # (1/K)⁵
                 -0.3131808e-19 ),   # (1/K)⁶

          O =  (  0.1097083e2   ,
                  0.6118742e-4  ,    # (1/K)
                 -0.1165003e-6  ,    # (1/K)²
                  0.9239354e-10 ,    # (1/K)³
                 -0.3490739e-13 ,    # (1/K)⁴
                  0.5116298e-17 ,    # (1/K)⁵
                  0.0           ),   # (1/K)⁶

          Ar = (  0.8049405e1   ,
                  0.2382822e-2  ,    # (1/K)
                 -0.3391366e-5  ,    # (1/K)²
                  0.2909714e-8  ,    # (1/K)³
                 -0.1481702e-11 ,    # (1/K)⁴
                  0.4127600e-15 ,    # (1/K)⁵
                 -0.4837461e-19 ),   # (1/K)⁶

          He = (  0.7646886e1   ,
                 -0.4383486e-3  ,    # (1/K)
                  0.4694319e-6  ,    # (1/K)²
                 -0.2894886e-9  ,    # (1/K)³
                  0.9451989e-13 ,    # (1/K)⁴
                 -0.1270838e-16 ,    # (1/K)⁵
                  0.0           ),   # (1/K)⁶

          H =  (
                 0.0,
                 0.0,                # (1/K)
                 0.0,                # (1/K)²
                 0.0,                # (1/K)³
                 0.0,                # (1/K)⁴
                 0.0,                # (1/K)⁵
                 0.0 )               # (1/K)⁶
        )
end

"""
Constants for the Jacchia-Roberts 1971 Atmospheric Model.

"""
const _jr1971_constants = JR1971_CONSTANTS{Float64}()

"""
Index of the species for the Jacchia-Roberts 1971 Atmospheric Model.

"""
const _jr1971_id = ( N₂ = 1, O₂ = 2, O  = 3, Ar = 4, He = 5, H  = 6 )

################################################################################
#                              Private Functions
################################################################################

"""
    function _jr1971_M(z::R) where R

Compute the mean molecular mass at altitude `z` [km] using the empirical profile
in eq. 1 [3,4].

"""
@inline function _jr1971_M(z::R) where R
    !(90 <= z <= 100.0) && @warn "The empirical model for the mean molecular mass is valid only for 90 <= z <= 100 km."
    @unpack Aa = _jr1971_constants
    @evalpoly(z, Aa[1], Aa[2], Aa[3], Aa[4], Aa[5], Aa[6], Aa[7])
end

"""
    function _jr1971_roots(p::Polynomial{R}) where R

Compute the roots of the polynomial `p` necessary to compute the density below
125 km. It returns the value `r₁`, `r₂`, `x`, and `y`.

"""
function _jr1971_roots(p::Vector{R}) where R
    # Compute the roots with a first guess.
    r = roots(p, Complex{R}[166.10; 61.32; 9.91 + 1.31im; 9.91 - 1.31im];
              polish = true)

    # Take the real roots and sort them.
    r_real = sort( r[ @. abs(imag(r)) < 1e-10 ];
                   lt = (a,b) -> real(a) < real(b) )

    # There must be 2 real roots.
    (length(r_real) != 2) && error("Internal error, there are not two real roots!")

    r₁ = real(r_real[2])
    r₂ = real(r_real[1])

    # Take the complex roots and sort them.
    r_complex = sort( r[ @. abs(imag(r)) >= 1e-10 ];
                      lt = (a,b) -> imag(a) < imag(b) )

    # There must be 2 complex roots.
    (length(r_complex) != 2) && error("Internal error, there are not two complex roots!")

    x = real(r_complex[2])
    y = imag(r_complex[2])

    r₁, r₂, x, y
end

"""
    function _jr1971_T(z::R, Tx::R, T∞::R) where R<:Number

Compute the temperature [K] at height `z` [km] given the temperature `Tx` [K] at
the inflection point, and the exospheric temperature `T∞` [K] according to the
theory of the model Jacchia-Roberts 1971 [1,3,4].

The inflection point is considered to by `z = 125 km`.

"""
function _jr1971_T(z::R, Tx::R, T∞::R) where R<:Number
    @unpack T₁, z₁, zx = _jr1971_constants

    # Check the parameters
    # ==========================================================================

    (z  < z₁) && error("The altitude must not be lower than $(z₁) km.")
    (T∞ < 0 ) && error("The exospheric temperature must be positive.")

    # Compute the temperature at desire altitude
    # ==========================================================================

    if z <= zx
        @unpack Ca = _jr1971_constants

        aux = @evalpoly(z, Ca[1], Ca[2], Ca[3], Ca[4], Ca[5])
        T   = Tx + (Tx - T₁)/35^4*aux
    else
        @unpack Ra, la = _jr1971_constants

        l = @evalpoly(T∞, la[1], la[2], la[3], la[4], la[5])
        T = T∞ - ( T∞ - Tx )*exp( -( Tx - T₁ )/( T∞ - Tx )*
                                   (  z - zx )/( zx - z₁ )*
                                   (    l    )/( Ra + z  ) )
    end

    T
end
