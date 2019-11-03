# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   SGP4 orbit propagator model.
#
#   This is a independent implementation of the algorithm presented in [1].
#   Notice that the readability of the code was the major concern about the
#   implementation here. Algorithms with better performance can be found at
#   Vallado's repository in:
#
#       https://celestrak.com/software/vallado-sw.asp
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Hoots, F. R., Roehrich, R. L (1980). Models for Propagation of NORAD
#       Elements Set. Spacetrack Report No. 3.
#
#   [2] Vallado, D. A., Crawford, P., Hujsak, R., Kelso, T. S (2006). Revisiting
#       Spacetrack Report #3: Rev1. AIAA.
#
#   [3] SGP4 Source code of STRF: https://github.com/cbassa/strf
#       The SGP4 C code available on STRF was converted by Paul. S. Crawford and
#       Andrew R. Brooks.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export sgp4_gc_wgs72, sgp4_gc_wgs84
export sgp4_init, sgp4!
export dsinit, dsper!, dssec!

################################################################################
#                                  Overloads
################################################################################

# Copy for SGP4_Structure.
Base.copy(m::SGP4_Structure) = SGP4_Structure([ getfield(m, k) for k = 1:length(fieldnames(m)) ]...)

# Deepcopy for SGP4_Structure.
Base.deepcopy(m::SGP4_Structure) = SGP4_Structure([ deepcopy(getfield(m, k)) for k = 1:length(fieldnames(m)) ]...)

################################################################################
#                                  Constants
################################################################################

# WGS-84 / EGM-08 Gravitational constants.
const sgp4_gc_wgs84 = SGP4_GravCte{Float64}(
        R0/1000,
        60.0/sqrt(6378.137^3/398600.5),
         0.00108262998905,
        -0.00000253215306,
        -0.00000161098761
       )

# WGS-72 Gravitational constants.
const sgp4_gc_wgs72 = SGP4_GravCte{Float64}(
        6378.135,
        60.0/sqrt(6378.135^3/398600.8),
         0.001082616,
        -0.00000253881,
        -0.00000165597
       )

# Constants realted to the Deep Space perturbations.
# ==================================================

const STEP          = 720.0
const MAX_INTEGRATE = STEP * 10000
const ZNS           = 1.19459E-5
const C1SS          = 2.9864797e-6
const ZES           = 0.01675
const ZNL           = 1.5835218e-4
const ZEL           = 0.05490
const C1L           = 4.7968065e-7
const ZSINIS        = 0.39785416
const ZCOSIS        = 0.91744867
const ZCOSGS        = 0.1945905
const ZSINGS        = -0.98088458
const Q22           = 1.7891679e-6
const Q31           = 2.1460748e-6
const Q33           = 2.2123015e-7
const G22           = 5.7686396
const G32           = 0.95240898
const G44           = 1.8014998
const G52           = 1.0508330
const G54           = 4.4108898
const ROOT22        = 1.7891679e-6
const ROOT32        = 3.7393792e-7
const ROOT44        = 7.3636953e-9
const ROOT52        = 1.1428639e-7
const ROOT54        = 2.1765803e-9
const THDT          = 4.37526908801129966e-3

################################################################################
#                                  Functions
################################################################################

"""
    function sgp4_init(spg4_gc::SGP4_GravCte{T}, epoch::Number, n_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, M_0::Number, bstar::Number) where T
    function sgp4_init(tle::TLE, sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84)

Initialize the data structure of SGP4 algorithm.

# Args

* `spg4_gc`: SPG4 gravitational constants (see `SGP4_GravCte`).
* `epoch`: Epoch of the orbital elements [Julian Day].
* `n_0`: SGP type "mean" mean motion at epoch [rad/min].
* `e_0`: "Mean" eccentricity at epoch.
* `i_0`: "Mean" inclination at epoch [rad].
* `Ω_0`: "Mean" longitude of the ascending node at epoch [rad].
* `ω_0`: "Mean" argument of perigee at epoch [rad].
* `M_0`: "Mean" mean anomaly at epoch [rad].
* `bstar`: Drag parameter (B*).
* `tle`: TLE to initialize the SPG4 (see `TLE`).

# Returns

The structure `SGP4_Structure` with the initialized parameters.

"""
function sgp4_init(tle::TLE, sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84)
    return sgp4_init(sgp4_gc, tle.epoch, tle.n*2*pi/(24*60), tle.e,
                     tle.i*pi/180, tle.Ω*pi/180, tle.ω*pi/180, tle.M*pi/180,
                     tle.bstar)
end

function sgp4_init(sgp4_gc::SGP4_GravCte{T},
                   epoch::Number,
                   n_0::Number,
                   e_0::Number,
                   i_0::Number,
                   Ω_0::Number,
                   ω_0::Number,
                   M_0::Number,
                   bstar::Number) where T

    # Unpack the gravitational constants to improve code readability.
    @unpack_SGP4_GravCte sgp4_gc

    # Constants
    # =========
    #
    # Note: [er] = Earth radii.

    # Distance units / Earth radii.
    AE = T(1)

    k_2  = +1/2*J2*AE^2
    k_4  = -3/8*J4*AE^4
    A_30 = -J3*AE^3

    # Kilometers / Earth radii.
    XKMPER = R0

    # Parameters for the SGP4 density function.
    s   =  78/XKMPER + 1
    q_0 = 120/XKMPER + 1

    # (q_0-s)^4 [er]^4
    QOMS2T = (q_0-s)^4

    # ==========================================================================

    # Auxiliary variables to improve the performance.
    # ===============================================

    e_0² = e_0^2

    sin_i_0, θ = sincos(i_0)
    θ²         = θ^2
    θ³         = θ^3
    θ⁴         = θ^4

    # ==========================================================================

    # Recover the original mean motion (nll_0) and semi-major axis (all_0) from
    # the input elements.

    aux = (3*θ²-1)/(1-e_0²)^(3/2)

    a_1 = (XKE/n_0)^(2/3)
    δ_1 = 3/2*k_2/a_1^2*aux
    a_0 = a_1*(1 - 1/3*δ_1 - δ_1^2 - 134/81*δ_1^3)
    δ_0 = 3/2*k_2/a_0^2*aux

    nll_0 = n_0/(1 + δ_0)

    # Vallado's implementation of SGP4 [2] compute the semi-major axis
    # considering the new angular velocity, which is called `no_unkozai`. In the
    # original SGP4 technical report [1], the semi-major axis was computed
    # considering:
    #
    #   all_0 = a_0/(1 - δ_0)
    #
    all_0 = (XKE/nll_0)^(2/3)

    # Initialization
    # ==============

    # Compute the orbit perigee [ER].
    perigee = (all_0*(1-e_0)-AE)*XKMPER

    # For perigee below 156 km, the values of S and QOMS2T are altered.
    if perigee < 156
        if perigee < 98
            s = 20/XKMPER + AE
        # Perigee between 98km and 156km.
        else
            s = all_0*(1-e_0) - s + AE
        end

        QOMS2T = (q_0 - s)^4
    end

    # Compute SGP4 constants.
    ξ    = 1/(all_0-s)
    β_0  = sqrt(1-e_0^2)
    η    = all_0*e_0*ξ

    # Vallado's implementation of SGP4 [2] considers the absolute value of
    # (1-η^2) here and in the C2 and C4 computation. Notice that, if (1-η^2) <
    # 0, then aux1 cannot be computed. The original SGP4 technical report [1]
    # does not mention anything about this.

    aux0 = abs(1-η^2)
    aux1 = aux0^(-7/2)
    aux2 = ξ^4*all_0*β_0^2*aux1

    C2 = QOMS2T*ξ^4*nll_0*aux1*
         ( all_0*( 1 + (3/2)η^2 + 4e_0*η + e_0*η^3) +
           3/2*(k_2*ξ)/aux0*(-1/2 + (3/2)θ²)*(8 + 24η^2 + 3η^4) )

    C1 = bstar*C2

    C3 = (e_0 > 1e-4) ? QOMS2T*ξ^5*A_30*nll_0*AE*sin_i_0/(k_2*e_0) : T(0)

    C4 = 2nll_0*QOMS2T*aux2*
         ( 2η*(1+e_0*η) + (1/2)e_0 + (1/2)η^3 -
           2k_2*ξ/(all_0*aux0)*
                (3(1-3θ²)*(1 + (3/2)η^2 - 2e_0*η - (1/2)e_0*η^3) +
                3/4*(1-θ²)*(2η^2 - e_0*η - e_0*η^3)*cos(2ω_0)))

    C5 = 2QOMS2T*aux2*(1 + (11/4)η*(η+e_0) + e_0*η^3)

    D2 = 4all_0*ξ*C1^2

    D3 = 4/3*all_0*ξ^2*( 17all_0 +   s)*C1^3

    # Vallado's implementation of SGP4 [2] uses all_0^2, instead of only all_0
    # that is seen in the original SGP4 Technical Report [1].
    D4 = 2/3*all_0^2*ξ^3*(221all_0 + 31s)*C1^4

    # Compute the time-derivative of some orbital elements.
    dotM = ( 1 + 3k_2*(-1 + 3θ²)/(2all_0^2*β_0^3) +
             3k_2^2*(13 - 78θ² + 137θ⁴)/(16all_0^4*β_0^7) )*nll_0

    dotω = ( -3k_2*(1-5θ²)/(2all_0^2*β_0^4) +
              3k_2^2*(7 - 114θ² + 395θ⁴)/(16all_0^4*β_0^8) +
              5k_4*(3 - 36θ² + 49θ⁴)/(4all_0^4*β_0^8) )*nll_0

    dotΩ1 = -3k_2*θ/(all_0^2*β_0^4)*nll_0

    dotΩ  = dotΩ1 + ( 3k_2^2*(4θ - 19θ³)/(2all_0^4*β_0^8) +
                        5k_4*θ*(3 - 7θ²)/(2all_0^4*β_0^8) )*nll_0

    # The current orbital parameters are obtained from the TLE.
    a_k = all_0
    e_k = e_0
    i_k = i_0
    Ω_k = Ω_0
    ω_k = ω_0
    M_k = M_0
    n_k = nll_0

    sgp4_ds::SGP4_DeepSpace{T} = SGP4_DeepSpace{T}()

    # If the orbit period is higher than 225 min., then we must consider the
    # deep space perturbations. This is indicated by selecting the algorithm
    # `:sdp4`.
    if 2π / n_0 >= 225.0
        algorithm = :sdp4

        # Initialize the values for the SDP4 (deep space) algorithm.
        sgp4_ds = dsinit(epoch,
                         nll_0, all_0, e_0, i_0, Ω_0, ω_0, M_0,
                         dotM, dotω, dotΩ)
    else
        # For perigee lower than 220 km, the equations are truncated to a linear
        # variation in `sqrt(a)` and quadratic variation in mean anomaly. Also,
        # the C5 term, the δω term, and the δM term are dropped. This is
        # indicated by selecting the algorithm `:sgp4_lowper`. Otherwise, if
        # perigee is higher or equal 220 km and the orbit period is lower than
        # 225 min., then we use the normal SGP4 algorithm by selecting `:sgp4`.
        algorithm = (perigee/AE >= (220 + (AE-1)*XKMPER)) ? :sgp4 : :sgp4_lowper
    end

    # Create the output structure with the data.
    SGP4_Structure{T}(

        epoch, n_0, e_0, i_0, Ω_0, ω_0, M_0, bstar, 0, a_k, e_k, i_k, Ω_k, ω_k,
        M_k, n_k, all_0, nll_0, AE, QOMS2T, β_0, ξ, η, sin_i_0, θ, θ², A_30,
        k_2, k_4, C1, C3, C4, C5, D2, D3, D4, dotM, dotω, dotΩ1, dotΩ,
        algorithm, sgp4_gc, sgp4_ds

       )
end

"""
    function sgp4!(sgp4d::SGP4_Structure{T}, t::Number) where T

Propagate the orbit defined in `sgp4d` (see `SGP4_Structure`) until the time `t`
[min]. Notice that the values in `sgp4d` will be modified.

# Returns

* The position vector represented in TEME frame at time `t` [km].
* The velocity vector represented in TEME frame at time `t` [km/s].

"""
function sgp4!(sgp4d::SGP4_Structure{T}, t::Number) where T
    # Unpack variables.
    @unpack_SGP4_Structure sgp4d
    @unpack_SGP4_GravCte   sgp4_gc

    # After unpacking sgp4d, we have two sets of orbit elements:
    #
    #   (n_0, e_0, i_0, Ω_0, ω_0, M_0),
    #
    # and
    #
    #   (n_k, e_k, i_k, Ω_k, ω_k, M_k).
    #
    # The first are those initial elements from the orbit defined in `sgp4_init`
    # function. The second are the current elements. During this functions, the
    # second set is updated by adding the many effects considered in SGP4.

    # Time elapsed since epoch.
    #
    # We convert to `T` to avoid numerical problems with very big numbers as
    # pointed out in:
    #
    #   https://github.com/JuliaLang/julia/issues/27355
    Δt = T(t)

    # Initialization of the current elements with the values of the epoch.
    n_k = nll_0
    a_k = all_0
    e_k = e_0
    i_k = i_0
    Ω_k = Ω_0
    ω_k = ω_0
    M_k = M_0

    # Auxiliary variables to improve code performance.
    sin_i_k = sin_i_0

    # Secular effects of atmospheric drag and gravitation.
    # ====================================================

    M_k = M_0 + dotM*Δt
    Ω_k = Ω_0 + dotΩ*Δt - 21/2*(nll_0*k_2*θ)/(all_0^2*β_0^2)*C1*Δt^2
    ω_k = ω_0 + dotω*Δt

    # Check if we need to use SDP4 (deep space) algorithm.
    if algorithm == :sdp4
        # Compute the elements perturbed by the secular effects.
        (n_k, e_k, i_k, Ω_k, ω_k, M_k) = dssec!(sgp4_ds,
                                                nll_0,
                                                e_0,
                                                i_0,
                                                ω_0,
                                                Ω_k,
                                                ω_k,
                                                M_k,
                                                dotω,
                                                Δt)

        a_k  = (XKE / n_k)^(2/3) * ( 1 - C1*Δt )^2
        e_k += -bstar*C4*Δt
        M_k += nll_0*(1.5C1*Δt^2)

    # Check if perigee is above 220 km.
    elseif algorithm == :sgp4

        sin_M_0, cos_M_0 = sincos(M_0)

        δω  = bstar*C3*cos(ω_0)*Δt

        δM  = (e_0 > 1e-4) ?
            -2/3*QOMS2T*bstar*ξ^4*AE/(e_0*η)*( (1 + η*cos(M_k))^3 -
                                               (1 + η*cos_M_0 )^3 ) : T(0)

        M_k += +δω + δM

        ω_k += -δω - δM

        e_k = e_0 - bstar*C4*Δt - bstar*C5*(sin(M_k) - sin_M_0)

        a_k = all_0*( @evalpoly(Δt, 1, -C1, -D2, -D3, -D4) )^2

        IL  = M_k + ω_k + Ω_k + nll_0*
              ( @evalpoly(Δt, 0, 0, 3/2*C1, +(D2 + 2C1^2),
                          +(1/4)*(3D3 + 12C1*D2 + 10C1^3),
                          +(1/5)*(3D4 + 12C1*D3 + 6D2^2 + 30C1^2*D2 + 15C1^4) ) )

    elseif algorithm == :sgp4_lowper
        # If so, then
        #     1. Drop all terms after C1 in `a` and `IL`.
        #     2. Drop all terms involving C5.
        #     3. Drop δω.
        #     4. Drop δM.
        e_k = e_0 - bstar*C4*Δt

        a_k = all_0*(1 - C1*Δt)^2

        IL  = M_k + ω_k + Ω_k + nll_0*(3/2)*C1*Δt^2
    else
        error("Unknown algorithm :$algorithm. Possible values are :sgp4, :sgp4_lowper, :sdp4.")
    end

    # TODO: Vallado's implementation [2] apply this normalization to the mean
    # anomaly. It is necessary to verify the reason for that.
    M_k_aux = M_k + ω_k + Ω_k
    Ω_k     = rem(Ω_k, 2π)
    ω_k     = rem(ω_k, 2π)
    M_k_aux = rem(M_k_aux, 2π)
    M_k     = rem(M_k_aux - ω_k - Ω_k, 2π)

    # Lunar-Solar Periodics for Deep Space Orbits
    # ===========================================

    # This is only necessary if we are using SDP4 algorithm.
    if algorithm == :sdp4
        # Compute the elements perturbed by the Lunar-Solar periodics.
        (e_k, i_k, Ω_k, ω_k, M_k) = dsper!(sgp4_ds,
                                           e_k,
                                           i_k,
                                           Ω_k,
                                           ω_k,
                                           M_k,
                                           Δt)

        IL = M_k + ω_k + Ω_k

        # Make sure that the inclination is always positive.
        if i_k < 0
            i_k = -i_k
            Ω_k += π
            ω_k -= π
        end

        # The inclination was changed, hence some auxiliary variables must be
        # recomputed.
        sin_i_k, θ = sincos(i_k)
        θ²         = θ^2
    end

    # Vallado's code does not let the eccentricity to be smaller than 1e-6.
    #
    # TODO: Verify why this is necessary. I did not find any reason for that.
    (e_k < 1e-6) && (e_k = T(1e-6))

    β = sqrt(1-e_k^2)

    # Compute the angular velocity [rad/min].
    n_k = XKE/a_k^(3/2)

    # Long-period periodic terms.
    # ===========================

    sin_ω_k, cos_ω_k = sincos(ω_k)

    a_xN = e_k*cos_ω_k

    # TODO: Vallado's implementation of SGP4 uses another equation here.
    # However, both produces the same result. Verify which one is better.
    #
    a_yNL = A_30*sin_i_k/(4*k_2*a_k*β^2)
    a_yN = e_k*sin_ω_k + a_yNL

    IL_L =  (1/2)a_yNL*a_xN*(3 + 5θ)/(1 + θ)

    IL_T = IL + IL_L

    # Solve Kepler's equation for (E + ω).
    # ====================================

    U = rem(IL_T - Ω_k, 2π)

    E_ω = U

    # Define the following variables that will be modified inside the loop so
    # that we can use them after the loop.
    sin_E_ω = T(0)
    cos_E_ω = T(0)

    for k = 1:10
        sin_E_ω, cos_E_ω = sincos(E_ω)

        ΔE_ω = (U - a_yN*cos_E_ω + a_xN*sin_E_ω - E_ω)/
               (1 - a_yN*sin_E_ω - a_xN*cos_E_ω)

        # Vallado proposes to limit the maximum increment.
        (abs(ΔE_ω) >= 0.95) && (ΔE_ω = sign(ΔE_ω)*0.95)

        E_ω += ΔE_ω

        # If the increment is less than a threshold, break the loop.
        #
        # Vallado proposes a threshold of 10^-12 instead of 10^-6.
        (abs(ΔE_ω) < 1e-12) && (break)
    end

    # Short-term periodic terms.
    # ==========================

    # Auxiliary variables.
    #
    # Note: the sine and cosine of E+ω was already computed in the previous
    # loop.

    e_cos_E = a_xN*cos_E_ω + a_yN*sin_E_ω
    e_sin_E = a_xN*sin_E_ω - a_yN*cos_E_ω

    e_L     = sqrt(a_xN^2 + a_yN^2)

    p_L     = a_k*(1-e_L^2)

    r       = a_k*(1-e_cos_E)

    dot_r   = XKE*sqrt(a_k)*e_sin_E/r

    r_dot_f = XKE*sqrt(p_L)/r

    aux     = e_sin_E/(1 + sqrt(1 - e_L^2))
    cos_u   = a_k/r*( cos_E_ω - a_xN + a_yN*aux )
    sin_u   = a_k/r*( sin_E_ω - a_yN - a_xN*aux )
    cos_2u  = 1 - 2*sin_u^2
    sin_2u  = 2cos_u*sin_u
    u       = atan(sin_u, cos_u)

    # Short-term periodic terms.

    Δr       = +k_2/(2p_L)*(1 - θ²)*cos_2u

    Δu       = -k_2/(4p_L^2)*(7θ² - 1)*sin_2u

    ΔΩ       = +3k_2*θ/(2p_L^2)*sin_2u

    Δi       = +3k_2*θ/(2p_L^2)*sin_i_k*cos_2u

    Δdot_r   = -k_2*n_k/p_L*(1 - θ²)*sin_2u

    Δr_dot_f = +k_2*n_k/p_L*( (1 - θ²)*cos_2u - 3/2*(1 - 3θ²) )

    # The short-term periodics are added to give the osculating quantities.

    r_k       = r*( 1 - (3/2)k_2*sqrt(1 - e_L^2)/p_L^2*(3θ² - 1) ) + Δr

    u_k       = u + Δu

    Ω_k       = Ω_k + ΔΩ

    i_k       = i_k + Δi

    dot_r_k   = dot_r + Δdot_r

    r_dot_f_k = r_dot_f + Δr_dot_f

    # Orientation vectors.
    sin_Ω_k, cos_Ω_k = sincos(Ω_k)
    sin_i_k, cos_i_k = sincos(i_k)
    sin_u_k, cos_u_k = sincos(u_k)

    M = SVector{3}(-sin_Ω_k*cos_i_k,
                   +cos_Ω_k*cos_i_k,
                            sin_i_k )

    N = SVector{3}(+cos_Ω_k,
                   +sin_Ω_k,
                   +T(0)   )

    Uv = M*sin_u_k + N*cos_u_k
    Vv = M*cos_u_k - N*sin_u_k

    r_TEME = r_k*Uv*R0
    v_TEME = (dot_r_k*Uv + r_dot_f_k*Vv)*R0/60

    # Update the variables.
    sgp4d.Δt  = Δt
    sgp4d.a_k = a_k
    sgp4d.e_k = e_k
    sgp4d.i_k = i_k
    sgp4d.Ω_k = Ω_k
    sgp4d.ω_k = ω_k
    sgp4d.M_k = M_k
    sgp4d.n_k = n_k

    (r_TEME, v_TEME)
end

################################################################################
#                             Deep Space Functions
################################################################################

"""
    function dsinit(epoch::T, nll_0::T, all_0::T, e_0::T, i_0::T, Ω_0::T, ω_0::T, M_0::T, dotM::T, dotω::T, dotΩ::T) where T<:Number

Initialize the deep space structure. This function performs the initial
computations and save the values at an instance of the structure
`SGP4_DeepSpace`. Those will be used when calling the functions `dsper!` and
`dpsec!`.

# Args

* `epoch`: Epoch of the initial orbit [Julian Day].
* `nll_0`: Initial mean motion [rad/min].
* `all_0`: Initial semi-major axis [ER].
* `e_0`: Initial eccentricity.
* `i_0`: Initial inclination [rad].
* `Ω_0`: Initial right ascencion of the ascending node [rad].
* `ω_0`: Initial argument of perigee [rad].
* `M_0`: Initial mean motion [rad].
* `dotM`: Time-derivative of the mean motion [rad/min].
* `dotω`: Time-derivative of the argument of perigee [rad/min].
* `dotΩ`: Time-derivative of the RAAN [rad/min].

# Returns

An instance of the structure `SGP4_DeepSpace` with the initalized values.

"""
function dsinit(epoch::T,
                nll_0::T,
                all_0::T,
                e_0::T,
                i_0::T,
                Ω_0::T,
                ω_0::T,
                M_0::T,
                dotM::T,
                dotω::T,
                dotΩ::T) where T<:Number

    sgp4_ds::SGP4_DeepSpace{T} = SGP4_DeepSpace{T}()
    @unpack_SGP4_DeepSpace sgp4_ds

    #                         Auxiliary Variables
    # ==========================================================================
    e_0²        = e_0^2
    e_0³        = e_0 * e_0²
    sqrt_1_e_0² = sqrt((1 - e_0²))
    inv_all_0   = 1/all_0
    inv_nll_0   = 1/nll_0
    se          = T(0)
    si          = T(0)
    sl          = T(0)
    sgh         = T(0)
    shdq        = T(0)

    sin_i_0, cos_i_0 = sincos(i_0)
    sin_Ω_0, cos_Ω_0 = sincos(Ω_0)
    sin_ω_0, cos_ω_0 = sincos(ω_0)

    sin_i_0²    = sin_i_0^2
    cos_i_0²    = cos_i_0^2
    xpidot      = dotω + dotΩ

    #                        Initial Configuration
    # ==========================================================================

    # Drop terms if inclination is smaller than 3 deg.
    ishq = (i_0 >= 3π/180) ? true : false

    # Do not let `sin_i_0` be 0.
    (abs(sin_i_0) < 1e-12) && (sin_i_0 = sign(sin_i_0)*1e-12)

    # Compute the Greenwhich Mean Sidereal Time at epoch.
    gmst = JDtoGMST(epoch)

    #                      Initialize Lunar Solar Terms
    # ==========================================================================

    # `day` is the number of days since Jan 0, 1900 at 12h.
    day    = epoch - (datetime2julian(DateTime(1900,1,1,12,0,0))-1)

    xnodce     = mod(4.5236020 - 9.2422029e-4day, 2π)
    stem, ctem = sincos(xnodce)

    zcosil = 0.91375164 - 0.03568096ctem
    zsinil = sqrt(1 - zcosil^2)

    zsinhl = 0.089683511stem/zsinil
    zcoshl = sqrt(1 - zsinhl^2)

    gam    = 5.8351514 + 0.0019443680day

    zx     = 0.39785416stem/zsinil
    zy     = zcoshl*ctem + 0.91744867zsinhl*stem
    zx     = atan(zx,zy)
    zx     = gam + zx - xnodce

    zsingl, zcosgl = sincos(zx)

    zmol   = mod(4.7199672 + 0.22997150day - gam, 2π)
    zmos   = mod(6.2565837 + 0.017201977day,      2π)

    #                            Do Solar Terms
    # ==========================================================================
    zcosg = ZCOSGS
    zsing = ZSINGS
    zcosi = ZCOSIS
    zsini = ZSINIS
    zcosh = cos_Ω_0
    zsinh = sin_Ω_0
    cc    = C1SS
    zn    = ZNS
    ze    = ZES
    zmo   = zmos

    for ls = 0:1
        a1  = +zcosg*zcosh + zsing*zcosi*zsinh
        a3  = -zsing*zcosh + zcosg*zcosi*zsinh
        a7  = -zcosg*zsinh + zsing*zcosi*zcosh
        a8  = +zsing*zsini
        a9  = +zsing*zsinh + zcosg*zcosi*zcosh
        a10 = +zcosg*zsini
        a2  = +cos_i_0*a7  + sin_i_0*a8
        a4  = +cos_i_0*a9  + sin_i_0*a10
        a5  = -sin_i_0*a7  + cos_i_0*a8
        a6  = -sin_i_0*a9  + cos_i_0*a10

        x1 = +a1*cos_ω_0 + a2*sin_ω_0
        x2 = +a3*cos_ω_0 + a4*sin_ω_0
        x3 = -a1*sin_ω_0 + a2*cos_ω_0
        x4 = -a3*sin_ω_0 + a4*cos_ω_0
        x5 = +a5*sin_ω_0
        x6 = +a6*sin_ω_0
        x7 = +a5*cos_ω_0
        x8 = +a6*cos_ω_0

        z31 = 12x1^2    - 3x3^2
        z32 = 24x1 * x2 - 6x3 * x4
        z33 = 12x2^2    - 3x4^2
        z1  = 3(   a1^2 + a2^2   ) + z31 * e_0²
        z2  = 6(a1 * a3 + a2 * a4) + z32 * e_0²
        z3  = 3(   a3^2 + a4^2   ) + z33 * e_0²
        z11 = -6a1 * a5 + e_0² * (-24x1 * x7 - 6x3 * x5)
        z12 = -6(a1 * a6 + a3 * a5) + e_0² * ( -24(x2 * x7 + x1 * x8) -
                                                 6(x3 * x6 + x4 * x5) )
        z13 = -6a3 * a6 + e_0² * (-24x2 * x8 - 6x4 * x6)
        z21 = +6a2 * a5 + e_0² * (+24x1 * x5 - 6x3 * x7)
        z22 = +6(a4 * a5 + a2 * a6) + e_0² * (24(x2 * x5 + x1 * x6) -
                                               6(x4 * x7 + x3 * x8) )
        z23 = +6a4 * a6 + e_0² * (24x2 * x6 - 6x4 * x8)
        z1  = +2z1 + (1 - e_0²) * z31
        z2  = +2z2 + (1 - e_0²) * z32
        z3  = +2z3 + (1 - e_0²) * z33
        s3  = +cc * inv_nll_0
        s2  = -0.5s3 / sqrt_1_e_0²
        s4  = +s3 * sqrt_1_e_0²
        s1  = -15e_0 * s4
        s5  = +x1 * x3 + x2 * x4
        s6  = +x2 * x3 + x1 * x4
        s7  = +x2 * x4 - x1 * x3
        se  = +s1 * zn * s5
        si  = +s2 * zn * (z11 + z13)
        sl  = -zn * s3 * (z1 + z3 - 14 - 6e_0²)
        sgh = +s4 * zn * (z31 + z33 - 6)

        shdq = zero(T)

        if ishq
            sh   = -zn * s2 * (z21 + z23);
			shdq = sh / sin_i_0;
        end

        ee2  =  +2s1 * s6
        e3   =  +2s1 * s7
        xi2  =  +2s2 * z12
        xi3  =  +2s2 * (z13 - z11)
        xl2  =  -2s3 * z2
        xl3  =  -2s3 * (z3 - z1)
        xl4  =  -2s3 * (-21 - 9e_0²) * ze
        xgh2 =  +2s4 * z32
        xgh3 =  +2s4 * (z33 - z31)
        xgh4 = -18s4 * ze
        xh2  =  -2s2 * z22
        xh3  =  -2s2 * (z23 - z21)

        (ls == 1) && break

        #                        Do Lunar Terms
        # ======================================================================
        sse   = se
        ssi   = si
        ssl   = sl
        ssh   = shdq
        ssg   = sgh - cos_i_0 * ssh
        se2   = ee2
        si2   = xi2
        sl2   = xl2
        sgh2  = xgh2
        sh2   = xh2
        se3   = e3
        si3   = xi3
        sl3   = xl3
        sgh3  = xgh3
        sh3   = xh3
        sl4   = xl4
        sgh4  = xgh4
        zcosg = zcosgl
        zsing = zsingl
        zcosi = zcosil
        zsini = zsinil
        zcosh = cos_Ω_0 * zcoshl + sin_Ω_0 * zsinhl
        zsinh = sin_Ω_0 * zcoshl - cos_Ω_0 * zsinhl
        zn    = ZNL
        cc    = C1L
        ze    = ZEL
        zmo   = zmol
    end

    sse += se
    ssi += si
    ssl += sl
    ssg += sgh - cos_i_0 * shdq
    ssh += shdq

    if (nll_0 < 0.0052359877) && (nll_0 > 0.0034906585)
        #        24h Synchronous Resonance Terms Initialization
        # ======================================================================

        iresfl = true;
        isynfl = true;

        g200    = e_0² * (0.8125e_0² - 2.5) + 1
        g310    = 2e_0² + 1
        g300    = e_0² * (6.60937e_0² - 6) + 1
        f220    = 0.75(cos_i_0 + 1)^2
        f311    = 0.9375(3.0cos_i_0 + 1)*sin_i_0^2 - 0.75(cos_i_0 + 1)
        f330    = 1.875 * (cos_i_0 + 1)^3
        del1    = 3(nll_0^2 * inv_all_0^2)
        del2    = 2del1 * f220 * g200 * Q22
        del3    = 3del1 * f330 * g300 * Q33 * inv_all_0
        del1    =  del1 * f311 * g310 * Q31 * inv_all_0
        fasx2   = 0.13130908
        fasx4   = 2.8843198
        fasx6   = 0.37448087
        xlamo   = mod(M_0 + Ω_0 + ω_0 - gmst, 2pi)
        bfact   = dotM + xpidot - THDT + ssl + ssg + ssh

    elseif (nll_0 >= 0.00826) && (nll_0 <= 0.00924) && (e_0 >= 0.5)
        #   Geopotential Resonance Initialization for 12 Hour Orbits
        # ======================================================================

        iresfl = true
        isynfl = false

        g201 = -0.306 - 0.44(e_0 - 0.64)

        if e_0 <= 0.65
            g211 = @evalpoly(e_0, + 3.6160, - 13.2470, + 16.29000)
            g310 = @evalpoly(e_0, -19.3020, +117.3900, -228.4190, +156.5910)
            g322 = @evalpoly(e_0, -18.9068, +109.7927, -214.6334, +146.5816)
            g410 = @evalpoly(e_0, -41.1220, +242.6940, -471.0940, +313.9530)
            g422 = @evalpoly(e_0, -146.407, +841.8800, -1629.014, +1083.435)
            g520 = @evalpoly(e_0, -532.114, +3017.977, -5740.032, +3708.276)
        else
            g211 = @evalpoly(e_0, -  72.099, +331.8190, -508.7380, +266.7240)
            g310 = @evalpoly(e_0, - 346.844, +1582.851, -2415.925, +1246.113)
            g322 = @evalpoly(e_0, - 342.585, +1554.908, -2366.899, +1215.972)
            g410 = @evalpoly(e_0, -1052.797, +4758.686, -7193.992, +3651.957)
            g422 = @evalpoly(e_0, -3581.690, +16178.11, -24462.77, +12422.52)

            if e_0 <= 0.715
                g520 = @evalpoly(e_0, +1464.74, -4664.75, +3763.64)
            else
                g520 = @evalpoly(e_0, -5149.66, +29936.92, -54087.36, +31324.56)
            end
        end

        if e_0 < 0.7
            g533 = @evalpoly(e_0, -919.22770, +4988.6100, -9064.7700, +5542.210)
            g521 = @evalpoly(e_0, -822.71072, +4568.6173, -8491.4146, +5337.524)
            g532 = @evalpoly(e_0, -853.66600, +4690.2500, -8624.7700, +5341.400)
        else
            g533 = @evalpoly(e_0, -37995.780, +161616.52, -229838.20, +109377.94)
            g521 = @evalpoly(e_0, -51752.104, +218913.95, -309468.16, +146349.42)
            g532 = @evalpoly(e_0, -40023.880, +170470.89, -242699.48, +115605.82)
        end

        f220 = +0.75(1 + 2cos_i_0 + cos_i_0²)
        f221 = +1.5sin_i_0²
        f321 = +1.875sin_i_0 * (1 - 2cos_i_0 - 3cos_i_0²)
        f322 = -1.875sin_i_0 * (1 + 2cos_i_0 - 3cos_i_0²)
        f441 = +35sin_i_0² * f220
        f442 = +39.375sin_i_0²^2
        f522 = +9.84375sin_i_0 * (sin_i_0² * (+1 - 2cos_i_0 - 5cos_i_0²) +
                                   0.33333333(-2 + 4cos_i_0 + 6cos_i_0²) )
        f523 = sin_i_0 * (4.92187512sin_i_0² * (-2 - 4cos_i_0 + 10cos_i_0²) +
                                     6.56250012(+1 + 2cos_i_0 -  3cos_i_0²) )
        f542 = 29.53125sin_i_0 * (+2 - 8cos_i_0 + cos_i_0² *
                                    (-12 + 8cos_i_0 + 10cos_i_0²))
        f543 = 29.53125sin_i_0 * (-2 - 8cos_i_0 + cos_i_0² *
                                    (+12 + 8cos_i_0 - 10cos_i_0²))

        temp1   = 3(nll_0 * inv_all_0)^2
        temp0   = temp1 * ROOT22
        d2201   = temp0 * f220 * g201
        d2211   = temp0 * f221 * g211
        temp1  *= inv_all_0
        temp0   = temp1 * ROOT32
        d3210   = temp0 * f321 * g310
        d3222   = temp0 * f322 * g322
        temp1  *= inv_all_0
        temp0   = 2temp1 * ROOT44
        d4410   = temp0 * f441 * g410
        d4422   = temp0 * f442 * g422
        temp1  *= inv_all_0
        temp0   = temp1 * ROOT52
        d5220   = temp0 * f522 * g520
        d5232   = temp0 * f523 * g532
        temp0   = 2temp1 * ROOT54
        d5421   = temp0 * f542 * g521
        d5433   = temp0 * f543 * g533
        xlamo   = mod(M_0 + 2Ω_0 - 2gmst, 2pi)
        bfact   = dotM + 2dotΩ - 2THDT + ssl + 2ssh
    else
        #                     Non Resonant Orbits
        # ======================================================================

        iresfl = false
        isynfl = false
    end

	if iresfl
        #                   Initialize the Integrator
        # ======================================================================
		xfact = bfact - nll_0
		xli   = xlamo
		atime = 0.0

        # TODO: Check if this variable can be removed from SGP4_DeepSpace.
        xni   = nll_0

        # Compute the "dot" terms.
        # ========================

        if (isynfl)
            sin_1, cos_1 = sincos(  (xli - fasx2) )
            sin_2, cos_2 = sincos( 2(xli - fasx4) )
            sin_3, cos_3 = sincos( 3(xli - fasx6) )

            xndot = del1 * sin_1 +  del2 * sin_2 +  del3 * sin_3
            xnddt = del1 * cos_1 + 2del2 * cos_2 + 3del3 * cos_3
        else
            ω     = ω_0 + dotω * atime

            sin_1,  cos_1  = sincos(2ω + xli  - G22)
            sin_2,  cos_2  = sincos(   + xli  - G22)
            sin_3,  cos_3  = sincos(+ω + xli  - G32)
            sin_4,  cos_4  = sincos(-ω + xli  - G32)
            sin_5,  cos_5  = sincos(+ω + xli  - G52)
            sin_6,  cos_6  = sincos(-ω + xli  - G52)
            sin_7,  cos_7  = sincos(2ω + 2xli - G44)
            sin_8,  cos_8  = sincos(     2xli - G44)
            sin_9,  cos_9  = sincos(+ω + 2xli - G54)
            sin_10, cos_10 = sincos(-ω + 2xli - G54)

            xndot = d2201 * sin_1 + d2211 * sin_2 + d3210 * sin_3 +
                    d3222 * sin_4 + d5220 * sin_5 + d5232 * sin_6 +
                    d4410 * sin_7 + d4422 * sin_8 + d5421 * sin_9 +
                    d5433 * sin_10

            xnddt =   d2201 * cos_1 + d2211 * cos_2 + d3210 * cos_3 +
                      d3222 * cos_4 + d5220 * cos_5 + d5232 * cos_6 +
                    2(d4410 * cos_7 + d4422 * cos_8 + d5421 * cos_9 +
                      d5433 * cos_10)
        end

        xldot  = xni + xfact
        xnddt *= xldot
    end

	# Set up for original mode (LS terms at epoch non-zero).
	pgh0 = ph0 = pe0 = pinc0 = pl0 = 0.0

    @pack! sgp4_ds = atime, xli, xni, xnq, xfact, ssl, ssg, ssh, sse, ssi,
                     xlamo, omegaq, omgdt, gmst, del1, del2, del3, fasx2, fasx4,
                     fasx6, d2201, d2211, d3210, d3222, d4410, d4422, d5220,
                     d5232, d5421, d5433, xnddt, xndot, xldot, zmos, se2, se3,
                     si2, si3, sl2, sl3, sl4, sgh2, sgh3, sgh4, sh2, sh3, zmol,
                     ee2, e3, xi2, xi3, xl2, xl3, xl4, xgh2, xgh3, xgh4, xh2,
                     xh3, pe, pinc, pgh, ph, pl, pgh0, ph0, pe0, pinc0, pl0,
                     isynfl, iresfl, ilsz

    return sgp4_ds
end

"""
    function dssec!(sgp4_ds::SGP4_DeepSpace{T}, nll_0::T, e_0::T, i_0::T, ω_0::T, Ω_k::T, ω_k::T, M_k::T, dotω::T, Δt::Number) where T<:Number

Compute the secular effects.

Notice that the values in the structure `sgp4_ds` **will be modified**.

# Args

* `sgp4_ds`: Deep space structure (see `SGP4_DeepSpace`).
* `nll_0`: Initial mean motion [rad/min].
* `e_0`: Initial eccentricity.
* `i_0`: Initial inclination [rad].
* `ω_0`: Initial argument of perigee [rad].
* `Ω_k`: Current right ascension of the ascending node [rad].
* `ω_k`: Current argument of perigee [rad].
* `M_k`: Current mean anomaly [rad].
* `dotω`: Time-derivative of the argument of perigee [rad/min].
* `Δt`: Time interval since the epoch [min].

# Returns

The following elements perturbed by the secular effects:

* Mean motion [rad/min].
* Eccentricity.
* Inclination [rad].
* Right ascension of the ascending node [rad].
* Argument of perigee [rad].
* Mean anomaly [rad].

"""
function dssec!(sgp4_ds::SGP4_DeepSpace{T},
                nll_0::T,
                e_0::T,
                i_0::T,
                ω_0::T,
                Ω_k::T,
                ω_k::T,
                M_k::T,
                dotω::T,
                Δt::Number) where T<:Number
    # Unpack variables.
    @unpack_SGP4_DeepSpace sgp4_ds

    M_sec = M_k + ssl * Δt
    e_sec = e_0 + sse * Δt
    i_sec = i_0 + ssi * Δt
    Ω_sec = Ω_k + ssh * Δt
    ω_sec = ω_k + ssg * Δt

    # TODO: Verify what this variable means. This is found in `dspace.m` of
    # Vallado's implementation [2].
    θ = mod(gmst + THDT * Δt, 2π)

    # If the orbit is not resonant, then nothing more should be computed.
    (!iresfl) && return (nll_0, e_sec, i_sec, Ω_sec, ω_sec, M_sec)

    #   Update Resonances using Numerical (Euler-Maclaurin) Integration
    # ==========================================================================

    # Epoch restart
    # =============

    # This verification is different between Vallado's [2] and [3]. We will use
    # [2] since it seems more recent.
    if  (atime == 0) || (Δt * atime <= 0) || (abs(Δt) < abs(atime))
        atime = T(0)
        xni   = nll_0
        xli   = xlamo
    end

    # Integration
    # ===========

    ft = Δt - atime

    # In [3], the integration process is performed only if `ft` is larger than
    # `STEP`. However, Vallado's implementation [2] does not verify this and the
    # integration is performed every time. This behavior was chose because it
    # seems that [3] is a more recent version of the algorithm.

    # Check integration direction.
    delt = (Δt >= atime) ? STEP : -STEP

    # Perform the integration with step `delt` until the difference between
    # the time `Δt` and `atime` is less then `STEP`.
    while true
        # Compute the dot terms.
        if (isynfl)

            sin_1, cos_1 = sincos(  (xli - fasx2) )
            sin_2, cos_2 = sincos( 2(xli - fasx4) )
            sin_3, cos_3 = sincos( 3(xli - fasx6) )

            xndot = del1 * sin_1 +  del2 * sin_2 +  del3 * sin_3
            xnddt = del1 * cos_1 + 2del2 * cos_2 + 3del3 * cos_3
        else
            ω     = ω_0 + dotω * atime

            sin_1,  cos_1  = sincos(2ω + xli  - G22)
            sin_2,  cos_2  = sincos(   + xli  - G22)
            sin_3,  cos_3  = sincos(+ω + xli  - G32)
            sin_4,  cos_4  = sincos(-ω + xli  - G32)
            sin_5,  cos_5  = sincos(+ω + xli  - G52)
            sin_6,  cos_6  = sincos(-ω + xli  - G52)
            sin_7,  cos_7  = sincos(2ω + 2xli - G44)
            sin_8,  cos_8  = sincos(     2xli - G44)
            sin_9,  cos_9  = sincos(+ω + 2xli - G54)
            sin_10, cos_10 = sincos(-ω + 2xli - G54)

            xndot = d2201 * sin_1 + d2211 * sin_2 + d3210 * sin_3 +
                    d3222 * sin_4 + d5220 * sin_5 + d5232 * sin_6 +
                    d4410 * sin_7 + d4422 * sin_8 + d5421 * sin_9 +
                    d5433 * sin_10

            xnddt =   d2201 * cos_1 + d2211 * cos_2 + d3210 * cos_3 +
                      d3222 * cos_4 + d5220 * cos_5 + d5232 * cos_6 +
                    2(d4410 * cos_7 + d4422 * cos_8 + d5421 * cos_9 +
                      d5433 * cos_10)
        end

        xldot  = xni + xfact
        xnddt *= xldot

        ft = Δt - atime
        (abs(ft) < STEP) && break

        # In Vallado's implementation [2], this is in the final of the loop
        # instead of at the beginning.
        xli   += delt * (xldot + delt * xndot / 2)
        xni   += delt * (xndot + delt * xnddt / 2)
        atime += delt
    end

    xl    = xli + ft * (xldot + ft * xndot / 2)
    n_sec = xni + ft * (xndot + ft * xnddt / 2)
    M_sec = !isynfl ? xl - 2Ω_sec + 2θ : xl - Ω_sec - ω_sec + θ

    @pack! sgp4_ds = atime, xni, xli, xnddt, xndot, xldot

    (n_sec, e_sec, i_sec, Ω_sec, ω_sec, M_sec)
end

"""
    function dsper!(sgp4_ds::SGP4_DeepSpace{T}, e_k::T, i_k::T, Ω_k::T, ω_k::T, M_k::T, Δt:Number) where T<:Number

Compute the effects caused by Lunar-Solar periodics.

Notice that the values in the structure `sgp4_ds` **will be modified**.

# Args

* `sgp4_ds`: Deep space structure (see `SGP4_DeepSpace`).
* `e_k`: Current eccentricity.
* `i_k`: Current inclination [rad].
* `Ω_k`: Current right ascension of the ascending node [rad].
* `ω_k`: Current argument of perigee [rad].
* `M_k`: Current mean anomaly [rad].
* `Δt`: Time interval since the epoch [min].

# Returns

The following elements perturbed by lunar-solar periodics.

* Eccentricity.
* Inclination [rad].
* Right ascension of the ascending node [rad].
* Argument of perigee [rad].
* Mean anomaly [rad].

"""
function dsper!(sgp4_ds::SGP4_DeepSpace{T},
                e_k::T,
                i_k::T,
                Ω_k::T,
                ω_k::T,
                M_k::T,
                Δt::Number) where T<:Number
    # Unpack variables.
    @unpack_SGP4_DeepSpace sgp4_ds

    #                          Update Solar Terms
    # ==========================================================================

	zm    = zmos +  ZNS * Δt
	zf    = zm   + 2ZES * sin(zm)

    sinzf, coszf = sincos(zf)

	f2    = +sinzf * sinzf / 2 - 0.25
	f3    = -sinzf * coszf / 2
	ses   = se2 * f2 + se3 * f3
	sis   = si2 * f2 + si3 * f3
	sls   = sl2 * f2 + sl3 * f3 + sl4 * sinzf

	sghs  = sgh2 * f2 + sgh3 * f3 + sgh4 * sinzf
	shs   = sh2  * f2 + sh3  * f3

    #                          Update Lunar Terms
    # ==========================================================================

	zm    = zmol +  ZNL * Δt
	zf    = zm   + 2ZEL * sin(zm)

    sinzf, coszf = sincos(zf)

	f2    = +sinzf * sinzf / 2 - 0.25
	f3    = -sinzf * coszf / 2
	sel   = ee2 * f2 + e3 * f3
	sil   = xi2 * f2 + xi3 * f3
	sll   = xl2 * f2 + xl3 * f3 + xl4 * sinzf

	sghl  = xgh2 * f2 + xgh3 * f3 + xgh4 * sinzf
	shl   = xh2  * f2 + xh3  * f3

    #                         Save computed values
    # ==========================================================================

	pgh  = sghs + sghl
	ph   = shs  + shl
	pe   = ses  + sel
	pinc = sis  + sil
	pl   = sls  + sll

    # Update inclination and eccentricity.
    e_per = e_k + pe
    i_per = i_k + pinc

    sinis, cosis = sincos(i_per)

    # The original algorithm considered the original inclination to select the
    # Lyddane Lunar-Solar perturbations algorithm. However, Vallado's
    # implementation [2] test the perturbed inclination to select this. It is
    # mentioned that this is the behavior selected in GSFC source code.
    if i_per >= 0.2
        tmp_ph = ph / sinis;
		ω_per  = ω_k + pgh - cosis * tmp_ph;
		Ω_per  = Ω_k + tmp_ph;
		M_per  = M_k + pl;
    else
        sinok   = sin(Ω_k)
        cosok   = cos(Ω_k)

        #                       |----------    dalf     ----------|
		alfdp   = sinis * sinok + ph * cosok + pinc * cosis * sinok
        #                       |----------    dbet     ----------|
		betdp   = sinis * cosok - ph * sinok + pinc * cosis * cosok

        # For the following computation, in which `Ω_per` is used without a
        # trigonometric function, it is advisable to make sure that it stays in
        # the interval [0, 2π].
        Ω_per   = mod(Ω_k, 2π)

        #                                   |----------    dls    ----------|
		xls     = M_k + ω_k + cosis * Ω_per + pl + pgh - pinc * Ω_per * sinis
        Ω_aux   = Ω_per
        Ω_per   = mod(atan(alfdp, betdp), 2π)

        if abs(Ω_aux - Ω_per) > π
            Ω_per = (Ω_per < Ω_aux) ? Ω_per + 2π : Ω_per - 2π
        end

		M_per   = M_k + pl;
		ω_per   = xls - M_per - cosis * Ω_per
    end

    @pack! sgp4_ds = pgh, ph, pe, pinc, pl

    (e_per, i_per, Ω_per, ω_per, M_per)
end
