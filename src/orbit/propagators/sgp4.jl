#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#
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
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2017-08-08: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export SGP4_GravCte, sgp4_gc_wgs84
export sgp4_init, sgp4!

################################################################################
#                             Types and Structures
################################################################################

# Gravitational constants for SGP4.
immutable SGP4_GravCte
    R0   # Earth equatorial radius [km].
    XKE  # sqrt(GM) [er/min]^(3/2).
    J2   # The second gravitational zonal harmonic of the Earth.
    J3   # The third  gravitational zonal harmonic of the Earth.
    J4   # The fourth gravitational zonal harmonic of the Earth.
end

# Serialization of the arguments in SGP4_GravCte.
function getindex(sgp4_gc::SGP4_GravCte, ::Colon)
    sgp4_gc.R0, sgp4_gc.XKE, sgp4_gc.J2, sgp4_gc.J3, sgp4_gc.J4
end

# SPG4 structure.
type SGP4_Structure
    # TLE parameters.
    t_0::Float64
    n_0::Float64
    e_0::Float64
    i_0::Float64
    Ω_0::Float64
    ω_0::Float64
    M_0::Float64
    bstar::Float64
    # Current parameters.
    a_k::Float64
    e_k::Float64
    i_k::Float64
    Ω_k::Float64
    ω_k::Float64
    M_k::Float64
    n_k::Float64
    # Parameters related with the orbit.
    all_0::Float64
    nll_0::Float64
    # Useful constants to decrease the computational burden.
    AE::Float64
    QOMS2T::Float64
    β_0::Float64
    ξ::Float64
    η::Float64
    sin_i_0::Float64
    θ::Float64
    θ2::Float64
    θ3::Float64
    θ4::Float64
    A_30::Float64
    k_2::Float64
    k_4::Float64
    C1::Float64
    C3::Float64
    C4::Float64
    C5::Float64
    D2::Float64
    D3::Float64
    D4::Float64
    # Others.
    isimp::Bool
    # SGP4 gravitational constants.
    sgp4_gc::SGP4_GravCte
end

# Serialization of the arguments in SGP4_Structure.
function getindex(sgp4d::SGP4_Structure, ::Colon)

    sgp4d.t_0, sgp4d.n_0, sgp4d.e_0, sgp4d.i_0, sgp4d.Ω_0, sgp4d.ω_0, sgp4d.M_0,
    sgp4d.bstar,sgp4d.a_k, sgp4d.e_k, sgp4d.i_k, sgp4d.Ω_k, sgp4d.ω_k,
    sgp4d.M_k, sgp4d.n_k, sgp4d.all_0, sgp4d.nll_0, sgp4d.AE, sgp4d.QOMS2T,
    sgp4d.β_0, sgp4d.ξ, sgp4d.η, sgp4d.sin_i_0, sgp4d.θ, sgp4d.θ2, sgp4d.θ3,
    sgp4d.θ4, sgp4d.A_30, sgp4d.k_2, sgp4d.k_4, sgp4d.C1, sgp4d.C3, sgp4d.C4,
    sgp4d.C5, sgp4d.D2, sgp4d.D3, sgp4d.D4, sgp4d.isimp, sgp4d.sgp4_gc

end

################################################################################
#                                  Constants
################################################################################

# WGS-84 / EGM-08 Gravitational constants.
sgp4_gc_wgs84 = SGP4_GravCte(
        R0/1000,
        60.0/sqrt(6378.137^3/398600.5),
         0.00108262998905,
        -0.00000253215306,
        -0.00000161098761
       )

################################################################################
#                                  Functions
################################################################################

"""
### function sgp4_init(spg4_gc::SGP4_GravCte, t_0::Number, n_0::Number, e_0::Number, i_0::Number, Ω_0::Number, ω_0::Number, M_0::Number, bstar::Number)

Initialize the data structure of SGP4 algorithm.

##### Args

* spg4_gc: SPG4 gravitational constants (see `SGP4_GravCte`).
* t_0: Epoch of the orbital elements [min].
* n_0: SGP type "mean" mean motion at epoch [rad/min].
* e_0: "Mean" eccentricity at epoch.
* i_0: "Mean" inclination at epoch [rad].
* Ω_0: "Mean" longigute of the ascending node at epoch [rad].
* ω_0: "Mean" argument of perigee at epoch [rad].
* M_0: "Mean" mean anomally at epoch [rad].
* bstar: Drag parameter (B*).

##### Returns

The structure `SGP4_Structure` with the intialized parameters.

"""

function sgp4_init(sgp4_gc::SGP4_GravCte,
                   t_0::Number,
                   n_0::Number,
                   e_0::Number,
                   i_0::Number,
                   Ω_0::Number,
                   ω_0::Number,
                   M_0::Number,
                   bstar::Number)

    # Unpack the gravitational constants to improve code readability.
    R0, XKE, J2, J3, J4 = sgp4_gc[:]

    # Constants
    # =========
    #
    # Note: [er] = Earth radii.

    # Distance units / Earth radii.
    AE = 1.0

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

    e2_0 = e_0^2

    sin_i_0 = sin(i_0)

    θ   = cos(i_0)
    θ2  = θ^2
    θ3  = θ^3
    θ4  = θ^4

    # ==========================================================================

    # Recover the original mean motion (nll_0) and semi-major axis (all_0) from
    # the input elements.

    aux = (3*θ2-1)/(1-e2_0)^(3/2)

    a_1 = (XKE/n_0)^(2/3)
    δ_1 = 3/2*k_2/a_1^2*aux
    a_0 = a_1*(1 - 1/3*δ_1 - δ_1^2 - 134/81*δ_1^3)
    δ_0 = 3/2*k_2/a_0^2*aux

    nll_0 = n_0/(1 + δ_0)
    all_0 = a_0/(1 - δ_0)

    # Initialization
    # ==============

    # For perigee less than 220 km, the `isimp` flag is set and the equations
    # are truncated to a linear variation in sqrt(a) and quadratic variation in
    # mean anomaly. Also, the C5 term, the δω term, and the δM term are dropped.
    isimp = ( (all_0*(1-e_0)/AE) < (220/XKMPER+AE) )

    # For perigee below 156 km, the values of S and QOMS2T are altered.
    perigee = (all_0*(1-e_0)-AE)*XKMPER

    if perigee < 156
        if perigee < 98
            s = 20.0/XKMPER + AE
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

    aux1 = (1-η^2)^(-7/2)
    aux2 = ξ^4*all_0*β_0^2*aux1

    C2 = QOMS2T*ξ^4*nll_0*aux1*
         ( all_0*( 1 + (3/2)η^2 + 4e_0*η + e_0*η^3) +
           3/2*(k_2*ξ)/(1-η^2)*(-1/2 + (3/2)θ2)*(8 + 24η^2 + 3η^4) )

    C1 = bstar*C2

    C3 = QOMS2T*ξ^5*A_30*nll_0*AE*sin_i_0/(k_2*e_0)

    C4 = 2nll_0*QOMS2T*aux2*
         ( 2η*(1+e_0*η) + (1/2)e_0 + (1/2)η^3 -
           2*k_2*ξ/(all_0*(1-η^2))*
                (3*(1-3θ2)*(1 + (3/2)η^2 - 2*e_0*η - (1/2)e_0*η^3) +
                3/4*(1-θ2)*(2η^2 - e_0*η - e_0*η^3)*cos(2*ω_0)))

    C5 = 2*QOMS2T*aux2*(1 + (11/4)η*(η+e_0) + e_0*η^3)

    D2 = 4*all_0*ξ*C1^2

    D3 = 4/3*all_0*ξ^2*( 17all_0 +   s)*C1^3

    D4 = 2/3*all_0*ξ^3*(221all_0 + 31s)*C1^4

    # The current orbital parameters are obtained from the TLE.
    a_k = all_0
    e_k = e_0
    i_k = i_0
    Ω_k = Ω_0
    ω_k = ω_0
    M_k = M_0
    n_k = nll_0

    # Create the output structure with the data.
    SGP4_Structure(

        t_0, n_0, e_0, i_0, Ω_0, ω_0, M_0, bstar, a_k, e_k, i_k, Ω_k, ω_k, M_k,
        n_k, all_0, nll_0, AE, QOMS2T, β_0, ξ, η, sin_i_0, θ, θ2, θ3, θ4, A_30,
        k_2, k_4, C1, C3, C4, C5, D2, D3, D4, isimp, sgp4_gc

       )
end

"""
### function sgp4!(sgp4d::SGP4_Structure, t::Number)

Propagate the orbit defined in `sgp4d` until the time `t`. Notice that the
values in `sgp4d` will be modified.

##### Args

* sgp4d: SPG4 structure (see `SGP4_Structure`).
* t: Time that the elements will be propagated [min].

##### Returns

* The position vector at time `t` [km].
* The velocity vector at time `t` [km/s].

"""

function sgp4!(sgp4d::SGP4_Structure, t::Number)
    # Unpack variables.
    t_0, n_0, e_0, i_0, Ω_0, ω_0, M_0, bstar, a_k, e_k, i_k, Ω_k, ω_k, M_k, n_k,
    all_0, nll_0, AE, QOMS2T, β_0, ξ, η, sin_i_0, θ, θ2, θ3, θ4, A_30, k_2, k_4,
    C1, C3, C4, C5, D2, D3, D4, isimp, sgp4_gc = sgp4d[:]

    R0, XKE, J2, J3, J4 = sgp4_gc[:]

    # Time elapsed since epoch.
    Δt = t-t_0

    # Secular effects of atmospheric drag and gravitation.
    # ====================================================

    M_DF = M_0 + ( 1 + 3*k_2*(-1 + 3θ2)/(2all_0^2*β_0^3) +
                   3k_2^2*(13 - 78θ2 + 137θ4)/(16all_0^4*β_0^7) )*nll_0*Δt

    ω_DF = ω_0 + ( -3*k_2*(1-5θ2)/(2all_0^2*β_0^4) +
                    3*k_2^2*(7 - 114θ2 + 395θ4)/(16all_0^4*β_0^8) +
                    5*k_4*(3 - 36θ2 + 49θ4)/(4all_0^4*β_0^8) )*nll_0*Δt

    Ω_DF = Ω_0 + ( -3*k_2*θ/(all_0^2*β_0^4) +
                    3*k_2^2*(4θ - 19θ3)/(2all_0^4*β_0^8) +
                    5*k_4*θ*(3 - 7θ2)/(2all_0^4*β_0^8) )*nll_0*Δt

    Ω    = Ω_DF - 21/2*(nll_0*k_2*θ)/(all_0^2*β_0^2)*C1*Δt^2

    # Check  if perigee is below 220 km.
    if !isimp
        δω  = bstar*C3*cos(ω_0)*Δt

        δM  = -2/3*QOMS2T*bstar*ξ^4*AE/(e_0*η)*( (1 + η*cos(M_DF))^3 -
                                                 (1 + η*cos(M_0) )^3 )

        M_p = M_DF + δω + δM

        ω   = ω_DF - δω - δM

        e   = e_0 - bstar*C4*Δt - bstar*C5*(sin(M_p) - sin(M_0))

        a   = all_0*(1 - C1*Δt - D2*Δt^2 - D3*Δt^3 - D4*Δt^4)^2

        IL  = M_p + ω + Ω + nll_0*
              ( (3/2)C1*Δt^2 + (D2 + 2C1^2)*Δt^3 +
                (1/4)*(3D3 + 12C1*D2 + 10C1^3)*Δt^4 +
                (1/5)*(3D4 + 12C1*D3 + 6*D2^2 + 30*C1^2*D2 + 15*C1^4)*Δt^5 )
    else
        # If so, then
        #     1. Drop all terms after C1 in `a` and `IL`.
        #     2. Drop all terms involving C5.
        #     3. Drop δω.
        #     4. Drop δM.
        M_p = M_DF

        ω   = ω_DF

        e   = e_0 - bstar*C4*Δt

        a   = all_0*(1 - C1*Δt)^2

        IL  = M_p + ω + Ω + nll_0*(3/2)C1*Δt^2
    end

    β = sqrt(1-e_0^2)

    # Compute the angular velocity [rad/min].
    n = XKE/a^(3/2)

    # Keep the angles between [0, 2π].
    Ω  = mod(Ω,  2*pi)
    ω  = mod(ω,  2*pi)
    IL = mod(IL, 2*pi)

    # Long-period periodic terms.
    # ===========================

    a_xN = e*cos(ω)

    a_yNL = A_30*sin_i_0/(4*k_2*a*β^2)

    a_yN = e*sin(ω) + a_yNL

    IL_L =  (1/2)a_yNL*a_xN*(3 + 5θ)/(1 + θ)

    IL_T = IL + IL_L

    # Solve Kepler's equation for (E + ω).
    # ====================================

    U = mod(IL_T - Ω, 2*pi)

    E_ω = U

    # Define the following variables that will be modified inside the loop so
    # that we can use them after the loop.
    sin_E_ω = 0.0
    cos_E_ω = 0.0

    for k = 1:10
        sin_E_ω = sin(E_ω)
        cos_E_ω = cos(E_ω)

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

    p_L     = a*(1-e_L^2)

    r       = a*(1-e_cos_E)

    dot_r   = XKE*sqrt(a)*e_sin_E/r

    r_dot_f = XKE*sqrt(p_L)/r

    aux     = e_sin_E/(1 + sqrt(1 - e_L^2))
    cos_u   = a/r*( cos_E_ω - a_xN + a_yN*aux )
    sin_u   = a/r*( sin_E_ω - a_yN - a_xN*aux )
    cos_2u  = 1 - 2*sin_u^2
    sin_2u  = 2cos_u*sin_u
    u       = atan2(sin_u, cos_u)

    # Short-term periodic terms.

    Δr       = +k_2/(2p_L)*(1 - θ2)*cos_2u

    Δu       = -k_2/(4p_L^2)*(7θ2 - 1)*sin_2u

    ΔΩ       = +3k_2*θ/(2p_L^2)*sin_2u

    Δi       = +3k_2*θ/(2p_L^2)*sin_i_0*cos_2u

    Δdot_r   = -k_2*n/p_L*(1 - θ2)*sin_2u

    Δr_dot_f = +k_2*n/p_L*( (1 - θ2)*cos_2u - 3/2*(1 - 3θ2) )

    # The short-term periodics are added to give the osculating quantities.

    r_k       = r*( 1 - (3/2)k_2*sqrt(1 - e_L^2)/p_L^2*(3θ^2 - 1) ) + Δr

    u_k       = u + Δu

    Ω_k       = Ω + ΔΩ

    i_k       = i_0 + Δi

    dot_r_k   = dot_r + Δdot_r

    r_dot_f_k = r_dot_f + Δr_dot_f

    # Orientation vectors.
    sin_Ω_k = sin(Ω_k)
    cos_Ω_k = cos(Ω_k)
    sin_i_k = sin(i_k)
    cos_i_k = cos(i_k)
    sin_u_k = sin(u_k)
    cos_u_k = cos(u_k)

    M = [ -sin_Ω_k*cos_i_k;
           cos_Ω_k*cos_i_k;
                   sin_i_k; ]

    N = [ +cos_Ω_k;
          +sin_Ω_k;
          +0.0 ]

    U = M*sin_u_k + N*cos_u_k
    V = M*cos_u_k - N*sin_u_k

    r = r_k*U*R0
    v = (dot_r_k*U + r_dot_f_k*V)*R0/60.0

    # Update the variables.
    sgp4d.a_k = a
    sgp4d.e_k = e
    sgp4d.i_k = i_k
    sgp4d.Ω_k = Ω_k
    sgp4d.ω_k = ω
    sgp4d.M_k = M_p
    sgp4d.n_k = n

    (r, v)
end

