# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Functions to compute general values related to the orbit.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export orbital_angular_velocity
export orbital_angular_velocity_to_semimajor_axis
export orbital_period
export raan_time_derivative

############################################################################################
#                                        Functions
############################################################################################

"""
    orbital_angular_velocity(a::Number, e::Number, i::Number; kwargs...) -> T
    orbital_angular_velocity(orb::Orbit{Tepoch, T}; kwargs...) where {Tepoch<:Number, T<:Number} -> T

Compute the angular velocity [rad/s] of an object in an orbit with semi-major axis `a` [m],
eccentricity `e`, and inclination `i` [rad]. The orbit can also be specified by `orb` (see
`Orbit`).

!!! note
    The output type `T` in the first signature is obtained by promoting the inputs to a
    float type.

# Keyword

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used.
    (**Default**: `:J2`)

# Perturbations

The keyword argument `perturbation` can be used to select the perturbation terms that will
be considered in the computation. The possible values are:

- `:J0`: Consider a Keplerian orbit.
- `:J2`: Consider the perturbation terms up to J2.
- `:J4`: Consider the perturbation terms J2, J4, and J2².

If `perturbation` is omitted, it defaults to `:J2`.
"""
function orbital_angular_velocity(
    a::T1,
    e::T2,
    i::T3;
    perturbation::Symbol = :J2
) where {T1<:Number, T2<:Number, T3<:Number}
    T = float(promote_type(T1, T2, T3))

    # Unperturbed orbit period.
    n₀ = √(T(GM_EARTH) / T(a)^3)

    # Perturbation computed using a Keplerian orbit.
    if perturbation === :J0
        return n₀

    # Perturbation computed using perturbations terms up to J2.
    elseif perturbation === :J2
        # Auxiliary variables.
        cos_i² = cos(T(i))^2
        aux    = 1 - T(e)^2
        R₀     = T(WGS84_ELLIPSOID.a)

        # Semi-lactum rectum.
        p = T(a) * aux

        # Orbit period considering the perturbations (up to J2).
        return n₀ + 3R₀^2 * T(EGM08_J2) / (4p^2) * n₀ * (√aux * (3cos_i² - 1) + (5cos_i² - 1))

    # Perturbation computed using perturbations terms J2, J4, and J2².
    elseif perturbation === :J4

        # Auxiliary variables
        e²     = T(e)^2
        e⁴     = T(e)^4
        sin_i  = sin(T(i))
        sin_i² = sin_i^2
        sin_i⁴ = sin_i^4
        aux    = (1 - e²)
        saux   = √aux
        R₀     = T(WGS84_ELLIPSOID.a)
        p      = (T(a) / R₀) * aux
        p²     = p^2
        p⁴     = p^4
        T_J2   = T(EGM08_J2)
        T_J2²  = T_J2^2
        T_J4   = T(EGM08_J4)

        # Notice that:
        #            .   .
        #   n = n₀ + ω + M₀
        #
        # in which the time-derivatives are computed as in [1, p. 692].

        δω  = ( 3 //   4) * n₀ * T_J2  / p² * (4 - 5sin_i²) +
              ( 9 // 384) * n₀ * T_J2² / p⁴ * (     56e² + (760 -  36e²) * sin_i² - (890 +  45e²) * sin_i⁴) -
              (15 // 128) * n₀ * T_J4  / p⁴ * (64 + 72e² - (248 + 252e²) * sin_i² + (196 + 189e²) * sin_i⁴)

        δM₀ = ( 3 //   4) * n₀ * T_J2  / p² * saux * (2 - 3sin_i²) +
              ( 3 // 512) * n₀ * T_J2² / p⁴ / saux * (320e² - 280e⁴ + (1600 - 1568e² + 328e⁴) * sin_i² + (-2096 + 1072e² +  79e⁴) * sin_i⁴) -
              (45 // 128) * n₀ * T_J4  / p⁴ * saux * e² * (-8 + 40sin_i - 35sin_i²)

        return n₀ + δω + δM₀
    else
        throw(ArgumentError("The perturbation parameter :$perturbation is invalid."))
    end
end

function orbital_angular_velocity(orb::Orbit; perturbation::Symbol = :J2)
    # Convert first to Keplerian elements.
    k = convert(KeplerianElements, orb)
    return orbital_angular_velocity(k.a, k.e, k.i; perturbation = perturbation)
end

"""
    orbital_angular_velocity_to_semimajor_axis(n::Number, e::Number, i::Number; kwargs...) -> T

Compute the semi-major axis [m] that will provide an angular velocity `n` [rad / s] in an
orbit with eccentricity `e` and inclination `i` [rad].

Notice that the angular velocity `n` is related to the nodal period, *i.e.* the time between
two consecutive passages by the ascending node.

!!! note
    The output type `T` in the first signature is obtained by promoting the inputs to a
    float type.

# Keyword

- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default** = `GM_EARTH`)
- `max_iterations::Int`: Maximum number of iterations allowed in the Newton-Raphson
    algorithm. (**Default** = 20)
- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used.
    (**Default**: `:J2`)
- `tolerance::Number`: Tolerance to stop the Newton-Raphson algorithm. (**Default** = 1e-10)

# Perturbations

The keyword argument `perturbation` can be used to select the perturbation terms that will
be considered in the computation. The possible values are:

- `:J0`: Consider a Keplerian orbit.
- `:J2`: Consider the perturbation terms up to J2.
- `:J4`: Consider the perturbation terms J2, J4, and J2².

If `perturbation` is omitted, it defaults to `:J2`.
"""
function orbital_angular_velocity_to_semimajor_axis(
    n::T1,
    e::T2,
    i::T3;
    m0::Number = GM_EARTH,
    max_iterations::Int = 20,
    perturbation::Symbol = :J2,
    tolerance::Number = 1e-10
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    T = float(promote_type(T1, T2, T3))

    if perturbation === :J0

        a = (T(m0) / T(n)^2)^(1 // 3)
        return a

    elseif perturbation === :J2
        # Get the semi-major axis [m] that will provide the mean motion `n` using
        # perturbation terms up to J2.
        #
        # This can only be done using a numerical algorithm to solve the following equation
        # for `a`:
        #
        #          ┌                                                              ┐
        #          │      3         √(1 - e²) . (3cos²(i) - 1) + (5cos²(i) - 1))  │
        #   n = n₀ │ 1 + ─── . J₂ . ────────────────────────────────────────────  │
        #          │      4                         a^2.(1-e^2)^2                 │
        #          └                                                              ┘
        #
        #           √μ
        #   n₀ =  ──────
        #         √(a^3)
        #
        # To improve algorithm stability, we will compute the normalized semi-major axis,
        # i.e. `a / R₀`.

        # Auxiliary variables to solve for the semi-major axis.
        R₀     = T(WGS84_ELLIPSOID.a)
        sqrt_μ = √(T(m0) / T(R₀)^3)
        cos_i  = cos(T(i))
        K      = (3 // 4) * T(EGM08_J2) * (√(1 - T(e)^2) * (3cos_i^2 - 1) + (5cos_i^2 - 1)) / (1 - T(e)^2)^2

        # Initial guess using a non-perturbed orbit.
        a = (T(m0) / T(n)^2)^(1 // 3) / R₀

        # Newton-Raphson algorithm
        # ==================================================================================
        #
        # Notice that we will allow, at most, `max_iterations`.
        for k in 1:max_iterations
            # Auxiliary variables.
            ap3o2  = √(a^3)    # -> a^( 2 / 3)
            ap5o2  = ap3o2 * a # -> a^( 5 / 2)
            ap7o2  = ap5o2 * a # -> a^( 7 / 2)
            ap11o2 = ap7o2 * a # -> a^(11 / 2)

            # Compute the residue.
            res = n - sqrt_μ / ap3o2 - K * sqrt_μ / ap7o2

            # Compute the Jacobian of the function.
            df = (3 // 2) * sqrt_μ / ap5o2 + (7 // 2) * K * sqrt_μ / ap11o2

            # Compute the new estimate.
            a = a - res / df

            (abs(res) < tolerance) && break
        end

        return a * R₀

    elseif perturbation === :J4

        # Auxiliary variables
        R₀     = T(WGS84_ELLIPSOID.a)
        sqrt_μ = √(T(m0) / T(R₀)^3)
        sin_i  = sin(T(i))
        sin_i² = sin_i^2
        sin_i⁴ = sin_i^3
        e²     = T(e)^2
        e⁴     = T(e)^4
        aux    = 1 - e²
        saux   = √aux

        # Get the semi-major axis using J4 perturbation theory [er].
        #
        # This can only be done using a numerical algorithm to solve the following equation
        # for `a`:
        #
        #            .   .
        #   n = n₀ + ω + M₀ ,
        #
        #           √μ
        #   n₀ =  ──────
        #         √(a^3)
        #
        # and the time-derivatives of the argument of perigee and the initial mean anomaly
        # is compute considering the terms J2, J4, and J2².
        #
        # To improve algorithm stability, we will compute the normalized semi-major axis,
        # i.e. `a / R₀`.

        K₁ = (3 // 4) * T(EGM08_J2) / aux^2 * ((4 - 5sin_i²) + √(1 - e²) * (2 - 3sin_i²))

        K₂ = (3 // 512) * T(EGM08_J2)^2 / aux^4 * (
            224e² + (3040 - 144e²) * sin_i² - (3560 + 180e²) * sin_i⁴ + 1 / saux * (
                (         320e² - 280e⁴) +
                ( 1600 - 1568e² + 328e⁴) * sin_i² +
                (-2096 + 1072e² +  79e⁴) * sin_i⁴
            )
        )

        K₃ = (15 // 128) * T(EGM08_J4) / aux^4 * (
            64 + 72e² -
            (248 + 252e²) * sin_i² +
            (196 + 189e²) * sin_i⁴ +
            e² / saux * (-8 + 40sin_i - 35sin_i²)
        )

        # Initial guess using a non-perturbed orbit.
        a = (T(m0) / T(n)^2)^(1 // 3) / R₀

        # Newton-Raphson algorithm
        # ==================================================================================
        #
        # Notice that we will allow, at most, `max_iterations`.
        for k = 1:max_iterations
            # Auxiliary variables.
            ap3o2  = √(a^3)      # -> a^( 3 / 2)
            ap5o2  =  ap3o2 * a  # -> a^( 5 / 2)
            ap7o2  =  ap5o2 * a  # -> a^( 7 / 2)
            ap9o2  =  ap7o2 * a  # -> a^( 9 / 2)
            ap11o2 =  ap9o2 * a  # -> a^(11 / 2)
            ap13o2 = ap11o2 * a  # -> a^(13 / 2)

            # Compute the residue.
            res = n - sqrt_μ * (1 / ap3o2 + K₁ / ap7o2 + (K₂ + K₃) / ap11o2)

            # Compute the Jacobian of the function.
            df = sqrt_μ / 2 * (3ap5o2 + 7K₁ / ap9o2 + 11 * (K₂ + K₃) / ap13o2)

            # Compute the new estimate.
            a = a - res / df

            (abs(res) < tolerance) && break
        end

        return a * R₀
    else
        throw(ArgumentError("The perturbation parameter :$perturbation is invalid."))
    end
end

"""
    orbital_period(a::Number, e::Number, i::Number; kwargs...) -> T
    orbital_period(orb::Orbit{Tepoch, T}; kwargs...) where {Tepoch<:Number, T<:Number} -> T

Compute the orbital period [s] of an object in an orbit with semi-major axis `a` [m],
eccentricity `e`, and inclination `i` [rad]. The orbit can also be specified by `orb` (see
`Orbit`).

!!! note
    The output type `T` in the first signature is obtained by promoting the inputs to a
    float type.

# Keyword

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used.
    (**Default**: `:J2`)

# Perturbations

The keyword argument `perturbation` can be used to select the perturbation terms that will
be considered in the computation. The possible values are:

- `:J0`: Consider a Keplerian orbit.
- `:J2`: Consider the perturbation terms up to J2.
- `:J4`: Consider the perturbation terms J2, J4, and J2².

If `perturbation` is omitted, it defaults to `:J2`.
"""
function orbital_period(a::Number, e::Number, i::Number; perturbation::Symbol = :J2)
    n = orbital_angular_velocity(a, e, i; perturbation = perturbation)
    T = eltype(n)
    return T(2π) / n
end

function orbital_period(orb::Orbit; perturbation::Symbol = :J2)
    # Convert first to Keplerian elements.
    k = convert(KeplerianElements, orb)
    return orbital_period(k.a, k.e, k.i; perturbation = perturbation)
end

"""
    raan_time_derivative(a::Number, e::Number, i::Number; kwargs...) -> T

Compute the time derivative of the right ascension of the ascending node (RAAN) [rad / s] in
an orbit with semi-major axis `a` [m], eccentricity `e`, and inclination `i` [rad]. The
orbit can also be specified by `orb` (see `Orbit`).

!!! note
    The output type `T` in the first signature is obtained by promoting the inputs to a
    float type.

# Keyword

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used.
    (**Default**: `:J2`)

# Perturbations

The keyword argument `perturbation` can be used to select the perturbation terms that will
be considered in the computation. The possible values are:

- `:J0`: Consider a Keplerian orbit.
- `:J2`: Consider the perturbation terms up to J2.
- `:J4`: Consider the perturbation terms J2, J4, and J2².

If `perturbation` is omitted, it defaults to `:J2`.
"""
function raan_time_derivative(
    a::T1,
    e::T2,
    i::T3;
    perturbation::Symbol = :J2
) where {T1<:Number, T2<:Number, T3<:Number}
    T = float(promote_type(T1, T2, T3))

    # Perturbation computed using a Keplerian orbit.
    if perturbation === :J0
        return zero(T)

    # Perturbation computed using perturbations terms up to J2.
    elseif perturbation === :J2
        # Auxiliary variables.
        J₂ = T(EGM08_J2)
        R₀ = T(WGS84_ELLIPSOID.a)
        e² = e^2
        p  = (a / R₀) * (1 - e²)
        p² = p^2

        sin_i, cos_i = sincos(i)

        # Unperturbed mean motion.
        n₀ = orbital_angular_velocity(a, e, i; perturbation = :J0)

        # Perturbed orbit mean motion.
        n̄ = n₀ * (1 + (3 // 4) * J₂ / p² * √(1 - e²) * (2 - 3sin_i^2))

        # First-order time-derivative of the RAAN [rad / s].
        ∂Ω = -(3 // 2) * n̄ * J₂ / p² * cos_i

        return ∂Ω

    # Perturbation computed using perturbation terms J₂, J₂², and J₄.
    elseif perturbation === :J4
        # Auxiliary variables.
        J₂   = T(EGM08_J2)
        J₂²  = J₂^2
        J₄   = T(EGM08_J4)
        R₀   = T(WGS84_ELLIPSOID.a)
        e²   = e^2
        p    = (a / R₀) * (1 - e²)
        p²   = p^2
        p⁴   = p^4
        aux  = 1 - e²
        saux = √aux

        sin_i, cos_i = sincos(i)
        sin_i² = sin_i^2

        # Unperturbed orbit period.
        n₀ = orbital_angular_velocity(a, e, i; perturbation = :J0)

        # Perturbed mean motion.
        ā = a * (1 - (3 // 4) * J₂ / p² * saux * (2 - 3sin_i²))
        p̄ = (ā / R₀) * aux
        n̄ = n₀ * (1 + (3 // 4) * J₂ / p̄^2 * saux * (2 - 3sin_i²))

        # First-order time-derivative of the RAAN [rad / s].
        δΩ = -( 3 // 2 ) * n̄ * J₂  / p² * cos_i +
              ( 3 // 32) * n̄ * J₂² / p⁴ * cos_i * (-36 -  4e² + 48saux + (40 - 5e² - 72saux) * sin_i²) +
              (15 // 32) * n̄ * J₄  / p⁴ * cos_i * (  8 + 12e² - (14 + 21e²) * sin_i²)

        return δΩ
    else
        throw(ArgumentError("The perturbation parameter $perturbation is not defined."))
    end
end

function raan_time_derivative(orb::Orbit; perturbation::Symbol = :J2)
    # Convert first to Keplerian elements.
    k = convert(KeplerianElements, orb)
    return raan_time_derivative(k.a, k.e, k.i; perturbation = perturbation)
end
