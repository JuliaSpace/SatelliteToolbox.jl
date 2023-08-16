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

# Keywords

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used.
    (**Default**: `:J2`)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default** = `GM_EARTH`)
- `J2::Number`: J₂ perturbation term. (**Default** = EGM08_J2)
- `J4::Number`: J₄ perturbation term. (**Default** = EGM08_J4)
- `R0::Number`: Earth's equatorial radius [m]. (**Default** = EARTH_EQUATORIAL_RADIUS)

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
    perturbation::Symbol = :J2,
    # Constants.
    J2::Number = EGM08_J2,
    J4::Number = EGM08_J4,
    m0::Number = GM_EARTH,
    R0::Number = EARTH_EQUATORIAL_RADIUS
) where {T1<:Number, T2<:Number, T3<:Number}
    T = float(promote_type(T1, T2, T3))

    R₀ = T(R0)
    μ  = T(m0)
    J₂ = T(J2)
    J₄ = T(J4)

    # Unperturbed orbit period.
    n₀ = √(μ / T(a)^3)

    # Perturbation computed using a Keplerian orbit.
    if perturbation == :J0
        return n₀

    # Perturbation computed using perturbations terms up to J2.
    elseif perturbation == :J2
        # Convert the inputs to the correct type.
        a₀ = T(a)
        e₀ = T(e)
        i₀ = T(i)

        # Initial values and auxiliary variables.
        al₀ = a₀ / R₀          # ........................... Normalized semi-major axis [er]
        e₀² = e₀^2             # .................................. Eccentricity squared [ ]
        p₀  = al₀ * (1 - e₀²)  # .................................... Semi-latus rectum [er]
        p₀² = p₀^2             # ........................... Semi-latus rectum squared [er²]

        sin_i₀, cos_i₀ = sincos(T(i₀))
        sin_i₀² = sin_i₀^2

        # We defined the orbit angular velocity here based on the nodal period, i.e., the
        # time it takes for the satellite to cross the ascending node two consecutive times.
        # Hence, we can compute it by:
        #
        #             ∂M     ∂ω
        #   angvel = ──── + ────.
        #             ∂t     ∂t
        #
        # The expressions for those time-derivatives were obtained from the J2 orbit
        # propagator of SatelliteToolboxPropagators.jl package.

        n̄     = n₀ * (1 + (3 // 4) * J₂ / p₀² * √(1 - e₀²) * (2 - 3sin_i₀²))
        ∂M_∂t = n̄
        ∂ω_∂t = +(3 // 4) * n̄ * J₂ / p₀² * (4 - 5sin_i₀²)

        # Angular velocity.
        ang_vel = ∂M_∂t + ∂ω_∂t

        return ang_vel

    # Perturbation computed using perturbations terms J2, J4, and J2².
    elseif perturbation == :J4
        # Convert the inputs to the correct type.
        a₀ = T(a)
        e₀ = T(e)
        i₀ = T(i)

        # Initial values and auxiliary variables.
        e₀² = e₀^2
        β²  = (1 - e₀²)
        β   = √β²

        al₀ = a₀ / R₀    # ................................. Normalized semi-major axis [er]
        J₂² = J₂^2       # ............................................. J2 constant squared
        p₀  = al₀ * β²   # .......................................... Semi-latus rectum [er]
        p₀² = p₀^2       # ................................. Semi-latus rectum squared [er²]
        p₀⁴ = p₀^4       # ........................ Semi-latus rectum to the 4th power [er⁴]

        sin_i₀, cos_i₀ = sincos(T(i₀))

        sin_i₀² = sin_i₀^2
        sin_i₀⁴ = sin_i₀^4
        cos_i₀⁴ = cos_i₀^4

        # We defined the orbit angular velocity here based on the nodal period, i.e., the
        # time it takes for the satellite to cross the ascending node two consecutive times.
        # Hence, we can compute it by:
        #
        #             ∂M     ∂ω
        #   angvel = ──── + ────.
        #             ∂t     ∂t
        #
        # The expressions for those time-derivatives were obtained from the J2 orbit
        # propagator of SatelliteToolboxPropagators.jl package.

        n̄ = n₀ * (
            1 +
            ( 3 // 4  ) * J₂  / p₀² * β * (2 - 3sin_i₀²) +
            ( 3 // 128) * J₂² / p₀⁴ * β * (120 + 64β - 40β² + (-240 - 192β + 40β²) * sin_i₀² + (105 + 144β + 25β²) * sin_i₀⁴) -
            (45 // 128) * J₄  / p₀⁴ * β * e₀² * (-8 + 40sin_i₀² - 35sin_i₀⁴)
        )

        ∂M_∂t = n̄
        ∂ω_∂t = ( 3 // 4  ) * n̄  * J₂  / p₀² * (4 - 5sin_i₀²) +
                ( 3 // 128) * n̄  * J₂² / p₀⁴ * (384 + 96e₀² - 384β + (-824 - 116e₀² + 1056β) * sin_i₀² + (430 - 5e₀² - 720β) * sin_i₀⁴) -
                (15 // 16 ) * n₀ * J₂² / p₀⁴ * e₀² * cos_i₀⁴ -
                (15 // 128) * n₀ * J₄  / p₀⁴ * (64 + 72e₀² - (248 + 252e₀²) * sin_i₀² + (196 + 189e₀²) * sin_i₀⁴)

        # Angular velocity.
        ang_vel = ∂M_∂t + ∂ω_∂t

        return ang_vel
    else
        throw(ArgumentError("The perturbation parameter :$perturbation is invalid."))
    end
end

function orbital_angular_velocity(orb::Orbit; kwargs...)
    # Convert first to Keplerian elements.
    k = convert(KeplerianElements, orb)
    return orbital_angular_velocity(k.a, k.e, k.i; kwargs...)
end

"""
    orbital_angular_velocity_to_semimajor_axis(angvel::Number, e::Number, i::Number; kwargs...) -> T, Bool

Compute the semi-major axis [m] that will provide an angular velocity `angvel` [rad / s] in
an orbit with eccentricity `e` and inclination `i` [rad].

Notice that the angular velocity `angvel` is related to the nodal period, *i.e.* the time
between two consecutive passages by the ascending node.

!!! note
    The output type `T` in the first signature is obtained by promoting the inputs to a
    float type.

# Keywords

- `max_iterations::Int`: Maximum number of iterations allowed in the Newton-Raphson
    algorithm. (**Default** = 20)
- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used.
    (**Default**: `:J2`)
- `tolerance::Union{Nothing, Number}`: Residue tolerances to verify if the numerical method
    has converged. If it is `nothing`, `√eps(T)` will be used, where `T` is the internal
    type for the computations. Notice that the residue function unit is [deg / min].
    (**Default** = nothing)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default** = `GM_EARTH`)
- `J2::Number`: J₂ perturbation term. (**Default** = EGM08_J2)
- `J4::Number`: J₄ perturbation term. (**Default** = EGM08_J4)
- `R0::Number`: Earth's equatorial radius [m]. (**Default** = EARTH_EQUATORIAL_RADIUS)

# Returns

- `T`: Semi-major axis [m].
- `Bool`: `true` if the numerical method converged, `false` otherwise.

# Perturbations

The keyword argument `perturbation` can be used to select the perturbation terms that will
be considered in the computation. The possible values are:

- `:J0`: Consider a Keplerian orbit.
- `:J2`: Consider the perturbation terms up to J2.
- `:J4`: Consider the perturbation terms J2, J4, and J2².

If `perturbation` is omitted, it defaults to `:J2`.
"""
function orbital_angular_velocity_to_semimajor_axis(
    angvel::T1,
    e::T2,
    i::T3;
    max_iterations::Int = 20,
    perturbation::Symbol = :J2,
    tolerance::Union{Nothing, Number} = nothing,
    # Constants.
    J2::Number = EGM08_J2,
    J4::Number = EGM08_J4,
    m0::Number = GM_EARTH,
    R0::Number = EARTH_EQUATORIAL_RADIUS
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    T = float(promote_type(T1, T2, T3))

    R₀  = T(R0)
    μ   = T(m0)
    J₂  = T(J2)
    J₄  = T(J4)
    tol = isnothing(tolerance) ? √eps(T) : T(tolerance)

    if perturbation == :J0

        a = (μ / T(angvel)^2)^(1 // 3)
        return a, true

    elseif perturbation == :J2
        rs_to_dm = T(60 * 180 / π)

        # Convert the inputs to the correct type.
        e₀  = T(e)
        i₀  = T(i)
        ω_d = T(angvel) * rs_to_dm

        # Auxiliary variables.
        β² = (1 - T(e₀)^2)
        β  = √β²
        β³ = β² * β
        β⁴ = β² * β²

        sin_i₀ = sin(T(i₀))
        sin_i₀² = sin_i₀^2

        k₁ = (3 // 4) * J₂ * (2 - 3sin_i₀²) / β³
        k₂ = (3 // 4) * J₂ * (4 - 5sin_i₀²) / β⁴
        k₃ = √(μ / R₀^3) * rs_to_dm

        # Newton-Raphson Algorithm
        # ==================================================================================

        # We defined the orbit angular velocity here based on the nodal period, i.e., the
        # time it takes for the satellite to cross the ascending node two consecutive times.
        # Hence, we can compute it by:
        #
        #             ∂M         ∂ω
        #   angvel = ──── (a) + ──── (a) .
        #             ∂t         ∂t
        #
        # The expressions for those time-derivatives were obtained from the J2 orbit
        # propagator of SatelliteToolboxPropagators.jl package.
        #
        # Since we cannot analytical isolate `a`, we will use a Newton-Raphson algorithm to
        # find the semi-major axis `a` that provides the desired angular velocity.

        # Initial guess based on the unperturbed model. Notice that we will estimate
        # `1 / √(a / R₀)`.
        isqrt_ā = √(R₀ * ((ω_d / rs_to_dm)^2 / T(μ))^(1 // 3))

        # By setting the initial values of `f₁` to `10tol`, we assure that the loop will be
        # executed at least one time.
        f₁ = 10tol

        # Loop.
        it = 1
        converged = true

        while abs(f₁) > tol
            isqrt_ā²  = isqrt_ā   * isqrt_ā
            isqrt_ā³  = isqrt_ā²  * isqrt_ā
            isqrt_ā⁶  = isqrt_ā³  * isqrt_ā³
            isqrt_ā⁷  = isqrt_ā⁶  * isqrt_ā
            isqrt_ā¹⁰ = isqrt_ā⁷  * isqrt_ā³
            isqrt_ā¹¹ = isqrt_ā¹⁰ * isqrt_ā

            # Compute the residue and the derivative.
            f₁ = ω_d - k₃ * (isqrt_ā³ + (k₁ + k₂) * isqrt_ā⁷ + k₁ * k₂ * isqrt_ā¹¹)

            @debug """
            Iteration #$it
              Estimation :
                a  = $(R₀ / (isqrt_ā * isqrt_ā) / 1000) km
              Residue :
                f₁ = $(f₁) ° / min
            """

            # Compute the function derivative.
            ∂f₁_∂isqrt_ā = - k₃ * (3isqrt_ā² + 7 * (k₁ + k₂) * isqrt_ā⁶ + 11 * k₁ * k₂ * isqrt_ā¹⁰)

            # Compute the new estimate.
            isqrt_ā = isqrt_ā - f₁ / ∂f₁_∂isqrt_ā

            # If the maximum number of iterations allowed has been reached, indicate that
            # the solution did not converged and exit loop.
            if (it >= max_iterations)
                converged = false
                break
            end

            it += 1
        end

        # Convert `isqrt_ā` to semi-major axis.
        a = R₀ / isqrt_ā^2

        return a, converged

    elseif perturbation == :J4
        rs_to_dm = T(60 * 180 / π)

        # Convert the inputs to the correct type.
        e₀  = T(e)
        i₀  = T(i)
        ω_d = T(angvel) * rs_to_dm

        # Auxiliary variables.
        e₀² = e₀^2
        J₂² = J₂^2
        β²  = (1 - e₀²)
        β   = √β²
        β³  = β² * β
        β⁴  = β² * β²
        β⁷  = β⁴ * β³
        β⁸  = β⁴ * β⁴

        sin_i₀, cos_i₀ = sincos(i₀)

        sin_i₀² = sin_i₀^2
        sin_i₀⁴ = sin_i₀^4
        cos_i₀⁴ = cos_i₀^4

        k₁ = +( 3 // 4  ) * J₂  / β³ * (2 - 3sin_i₀²)
        k₂ = +( 3 // 128) * J₂² / β⁷ * (120 + 64β - 40β² + (-240 - 192β + 40β²) * sin_i₀² + (105 + 144β + 25β²) * sin_i₀⁴)
        k₃ = -(45 // 128) * J₄  / β⁷ * e₀² * (-8 + 40sin_i₀² - 35sin_i₀⁴)
        k₄ = +( 3 // 4  ) * J₂  / β⁴ * (4 - 5sin_i₀²)
        k₅ = +( 3 // 128) * J₂² / β⁸ * (384 + 96e₀² - 384β + (-824 - 116e₀² + 1056β) * sin_i₀² + (430 - 5e₀² - 720β) * sin_i₀⁴)
        k₆ = -(15 // 16 ) * J₂² / β⁸ * e₀² * cos_i₀⁴
        k₇ = -(15 // 128) * J₄  / β⁸ * (64 + 72e₀² - (248 + 252e₀²) * sin_i₀² + (196 + 189e₀²) * sin_i₀⁴)
        k₈ = √(μ / R₀^3) * rs_to_dm

        # Newton-Raphson Algorithm
        # ==================================================================================

        # We defined the orbit angular velocity here based on the nodal period, i.e., the
        # time it takes for the satellite to cross the ascending node two consecutive times.
        # Hence, we can compute it by:
        #
        #             ∂M         ∂ω
        #   angvel = ──── (a) + ──── (a) .
        #             ∂t         ∂t
        #
        # The expressions for those time-derivatives were obtained from the J2 orbit
        # propagator of SatelliteToolboxPropagators.jl package.
        #
        # Since we cannot analytical isolate `a`, we will use a Newton-Raphson algorithm to
        # find the semi-major axis `a` that provides the desired angular velocity.

        # Initial guess based on the unperturbed model. Notice that we will estimate
        # `1 / √(a / R₀)`.
        isqrt_ā = √(R₀ * ((ω_d / rs_to_dm)^2 / T(μ))^(1 // 3))

        # By setting the initial values of `f₁` to `10tol`, we assure that the loop will be
        # executed at least one time.
        f₁ = 10tol

        # Loop.
        it = 1
        converged = true

        while abs(f₁) > tol
            isqrt_ā²  = isqrt_ā   * isqrt_ā
            isqrt_ā³  = isqrt_ā²  * isqrt_ā
            isqrt_ā⁶  = isqrt_ā³  * isqrt_ā³
            isqrt_ā⁷  = isqrt_ā⁶  * isqrt_ā
            isqrt_ā¹⁰ = isqrt_ā⁷  * isqrt_ā³
            isqrt_ā¹¹ = isqrt_ā¹⁰ * isqrt_ā
            isqrt_ā¹⁴ = isqrt_ā⁷  * isqrt_ā⁷
            isqrt_ā¹⁵ = isqrt_ā¹⁴ * isqrt_ā
            isqrt_ā¹⁸ = isqrt_ā¹¹ * isqrt_ā⁷
            isqrt_ā¹⁹ = isqrt_ā¹⁸ * isqrt_ā

            # Compute the residue and the derivative.
            f₁ = ω_d - k₈ * (
                isqrt_ā³  +
                isqrt_ā⁷  * (k₁ + k₄) +
                isqrt_ā¹¹ * (k₁ * k₄ + k₂ + k₃ + k₅ + k₆ + k₇) +
                isqrt_ā¹⁵ * (k₄ * (k₂ + k₃) + k₁ * k₅) +
                isqrt_ā¹⁹ * (k₂ + k₃) * k₅
            )

            @debug """
            Iteration #$it
              Estimation :
                a  = $(R₀ / (isqrt_ā * isqrt_ā) / 1000) km
              Residue :
                f₁ = $(f₁) ° / min
            """

            # Compute the function derivative.
            ∂f₁_∂isqrt_ā = - k₈ * (
                3  * isqrt_ā²  +
                7  * isqrt_ā⁶  * (k₁ + k₄) +
                11 * isqrt_ā¹⁰ * (k₁ * k₄ + k₂ + k₃ + k₅ + k₆ + k₇) +
                15 * isqrt_ā¹⁴ * (k₄ * (k₂ + k₃) + k₁ * k₅) +
                19 * isqrt_ā¹⁸ * (k₂ + k₃) * k₅
            )

            # Compute the new estimate.
            isqrt_ā = isqrt_ā - f₁ / ∂f₁_∂isqrt_ā

            # If the maximum number of iterations allowed has been reached, indicate that
            # the solution did not converged and exit loop.
            if (it >= max_iterations)
                converged = false
                break
            end

            it += 1
        end

        # Convert `isqrt_ā` to semi-major axis.
        a = R₀ / isqrt_ā^2

        return a, converged
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

# Keywords

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used.
    (**Default**: `:J2`)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default** = `GM_EARTH`)
- `J2::Number`: J₂ perturbation term. (**Default** = EGM08_J2)
- `J4::Number`: J₄ perturbation term. (**Default** = EGM08_J4)
- `R0::Number`: Earth's equatorial radius [m]. (**Default** = EARTH_EQUATORIAL_RADIUS)

# Perturbations

The keyword argument `perturbation` can be used to select the perturbation terms that will
be considered in the computation. The possible values are:

- `:J0`: Consider a Keplerian orbit.
- `:J2`: Consider the perturbation terms up to J2.
- `:J4`: Consider the perturbation terms J2, J4, and J2².

If `perturbation` is omitted, it defaults to `:J2`.
"""
function orbital_period(a::Number, e::Number, i::Number; kwargs...)
    n = orbital_angular_velocity(a, e, i; kwargs...)
    T = eltype(n)
    return T(2π) / n
end

function orbital_period(orb::Orbit; kwargs...)
    # Convert first to Keplerian elements.
    k = convert(KeplerianElements, orb)
    return orbital_period(k.a, k.e, k.i; kwargs...)
end

"""
    raan_time_derivative(a::Number, e::Number, i::Number; kwargs...) -> T

Compute the time derivative of the right ascension of the ascending node (RAAN) [rad / s] in
an orbit with semi-major axis `a` [m], eccentricity `e`, and inclination `i` [rad]. The
orbit can also be specified by `orb` (see `Orbit`).

!!! note
    The output type `T` in the first signature is obtained by promoting the inputs to a
    float type.

# Keywords

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used.
    (**Default**: `:J2`)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default** = `GM_EARTH`)
- `J2::Number`: J₂ perturbation term. (**Default** = EGM08_J2)
- `J4::Number`: J₄ perturbation term. (**Default** = EGM08_J4)
- `R0::Number`: Earth's equatorial radius [m]. (**Default** = EARTH_EQUATORIAL_RADIUS)

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
    perturbation::Symbol = :J2,
    # Constants.
    J2::Number = EGM08_J2,
    J4::Number = EGM08_J4,
    m0::Number = GM_EARTH,
    R0::Number = EARTH_EQUATORIAL_RADIUS
) where {T1<:Number, T2<:Number, T3<:Number}
    T = float(promote_type(T1, T2, T3))

    R₀ = T(R0)
    μ  = T(m0)
    J₂ = T(J2)
    J₄ = T(J4)

    # Perturbation computed using a Keplerian orbit.
    if perturbation == :J0
        return zero(T)

    # Perturbation computed using perturbations terms up to J2.
    elseif perturbation == :J2
        # Convert the inputs to the correct type.
        a₀  = T(a)
        e₀  = T(e)
        i₀  = T(i)

        # Auxiliary variables.
        μm  = √(μ / R₀^3)
        al₀ = a₀ / R₀
        e₀² = e₀^2
        n₀  = μm / √(al₀^3)
        p₀  = al₀ * (1 - e₀²)
        p₀² = p₀^2

        sin_i₀, cos_i₀ = sincos(T(i₀))
        sin_i₀² = sin_i₀^2
        β²      = 1 - e₀²
        β       = √β²

        # Perturbed orbit mean motion.
        n̄ = n₀ * (1 + (3 // 4) * J₂ / p₀² * β * (2 - 3sin_i₀²))

        # First-order time-derivative of the RAAN [rad / s].
        ∂Ω = -(3 // 2) * n̄ * J₂ / p₀² * cos_i₀

        return ∂Ω

    # Perturbation computed using perturbation terms J₂, J₂², and J₄.
    elseif perturbation == :J4
        # Convert the inputs to the correct type.
        a₀  = T(a)
        e₀  = T(e)
        i₀  = T(i)

        # Auxiliary variables.
        μm  = √(μ / R₀^3)
        al₀ = a₀ / R₀
        e₀² = e₀^2
        p₀  = al₀ * (1 - e₀²)
        p₀² = p₀^2
        p₀⁴ = p₀^4
        n₀  = μm / √(al₀^3)
        J₂² = J₂^2

        sin_i₀, cos_i₀ = sincos(T(i₀))

        sin_i₀² = sin_i₀^2
        sin_i₀⁴ = sin_i₀^4
        β²      = (1 - e₀²)
        β       = √β²

        sin_i, cos_i = sincos(i)
        sin_i² = sin_i^2

        # Perturbed mean motion.
        kn₂  = J₂  / p₀² * β
        kn₂₂ = J₂² / p₀⁴ * β
        kn₄  = J₄  / p₀⁴ * β

        n̄ = n₀ * (
            1 +
            ( 3 // 4  ) * kn₂  * (2 - 3sin_i₀²) +
            ( 3 // 128) * kn₂₂ * (120 + 64β - 40β² + (-240 - 192β + 40β²) * sin_i₀² + (105 + 144β + 25β²) * sin_i₀⁴) -
            (45 // 128) * kn₄  * e₀² * (-8 + 40sin_i₀² - 35sin_i₀⁴)
        )

        # First-order time-derivative of the RAAN [rad / s].
        k̄₂  = n̄  * J₂  / p₀²
        k̄₂₂ = n̄  * J₂² / p₀⁴
        k₂₂ = n₀ * J₂² / p₀⁴
        k₄  = n₀ * J₄  / p₀⁴

        ∂Ω = -( 3 // 2 ) * k̄₂  * cos_i₀ +
              ( 3 // 32) * k̄₂₂ * cos_i₀ * (-36 -  4e₀² + 48β + (40 - 5e₀² - 72β) * sin_i₀²) +
              (15 // 32) * k₄  * cos_i₀ * (8 + 12e₀² - (14 + 21e₀²) * sin_i₀²)

        return ∂Ω
    else
        throw(ArgumentError("The perturbation parameter $perturbation is not defined."))
    end
end

function raan_time_derivative(orb::Orbit; kwargs...)
    # Convert first to Keplerian elements.
    k = convert(KeplerianElements, orb)
    return raan_time_derivative(k.a, k.e, k.i; kwargs...)
end
