## Description #############################################################################
#
#   Functions to compute general values related to the orbit.
#
## References ##############################################################################
#
# [1] Kozai, Y (1959). The Motion of a Close Earth Satellite. The Astronomical Journal,
#     v. 64, no. 1274, pp. 367 -- 377.
#
# [2] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorne, CA.
#
############################################################################################

export orbital_angular_velocity
export orbital_angular_velocity_to_semimajor_axis
export orbital_period
export raan_time_derivative

############################################################################################
#                                        Functions                                         #
############################################################################################

"""
    orbital_angular_velocity(a::Number, e::Number, i::Number; kwargs...) -> T
    orbital_angular_velocity(orb::Orbit; kwargs...) -> T

Compute the angular velocity [rad / s] of an object in an orbit with semi-major axis `a`
[m], eccentricity `e` [-], and inclination `i` [rad]. The orbit can also be specified by
`orb` (see `Orbit`). The inputs are validated and this function throws an error if they do
not describe a valid elliptical orbit.

The angular velocity is defined here based on the nodal period, *i.e.* the time between two
consecutive passages by the ascending node. Hence, it is the sum of the perturbed mean
motion and the argument of perigee time derivative.

!!! note

    The output type `T` in the first signature is obtained by promoting the inputs to a
    float type.

# Keywords

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used.
    (**Default**: `:J2`)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default**: `GM_EARTH`)
- `J2::Number`: J₂ perturbation term.
    (**Default**: `EGM_2008_J2`)
- `J4::Number`: J₄ perturbation term.
    (**Default**: `EGM_2008_J4`)
- `R0::Number`: Earth's equatorial radius [m].
    (**Default**: `EARTH_EQUATORIAL_RADIUS`)

# Perturbations

The keyword argument `perturbation` can be used to select the perturbation terms that will
be considered in the computation. The possible values are:

- `:J0`: Consider a Keplerian orbit.
- `:J2`: Consider the perturbation terms up to J2.
- `:J4`: Consider the perturbation terms J2, J4, and J2².

If `perturbation` is omitted, it defaults to `:J2`.

# Extended help

## Throws

- `ArgumentError`: If `perturbation` is not `:J0`, `:J2`, or `:J4`.
- `ArgumentError`: If the eccentricity `e` is not in the interval [0, 1).
- `ArgumentError`: If the perigee radius `a * (1 - e)` is not positive.
"""
function orbital_angular_velocity(
    a::T1,
    e::T2,
    i::T3;
    perturbation::Symbol = :J2,
    # Constants.
    J2::Number = EGM_2008_J2,
    J4::Number = EGM_2008_J4,
    m0::Number = GM_EARTH,
    R0::Number = EARTH_EQUATORIAL_RADIUS
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    T = float(promote_type(T1, T2, T3))

    n̄, ∂ω, _ = _secular_rates(perturbation, T(a), T(e), T(i), T(m0), T(R0), T(J2), T(J4))

    return n̄ + ∂ω
end

# NOTE: The method for `KeplerianElements` reads the fields directly, avoiding the
# conversion to the true anomaly (which solves Kepler's equation) that `convert` performs
# for the mean anomaly representation.
function orbital_angular_velocity(orb::KeplerianElements; kwargs...)
    return orbital_angular_velocity(
        orb.semi_major_axis,
        orb.eccentricity,
        orb.inclination;
        kwargs...
    )
end

function orbital_angular_velocity(orb::Orbit; kwargs...)
    # Convert first to Keplerian elements.
    return orbital_angular_velocity(convert(KeplerianElements, orb); kwargs...)
end

"""
    orbital_angular_velocity_to_semimajor_axis(
        angvel::Number,
        e::Number,
        i::Number;
        kwargs...
    ) -> T, Bool

Compute the semi-major axis [m] that will provide an angular velocity `angvel` [rad / s] in
an orbit with eccentricity `e` [-] and inclination `i` [rad]. The inputs are validated and
this function throws an error if they do not describe a valid elliptical orbit.

Notice that the angular velocity `angvel` is related to the nodal period, *i.e.* the time
between two consecutive passages by the ascending node.

!!! note

    The output type `T` is obtained by promoting the inputs to a float type.

# Keywords

- `max_iterations::Int`: Maximum number of iterations allowed in the Newton-Raphson
    algorithm.
    (**Default**: `20`)
- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used.
    (**Default**: `:J2`)
- `tolerance::Union{Nothing, Number}`: Relative tolerance to verify if the numerical method
    has converged. The algorithm converges when the absolute difference between the angular
    velocity of the estimated orbit and `angvel` is lower than or equal to
    `tolerance * angvel`. It must be greater than 0, otherwise this function throws an
    `ArgumentError`. If it is `nothing`, `√eps(T)` will be used, where `T` is the internal
    type for the computations.
    (**Default**: `nothing`)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default**: `GM_EARTH`)
- `J2::Number`: J₂ perturbation term.
    (**Default**: `EGM_2008_J2`)
- `J4::Number`: J₄ perturbation term.
    (**Default**: `EGM_2008_J4`)
- `R0::Number`: Earth's equatorial radius [m].
    (**Default**: `EARTH_EQUATORIAL_RADIUS`)

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

# Extended help

## Throws

- `ArgumentError`: If `perturbation` is not `:J0`, `:J2`, or `:J4`.
- `ArgumentError`: If the angular velocity `angvel` is not positive.
- `ArgumentError`: If the eccentricity `e` is not in the interval [0, 1).
- `ArgumentError`: If the keyword `tolerance` is not `nothing` and not positive.
"""
function orbital_angular_velocity_to_semimajor_axis(
    angvel::T1,
    e::T2,
    i::T3;
    max_iterations::Int = 20,
    perturbation::Symbol = :J2,
    tolerance::Union{Nothing, Number} = nothing,
    # Constants.
    J2::Number = EGM_2008_J2,
    J4::Number = EGM_2008_J4,
    m0::Number = GM_EARTH,
    R0::Number = EARTH_EQUATORIAL_RADIUS
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    T = float(promote_type(T1, T2, T3))

    angvel <= 0 && throw(ArgumentError("The angular velocity must be greater than 0."))

    !(0 <= e < 1) && throw(
        ArgumentError("The eccentricity must be in the interval [0, 1), but it is $e.")
    )

    !isnothing(tolerance) && (tolerance <= 0) && throw(
        ArgumentError("The keyword `tolerance` must be greater than 0.")
    )

    # Convert the inputs to the correct type.
    R₀ = T(R0)
    μ  = T(m0)
    J₂ = T(J2)
    J₄ = T(J4)
    e₀ = T(e)
    i₀ = T(i)
    ω  = T(angvel)

    tol = isnothing(tolerance) ? √eps(T) : T(tolerance)

    # Semi-major axis of the unperturbed orbit with the desired angular velocity [m]. It is
    # the solution for the Keplerian orbit and the initial guess for the perturbed models.
    a₀ = cbrt(μ / ω^2)

    perturbation == :J0 && return a₀, true

    # == Newton-Raphson Algorithm ==========================================================

    # We defined the orbit angular velocity here based on the nodal period, i.e., the time
    # it takes for the satellite to cross the ascending node two consecutive times. Hence,
    # we can compute it by:
    #
    #             ∂M         ∂ω
    #   angvel = ──── (a) + ──── (a) .
    #             ∂t         ∂t
    #
    # Since we cannot analytically isolate `a`, we will use a Newton-Raphson algorithm to
    # find the semi-major axis `a` that provides the desired angular velocity.
    #
    # If we define `x = 1 / √(a / R₀)` and `y = x⁴`, the angular velocity is the polynomial:
    #
    #   angvel(x) = √(μ / R₀³) ⋅ x³ ⋅ g(y),    g(y) = 1 + c₁ y + c₂ y² + c₃ y³ + c₄ y⁴ ,
    #
    # where the coefficients `cₖ` depend only on the perturbation model, the eccentricity,
    # and the inclination (see `_angular_velocity_polynomial_coefficients`). We estimate `x`
    # and evaluate `g` and its derivative using the Horner's method.
    c₁, c₂, c₃, c₄ = _angular_velocity_polynomial_coefficients(perturbation, e₀, i₀, J₂, J₄)

    # Normalized unperturbed mean motion [rad / s], i.e. the mean motion at `a = R₀`.
    k = √(μ / R₀^3)

    # Initial guess based on the unperturbed model.
    x = √(R₀ / a₀)

    # Newton-Raphson loop. The residue is evaluated at the top of the loop so that the
    # `converged` flag always describes the returned estimate.
    it = 0
    converged = false

    while true
        x² = x * x
        x³ = x² * x
        y  = x² * x²

        # Evaluate the polynomial and its derivative with respect to `y`.
        g  = 1 + y * (c₁ + y * (c₂ + y * (c₃ + y * c₄)))
        ∂g = c₁ + y * (2c₂ + y * (3c₃ + 4c₄ * y))

        # Compute the residue at the current estimate [rad / s].
        f = k * x³ * g - ω

        @debug """
        Iteration #$it
          Estimation :
            a = $(R₀ / x² / 1000) km
          Residue :
            f = $(f) rad / s
        """

        # If the residue at the current estimate is within the relative tolerance, indicate
        # that the solution converged and exit the loop.
        if abs(f) <= tol * ω
            converged = true
            break
        end

        # If the maximum number of iterations allowed has been reached, indicate that the
        # solution did not converge and exit the loop.
        (it >= max_iterations) && break

        # Compute the residue derivative with respect to `x`.
        ∂f = k * x² * (3g + 4y * ∂g)

        # Compute the new estimate.
        x -= f / ∂f

        it += 1
    end

    # Convert `x` to semi-major axis.
    a = R₀ / x^2

    return a, converged
end

"""
    orbital_period(a::Number, e::Number, i::Number; kwargs...) -> T
    orbital_period(orb::Orbit; kwargs...) -> T

Compute the orbital period [s] of an object in an orbit with semi-major axis `a` [m],
eccentricity `e` [-], and inclination `i` [rad]. The orbit can also be specified by `orb`
(see `Orbit`). The inputs are validated and this function throws an error if they do not
describe a valid elliptical orbit.

The period is defined here based on the nodal period, *i.e.* the time between two
consecutive passages by the ascending node.

!!! note

    The output type `T` in the first signature is obtained by promoting the inputs to a
    float type.

# Keywords

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used.
    (**Default**: `:J2`)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default**: `GM_EARTH`)
- `J2::Number`: J₂ perturbation term.
    (**Default**: `EGM_2008_J2`)
- `J4::Number`: J₄ perturbation term.
    (**Default**: `EGM_2008_J4`)
- `R0::Number`: Earth's equatorial radius [m].
    (**Default**: `EARTH_EQUATORIAL_RADIUS`)

# Perturbations

The keyword argument `perturbation` can be used to select the perturbation terms that will
be considered in the computation. The possible values are:

- `:J0`: Consider a Keplerian orbit.
- `:J2`: Consider the perturbation terms up to J2.
- `:J4`: Consider the perturbation terms J2, J4, and J2².

If `perturbation` is omitted, it defaults to `:J2`.

# Extended help

## Throws

- `ArgumentError`: If `perturbation` is not `:J0`, `:J2`, or `:J4`.
- `ArgumentError`: If the eccentricity `e` is not in the interval [0, 1).
- `ArgumentError`: If the perigee radius `a * (1 - e)` is not positive.
"""
function orbital_period(a::Number, e::Number, i::Number; kwargs...)
    n = orbital_angular_velocity(a, e, i; kwargs...)
    T = typeof(n)
    return T(2π) / n
end

# NOTE: The method for `KeplerianElements` reads the fields directly, avoiding the
# conversion to the true anomaly (which solves Kepler's equation) that `convert` performs
# for the mean anomaly representation.
function orbital_period(orb::KeplerianElements; kwargs...)
    return orbital_period(
        orb.semi_major_axis,
        orb.eccentricity,
        orb.inclination;
        kwargs...
    )
end

function orbital_period(orb::Orbit; kwargs...)
    # Convert first to Keplerian elements.
    return orbital_period(convert(KeplerianElements, orb); kwargs...)
end

"""
    raan_time_derivative(a::Number, e::Number, i::Number; kwargs...) -> T
    raan_time_derivative(orb::Orbit; kwargs...) -> T

Compute the time derivative of the right ascension of the ascending node (RAAN) [rad / s] in
an orbit with semi-major axis `a` [m], eccentricity `e` [-], and inclination `i` [rad]. The
orbit can also be specified by `orb` (see `Orbit`). The inputs are validated and this
function throws an error if they do not describe a valid elliptical orbit.

!!! note

    The output type `T` in the first signature is obtained by promoting the inputs to a
    float type.

# Keywords

- `perturbation::Symbol`: Symbol to select the perturbation terms that will be used.
    (**Default**: `:J2`)
- `m0::Number`: Standard gravitational parameter for Earth [m³ / s²].
    (**Default**: `GM_EARTH`)
- `J2::Number`: J₂ perturbation term.
    (**Default**: `EGM_2008_J2`)
- `J4::Number`: J₄ perturbation term.
    (**Default**: `EGM_2008_J4`)
- `R0::Number`: Earth's equatorial radius [m].
    (**Default**: `EARTH_EQUATORIAL_RADIUS`)

# Perturbations

The keyword argument `perturbation` can be used to select the perturbation terms that will
be considered in the computation. The possible values are:

- `:J0`: Consider a Keplerian orbit.
- `:J2`: Consider the perturbation terms up to J2.
- `:J4`: Consider the perturbation terms J2, J4, and J2².

If `perturbation` is omitted, it defaults to `:J2`.

# Extended help

## Throws

- `ArgumentError`: If `perturbation` is not `:J0`, `:J2`, or `:J4`.
- `ArgumentError`: If the eccentricity `e` is not in the interval [0, 1).
- `ArgumentError`: If the perigee radius `a * (1 - e)` is not positive.
"""
function raan_time_derivative(
    a::T1,
    e::T2,
    i::T3;
    perturbation::Symbol = :J2,
    # Constants.
    J2::Number = EGM_2008_J2,
    J4::Number = EGM_2008_J4,
    m0::Number = GM_EARTH,
    R0::Number = EARTH_EQUATORIAL_RADIUS
) where {T1 <: Number, T2 <: Number, T3 <: Number}
    T = float(promote_type(T1, T2, T3))

    _, _, ∂Ω = _secular_rates(perturbation, T(a), T(e), T(i), T(m0), T(R0), T(J2), T(J4))

    return ∂Ω
end

# NOTE: The method for `KeplerianElements` reads the fields directly, avoiding the
# conversion to the true anomaly (which solves Kepler's equation) that `convert` performs
# for the mean anomaly representation.
function raan_time_derivative(orb::KeplerianElements; kwargs...)
    return raan_time_derivative(
        orb.semi_major_axis,
        orb.eccentricity,
        orb.inclination;
        kwargs...
    )
end

function raan_time_derivative(orb::Orbit; kwargs...)
    # Convert first to Keplerian elements.
    return raan_time_derivative(convert(KeplerianElements, orb); kwargs...)
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

# The secular theory implemented here follows the J2 and J4 orbit propagators of
# SatelliteToolboxPropagators.jl, which are based on [1, 2]. Given the unperturbed mean
# motion `n₀` and the semi-latus rectum `p₀` normalized by the Earth's equatorial radius,
# the perturbed mean motion `n̄` and the first-order time derivatives of the argument of
# perigee `∂ω` and of the RAAN `∂Ω` are:
#
#   n̄  = n₀ ⋅ (1 + A / p₀² + B / p₀⁴) ,
#   ∂ω = n̄ ⋅ (C / p₀² + D / p₀⁴) + n₀ ⋅ E / p₀⁴ ,
#   ∂Ω = n̄ ⋅ (F / p₀² + G / p₀⁴) + n₀ ⋅ H / p₀⁴ ,
#
# where the coefficients `A` to `H` depend only on the perturbation model, the eccentricity,
# and the inclination. They are computed by `_secular_coefficients`, which is the single
# place in this package where the perturbation theory is written down.

"""
    _secular_coefficients(
        perturbation::Symbol,
        e::T,
        i::T,
        J₂::T,
        J₄::T
    ) where {T <: Number} -> NTuple{8, T}

    _secular_coefficients(
        ::Val{:J0},
        e::T,
        i::T,
        J₂::T,
        J₄::T
    ) where {T <: Number} -> NTuple{8, T}

    _secular_coefficients(
        ::Val{:J2},
        e::T,
        i::T,
        J₂::T,
        J₄::T
    ) where {T <: Number} -> NTuple{8, T}

    _secular_coefficients(
        ::Val{:J4},
        e::T,
        i::T,
        J₂::T,
        J₄::T
    ) where {T <: Number} -> NTuple{8, T}

Compute the coefficients `(A, B, C, D, E, F, G, H)` of the secular theory selected by
`perturbation` (`:J0`, `:J2`, or `:J4`) for an orbit with eccentricity `e` [-] and
inclination `i` [rad], using the zonal harmonics `J₂` and `J₄`.

The coefficients relate the unperturbed mean motion `n₀` and the normalized semi-latus
rectum `p₀` [er] to the perturbed mean motion `n̄`, the argument of perigee time derivative
`∂ω`, and the RAAN time derivative `∂Ω` as follows:

    n̄  = n₀ ⋅ (1 + A / p₀² + B / p₀⁴)
    ∂ω = n̄ ⋅ (C / p₀² + D / p₀⁴) + n₀ ⋅ E / p₀⁴
    ∂Ω = n̄ ⋅ (F / p₀² + G / p₀⁴) + n₀ ⋅ H / p₀⁴

# Extended help

## Throws

- `ArgumentError`: If `perturbation` is not `:J0`, `:J2`, or `:J4`.
"""
function _secular_coefficients(
    perturbation::Symbol,
    e::T,
    i::T,
    J₂::T,
    J₄::T
) where {T <: Number}
    perturbation == :J0 && return _secular_coefficients(Val(:J0), e, i, J₂, J₄)
    perturbation == :J2 && return _secular_coefficients(Val(:J2), e, i, J₂, J₄)
    perturbation == :J4 && return _secular_coefficients(Val(:J4), e, i, J₂, J₄)
    return throw(ArgumentError("The perturbation parameter :$perturbation is invalid."))
end

function _secular_coefficients(::Val{:J0}, e::T, i::T, J₂::T, J₄::T) where {T <: Number}
    z = zero(T)

    return (z, z, z, z, z, z, z, z)
end

function _secular_coefficients(::Val{:J2}, e::T, i::T, J₂::T, J₄::T) where {T <: Number}
    sin_i, cos_i = sincos(i)
    sin_i² = sin_i^2
    β = √(1 - e^2)

    # First-order secular terms, which depend only on J₂ [1].
    A = +(3//4) * J₂ * β * (2 - 3sin_i²)
    C = +(3//4) * J₂ * (4 - 5sin_i²)
    F = -(3//2) * J₂ * cos_i

    z = zero(T)

    return (A, z, C, z, z, F, z, z)
end

function _secular_coefficients(::Val{:J4}, e::T, i::T, J₂::T, J₄::T) where {T <: Number}
    sin_i, cos_i = sincos(i)

    sin_i² = sin_i^2
    sin_i⁴ = sin_i^4
    cos_i⁴ = cos_i^4
    e²     = e^2
    β²     = 1 - e²
    β      = √β²
    J₂²    = J₂^2

    # First-order secular terms, which depend only on J₂ [1].
    A = +(3//4) * J₂ * β * (2 - 3sin_i²)
    C = +(3//4) * J₂ * (4 - 5sin_i²)
    F = -(3//2) * J₂ * cos_i

    # Second-order secular terms, which depend on J₂² and J₄ [1].
    B = +(3//128) * J₂² * β * (
            120 + 64β - 40β² +
            (-240 - 192β + 40β²) * sin_i² +
            (105 + 144β + 25β²) * sin_i⁴
        ) -
        (45//128) * J₄ * β * e² * (-8 + 40sin_i² - 35sin_i⁴)

    D = +(3//128) * J₂² * (
            384 + 96e² - 384β +
            (-824 - 116e² + 1056β) * sin_i² +
            (430 - 5e² - 720β) * sin_i⁴
        )

    E = -(15//16) * J₂² * e² * cos_i⁴ -
        (15//128) * J₄ * (
            64 + 72e² -
            (248 + 252e²) * sin_i² +
            (196 + 189e²) * sin_i⁴
        )

    G = +(3//32) * J₂² * cos_i * (-36 - 4e² + 48β + (40 - 5e² - 72β) * sin_i²)

    H = +(15//32) * J₄ * cos_i * (8 + 12e² - (14 + 21e²) * sin_i²)

    return (A, B, C, D, E, F, G, H)
end

"""
    _secular_rates(
        perturbation::Symbol,
        a::T,
        e::T,
        i::T,
        μ::T,
        R₀::T,
        J₂::T,
        J₄::T
    ) where {T <: Number} -> T, T, T

Compute the perturbed mean motion [rad / s], the argument of perigee time derivative
[rad / s], and the RAAN time derivative [rad / s] of an orbit with semi-major axis `a` [m],
eccentricity `e` [-], and inclination `i` [rad] using the secular theory selected by
`perturbation` (`:J0`, `:J2`, or `:J4`), the standard gravitational parameter `μ` [m³ / s²],
the Earth's equatorial radius `R₀` [m], and the zonal harmonics `J₂` and `J₄`.

# Returns

- `T`: Perturbed mean motion `n̄` [rad / s].
- `T`: Argument of perigee time derivative `∂ω` [rad / s].
- `T`: RAAN time derivative `∂Ω` [rad / s].

# Extended help

## Throws

- `ArgumentError`: If `perturbation` is not `:J0`, `:J2`, or `:J4`.
- `ArgumentError`: If the eccentricity `e` is not in the interval [0, 1).
- `ArgumentError`: If the perigee radius `a * (1 - e)` is not positive.
"""
function _secular_rates(
    perturbation::Symbol,
    a::T,
    e::T,
    i::T,
    μ::T,
    R₀::T,
    J₂::T,
    J₄::T
) where {T <: Number}
    # The theory is only valid for elliptical orbits. Without these checks, the user would
    # get a `DomainError` from an internal square root or silently wrong results.
    !(0 <= e < 1) && throw(
        ArgumentError("The eccentricity must be in the interval [0, 1), but it is $e.")
    )

    a * (1 - e) <= 0 && throw(
        ArgumentError(
            "The perigee radius must be positive, but the semi-major axis is $a m and " *
            "the eccentricity is $e."
        )
    )

    A, B, C, D, E, F, G, H = _secular_coefficients(perturbation, e, i, J₂, J₄)

    # Auxiliary variables.
    n₀   = √(μ / a^3)         # .......................... Unperturbed mean motion [rad / s]
    p₀   = a / R₀ * (1 - e^2) # .......................... Normalized semi-latus rectum [er]
    ip₀² = 1 / p₀^2           # ............................................. 1 / p₀² [er⁻²]
    ip₀⁴ = ip₀²^2             # ............................................. 1 / p₀⁴ [er⁻⁴]

    # Perturbed mean motion [rad / s].
    n̄ = n₀ * (1 + A * ip₀² + B * ip₀⁴)

    # First-order time derivatives of the argument of perigee and of the RAAN [rad / s].
    ∂ω = n̄ * (C * ip₀² + D * ip₀⁴) + n₀ * E * ip₀⁴
    ∂Ω = n̄ * (F * ip₀² + G * ip₀⁴) + n₀ * H * ip₀⁴

    return n̄, ∂ω, ∂Ω
end

"""
    _angular_velocity_polynomial_coefficients(
        perturbation::Symbol,
        e::T,
        i::T,
        J₂::T,
        J₄::T
    ) where {T <: Number} -> NTuple{4, T}

Compute the coefficients `(c₁, c₂, c₃, c₄)` of the polynomial that provides the orbital
angular velocity as a function of `x = 1 / √(a / R₀)`, where `a` is the semi-major axis and
`R₀` is the Earth's equatorial radius, for an orbit with eccentricity `e` [-] and
inclination `i` [rad] using the secular theory selected by `perturbation` (`:J0`, `:J2`, or
`:J4`) and the zonal harmonics `J₂` and `J₄`.

The angular velocity is given by:

    angvel(x) = √(μ / R₀³) ⋅ x³ ⋅ (1 + c₁ y + c₂ y² + c₃ y³ + c₄ y⁴),    y = x⁴ .

# Extended help

## Throws

- `ArgumentError`: If `perturbation` is not `:J0`, `:J2`, or `:J4`.
"""
function _angular_velocity_polynomial_coefficients(
    perturbation::Symbol,
    e::T,
    i::T,
    J₂::T,
    J₄::T
) where {T <: Number}
    A, B, C, D, E, _, _, _ = _secular_coefficients(perturbation, e, i, J₂, J₄)

    # Since `p₀ = ā β²`, where `ā = a / R₀` and `β² = 1 - e²`, we have `1 / p₀² = y / β⁴`
    # and `1 / p₀⁴ = y² / β⁸`. Hence, we fold the eccentricity factors into the secular
    # coefficients.
    iβ⁴ = 1 / (1 - e^2)^2
    iβ⁸ = iβ⁴^2

    Ā = A * iβ⁴
    B̄ = B * iβ⁸
    C̄ = C * iβ⁴
    D̄ = D * iβ⁸
    Ē = E * iβ⁸

    # The angular velocity is:
    #
    #   n̄ + ∂ω = n₀ ⋅ [(1 + Ā y + B̄ y²) ⋅ (1 + C̄ y + D̄ y²) + Ē y²] .
    #
    # Expanding the product, we obtain the polynomial coefficients.
    c₁ = Ā + C̄
    c₂ = Ā * C̄ + B̄ + D̄ + Ē
    c₃ = Ā * D̄ + B̄ * C̄
    c₄ = B̄ * D̄

    return (c₁, c₂, c₃, c₄)
end
