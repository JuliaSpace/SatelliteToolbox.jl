#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions to compute general values related to the orbit.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export @check_orbit
export angvel, angvel_to_a, dargp, draan, period

################################################################################
#                                    Macros
################################################################################

"""
    @check_orbit(a, e)

Verify if the orbit with semi-major axis `a` [m] and eccentricity `e` is valid.
This macro throws an exception if the orbit is not valid.

Return `true` is the orbit is valid, and `false` otherwise.

"""
macro check_orbit(a, e)
    quote
        # Check if the arguments are valid.
        if $(esc(a)) < 0
            false
        elseif !( 0 <= $(esc(e)) < 1 )
            false
        else
            # Compute the orbit perigee and check if it is inside Earth.
            perigee = $(esc(a))*(1-$(esc(e)))

            if perigee < R0
                false
            else
                true
            end
        end
    end
end

################################################################################
#                                  Functions
################################################################################

"""
    angvel(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    angvel(orb::Orbit, pert::Symbol = :J2)

Compute the angular velocity [rad/s] of an object in an orbit with semi-major
axis `a` [m], eccentricity `e`, and inclination `i` [rad], using the
perturbation terms specified by the symbol `pert`. The orbit can also be
specified by `orb` (see [`Orbit`](@ref)).

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.
* `:J4`: Consider the perturbation terms J2, J4, and J2².

If `pert` is omitted, then it defaults to `:J2`.

"""
@inline function angvel(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    # Unperturbed orbit period.
    n0 = sqrt(m0/Float64(a)^3)

    # Perturbation computed using a Keplerian orbit.
    if pert == :J0
        return n0
    # Perturbation computed using perturbations terms up to J2.
    elseif pert == :J2
        # Auxiliary variables.
        cos_i² = cos(i)^2
        aux   = 1-e^2

        # Semi-lactum rectum.
        p = a*aux

        # Orbit period considering the perturbations (up to J2).
        return n0 + 3R0^2*J2/(4p^2)*n0*(sqrt(aux)*(3cos_i²-1) + (5cos_i²-1))

    # Perturbation computed using perturbations terms J2, J4, and J2².
    elseif pert == :J4

        # Auxiliary variables
        e²     = e^2
        e⁴     = e^4
        sin_i  = sin(i)
        sin_i² = sin_i^2
        sin_i⁴ = sin_i^4
        aux    = (1-e²)
        saux   = sqrt(aux)
        p      = (a/R0)*aux
        p²     = p^2
        p⁴     = p^4

        # Notice that:
        #            .   .
        #   n = n₀ + ω + M₀
        #
        # in which the time-derivatives are computed as in [1, p. 692].

        δω   =  3/  4*n0*J2  /p²*(4 - 5sin_i²) +
                9/384*n0*J2^2/p⁴*(     56e² + (760 -  36e²)*sin_i² - (890 +  45e²)*sin_i⁴) -
               15/128*n0*J4  /p⁴*(64 + 72e² - (248 + 252e²)*sin_i² + (196 + 189e²)*sin_i⁴)

        δM_0 =  3/  4*n0*J2  /p²*saux*(2 - 3sin_i²) +
                3/512*n0*J2^2/p⁴/saux*(320e² - 280e⁴ + (1600 - 1568e² + 328e⁴)*sin_i² + (-2096 + 1072e² +  79e⁴)*sin_i⁴) -
               45/128*n0*J4  /p⁴*saux*e²*(-8 + 40sin_i - 35sin_i²)

        return n0 + δω + δM_0
    else
        throw(ArgumentError("The perturbation parameter $pert is not defined."))
    end
end

@inline function angvel(orb::Orbit, pert::Symbol = :J2)
    # Convert first to Keplerian elements.
    k = convert(KeplerianElements, orb)
    return angvel(get_a(k), get_e(k), get_i(k), pert)
end

"""
    angvel_to_a(n::Number, e::Number, i::Number, pert::Symbol = :J2; μ::Number = m0, max_iter::Int = 20, tol::Number = 1e-10)

Compute the semi-major axis that will provide an angular velocity `n` [rad/s] in
an orbit with eccentricity `e` and inclination `i` [rad], using the perturbation
terms specified by the symbol `pert`.

Notice that the angular velocity `n` is related to the nodal period, *i.e.* the
time between two consecutive passages by the ascending node.

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.
* `:J4`: Consider the perturbation terms J2, J4, and J2².

If `pert` is omitted, then it defaults to `:J2`.

# Keyword

* `μ`: Standard gravitational parameter for Earth [m^3/s^2].
       (**Default** = `m0`)
* `max_iter`: Maximum number of iterations allowed in the Newton-Raphson
              algorithm. (**Default** = 20)
* `tol`: Tolerance to stop the Newton-Raphson algorithm. (**Default** = 1e-10)

"""
function angvel_to_a(n::Number, e::Number, i::Number, pert::Symbol = :J2;
                     μ::Number = m0, max_iter::Int = 20, tol::Number = 1e-10)

    if pert == :J0

        a = (μ/n^2)^(1/3)
        return a

    elseif pert == :J2
        # Get the semi-major axis [m] that will provide the mean motion `n`
        # using perturbation terms up to J2.
        #
        # This can only be done using a numerical algorithm to solve the
        # following equation for `a`:
        #
        #           -                                                      -
        #          |      3        sqrt(1-e²).(3cos²(i)-1) + (5cos²(i)-1))  |
        #   n = n₀ | 1 + ---. J2 .----------------------------------------- |
        #          |      4                    a^2.(1-e^2)^2                |
        #           -                                                      -
        #
        #         sqrt(μ)
        #   n₀ = ----------
        #         a^(3/2)
        #

        # To improve algorithm stability, we will compute the normalized
        # semi-major axis, i.e. `a/R0`.

        # Auxiliary variables to solve for the semi-major axis.
        sqrt_μ = sqrt(μ/R0^3)
        cos_i  = cos(i)
        K      = 3/4*J2*( sqrt(1-e^2)*(3cos_i^2-1) + (5cos_i^2-1) )/(1-e^2)^2

        # Initial guess using a non-perturbed orbit.
        a = (μ/n^2)^(1/3)/R0

        # Newton-Raphson algorithm.
        #
        # Notice that we will allow, at most, 20 iterations.
        for k = 1:20
            # Auxiliary variables.
            ap3o2  = a^(3/2)
            ap5o2  = ap3o2*a # -> a^(5/2)
            ap7o2  = ap5o2*a # -> a^(7/2)
            ap11o2 = ap7o2*a # -> a^(11/2)

            # Compute the residue.
            res = n - sqrt_μ/ap3o2 - K*sqrt_μ/ap7o2

            # Compute the Jacobian of the function.
            df = +3/2*sqrt_μ/ap5o2 + 7/2*K*sqrt_μ/ap11o2

            # Compute the new estimate.
            a = a - res/df

            (abs(res) < tol) && break
        end

        return a*R0

    elseif pert == :J4

        # Auxiliary variables
        sqrt_μ = sqrt(μ/R0^3)
        sin_i  = sin(i)
        sin_i² = sin_i^2
        sin_i⁴ = sin_i^3
        e²     = e^2
        e⁴     = e^4
        aux    = 1-e²
        saux   = sqrt(aux)

        # Get the semi-major axis using J4 perturbation theory [er].
        #
        # This can only be done using a numerical algorithm to solve the
        # following equation for `a`:
        #
        #            .   .
        #   n = n₀ + ω + M₀ ,
        #
        #         sqrt(μm)
        #   n₀ = ---------- ,
        #          a^(3/2)
        #
        # and the time-derivatives of the argument of perigee and the initial
        # mean anomaly is compute considering the terms J2, J4, and J2².

        # To improve algorithm stability, we will compute the normalized
        # semi-major axis, i.e. `a/R0`.

        K₁ = 3/4*J2/aux^2*( (4 - 5sin_i²) + sqrt(1-e²)*(2 - 3sin_i²) )

        K₂ = 3/512*J2^2/aux^4*(224e² + (3040 - 144e²)*sin_i² - (3560 + 180e²)*sin_i⁴ +
                               1/saux*((         320e² - 280e⁴) +
                                       ( 1600 - 1568e² + 328e⁴)*sin_i² +
                                       (-2096 + 1072e² +  79e⁴)*sin_i⁴))

        K₃ = 15/128*J4/aux^4*(64 + 72e² - (248 + 252e²)*sin_i² + (196 + 189e²)*sin_i⁴ +
                              e²/saux*(-8 + 40sin_i - 35sin_i²))

        # Initial guess using a non-perturbed orbit.
        a = (μ/n^2)^(1/3)/R0

        # Newton-Raphson algorithm.
        #
        # Notice that we will allow, at most, 20 iterations.
        for k = 1:20
            # Auxiliary variables.
            ap3o2  = a^(3/2)
            ap5o2  =  ap3o2*a  # -> a^(5/2)
            ap7o2  =  ap5o2*a  # -> a^(7/2)
            ap9o2  =  ap7o2*a  # -> a^(9/2)
            ap11o2 =  ap9o2*a  # -> a^(11/2)
            ap13o2 = ap11o2*a  # -> a^(13/2)

            # Compute the residue.
            res = n - sqrt_μ*(1/ap3o2 + K₁/ap7o2 + (K₂ + K₃)/ap11o2)

            # Compute the Jacobian of the function.
            df = sqrt_μ/2*( 3ap5o2 + 7K₁/ap9o2 + 11*(K₂ + K₃)/ap13o2)

            # Compute the new estimate.
            a = a - res/df

            (abs(res) < tol) && break
        end

        return a*R0
    else
        throw(ArgumentError("The perturbation parameter $pert is not defined."))
    end
end

"""
    dargp(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    dargp(orb::Orbit, pert::Symbol = :J2)

Compute the time-derivative of the argument of perigee [rad/s] of an orbit with
semi-major axis `a` [m], eccentricity `e`, and inclination `i` [rad], using the
perturbation terms specified by the symbol `pert`. The orbit can also be
specified by `orb` (see [`Orbit`](@ref)).

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.
* `:J4`: Consider the perturbation terms J2, J4, and J2².

If `pert` is omitted, then it defaults to `:J2`.

"""
@inline function dargp(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    # Perturbation computed using a Keplerian orbit.
    if pert == :J0
        return 0.0
    # Perturbation computed using perturbations terms up to J2.
    elseif pert == :J2
        # Semi-lactum rectum.
        p = a*(1-e^2)

        # Unperturbed orbit period.
        n0 = angvel(a, e, i, :J0)

        # Perturbation of the argument of perigee.
        return 3R0^2*J2/(4p^2)*n0*(5cos(i)^2-1)

    elseif pert == :J4
        # Auxiliary variables
        e²     = e^2
        e⁴     = e^4
        sin_i  = sin(i)
        sin_i² = sin_i^2
        sin_i⁴ = sin_i^4
        aux    = (1-e²)
        saux   = sqrt(aux)
        p      = (a/R0)*aux
        p²     = p^2
        p⁴     = p^4

        # Unperturbed orbit period.
        n0 = angvel(a, e, i, :J0)

        # Perturbation of the argument of perigee.
        δω =  3/  4*n0*J2  /p²*(4 - 5sin_i²) +
              9/384*n0*J2^2/p⁴*(     56e² + (760 -  36e²)*sin_i² - (890 +  45*e²)*sin_i⁴) -
             15/128*n0*J4  /p⁴*(64 + 72e² - (248 + 252e²)*sin_i² + (196 + 189*e²)*sin_i⁴)

        return δω
    else
        throw(ArgumentError("The perturbation parameter $pert is not defined."))
    end
end

function dargp(orb::Orbit, pert::Symbol = :J2)
    # Convert first to Keplerian elements.
    k = convert(KeplerianElements, orb)
    return dargp(get_a(k), get_e(k), get_i(k), pert)
end

"""
    draan(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    draan(orb::Orbit, pert::Symbol = :J2)

Compute the time-derivative of the right ascension of the ascending node [rad/s]
of an orbit with semi-major axis `a` [m], eccentricity `e`, and inclination `i`
[rad], using the perturbation terms specified by the symbol `pert`. The orbit
can also be specified by `orb` (see [`Orbit`](@ref)).

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.
* `:J4`: Consider the perturbation terms J2, J4, and J2².

If `pert` is omitted, then it defaults to `:J2`.

"""
@inline function draan(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    # Perturbation computed using a Keplerian orbit.
    if pert == :J0
        return 0.0
    # Perturbation computed using perturbations terms up to J2.
    elseif pert == :J2
        # Semi-lactum rectum.
        p = a*(1-e^2)

        # Unperturbed orbit period.
        n0 = angvel(a, e, i, :J0)

        # Perturbation of the right ascension of the ascending node.
        return -3/2*R0^2/(p^2)*n0*J2*cos(i)

    # Perturbation computed using perturbation terms J2, J4, and J2².
    elseif pert == :J4
        # Auxiliary variables
        sin_i, cos_i = sincos(i)
        sin_i²       = sin_i^2
        e²           = e^2
        p            = (a/R0)*(1-e²)
        p²           = p^2
        p⁴           = p^4

        # Unperturbed orbit period.
        n0 = angvel(a, e, i, :J0)

        # Perturbation of the right ascension of the ascending node.
        δΩ = -3/2 *n0*J2  /p²*cos_i +
              3/32*n0*J2^2/p⁴*cos_i*(12 -  4e² - (80 +  5e²)*sin_i²) +
             15/32*n0*J4  /p⁴*cos_i*( 8 + 12e² - (14 + 21e²)*sin_i²)

        return δΩ
    else
        throw(ArgumentError("The perturbation parameter $pert is not defined."))
    end
end

function draan(orb::Orbit, pert::Symbol = :J2)
    # Convert first to Keplerian elements.
    k = convert(KeplerianElements, orb)
    return draan(get_a(k), get_e(k), get_i(k), pert)
end

"""
    period(a::Number, e::Number, i::Number, pert::Symbol = :J2)

Compute the period [s] of an object in an orbit with semi-major axis `a` [m],
eccentricity `e`, and inclination `i` [rad], using the perturbation terms
specified by the symbol `pert`. The orbit can also be specified by `orb` (see
[`Orbit`](@ref)).

pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.
* `:J4`: Consider the perturbation terms J2, J4, and J2².

If `pert` is omitted, then it defaults to `:J2`.

"""
@inline function period(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    n = angvel(a, e, i, pert)
    return 2π/n
end

function period(orb::Orbit, pert::Symbol = :J2)
    # Convert first to Keplerian elements.
    k = convert(KeplerianElements, orb)
    return period(get_a(k), get_e(k), get_i(k), pert)
end
