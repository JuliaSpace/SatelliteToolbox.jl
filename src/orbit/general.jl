#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions to compute general values related to the orbit.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export @check_orbit
export angvel, angvel_to_a, dArgPer, dRAAN, period

################################################################################
#                                  Overloads
################################################################################

copy(orb::Orbit) = Orbit(orb.t, orb.a, orb.e, orb.i, orb.Ω, orb.ω, orb.f)

function display(orb::Orbit)
    d2r = 180/π

    # Definition of colors that will be used for printing.
    b = "\x1b[1m"
    d = "\x1b[0m"
    g = "\x1b[1m\x1b[32m"
    y = "\x1b[1m\x1b[33m"

    # Print the TLE information.
    println()
    print("                 $(g)Orbit$(d)\n")
    print("$(y)  ======================================$(d)\n")
    print("$(b)                  t = $(d)"); @printf("%14.5f\n",    orb.t)
    print("$(b)    Semi-major axis = $(d)"); @printf("%13.4f km\n", orb.a/1000)
    print("$(b)       Eccentricity = $(d)"); @printf("%15.6f\n",    orb.e)
    print("$(b)        Inclination = $(d)"); @printf("%13.4f ˚\n",  orb.i*d2r)
    print("$(b)               RAAN = $(d)"); @printf("%13.4f ˚\n",  orb.Ω*d2r)
    print("$(b)    Arg. of Perigee = $(d)"); @printf("%13.4f ˚\n",  orb.ω*d2r)
    print("$(b)       True Anomaly = $(d)"); @printf("%13.4f ˚\n",  orb.f*d2r)
end

################################################################################
#                                    Macros
################################################################################

"""
    macro check_orbit(a, e)

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
    function Orbit(a::Number, e::Number, i::Number, Ω::Number, ω::Number, f::Number)

Create an orbit with semi-major axis `a` [m], eccentricity `e`, inclination `i`
[rad], right ascension of the ascending node `Ω` [rad], argument of perigee `ω`
[rad], and true anomaly `f` [rad].

# Returns

An object of type `Orbit` with the specified orbit. The orbit epoch is defined
as 0.0.

"""
Orbit(a::Number, e::Number, i::Number, Ω::Number, ω::Number, f::Number) =
    Orbit(0.0, a, e, i, Ω, ω, f)

################################################################################
#                                  Functions
################################################################################

"""
    function angvel(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    function angvel(orb::Orbit, pert::Symbol = :J2)

Compute the angular velocity [rad/s] of an object in an orbit with semi-major
axis `a` [m], eccentricity `e`, and inclination `i` [rad], using the
perturbation terms specified by the symbol `pert`. The orbit can also be
specified by `orb`, which is an instance of the structure `Orbit`.

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.

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
    else
        throw(ArgumentError("The perturbation parameter $pert is not defined."))
    end

end

@inline angvel(orb::Orbit, pert::Symbol = :J2) = angvel(orb.a, orb.e, orb.i, pert)

"""
    function angvel_to_a(n::Number, e::Number, i::Number, pert::Symbol = :J2; μ::Number = m0)

Compute the semi-major axis that will provide an angular velocity `n` [rad/s] in
an orbit with eccentricity `e` and inclination `i` [rad], using the perturbation
terms specified by the symbol `pert`.

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.

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
    else
        throw(ArgumentError("The perturbation parameter $pert is not defined."))
    end
end

"""
    function dArgPer(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    function dArgPer(orb::Orbit, pert::Symbol = :J2)

Compute the time-derivative of the argument of perigee [rad/s] of an orbit with
semi-major axis `a` [m], eccentricity `e`, and inclination `i` [rad], using the
perturbation terms specified by the symbol `pert`. The orbit can also be
specified by `orb`, which is an instance of the structure `Orbit`.

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.

If `pert` is omitted, then it defaults to `:J2`.

"""
@inline function dArgPer(a::Number, e::Number, i::Number, pert::Symbol = :J2)
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
    else
        throw(ArgumentError("The perturbation parameter $pert is not defined."))
    end

end

@inline dArgPer(orb::Orbit, pert::Symbol = :J2) = dArgPer(orb.a, orb.e, orb.i, pert)

"""
    function dRAAN(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    function dRAAN(orb::Orbit, pert::Symbol = :J2)

Compute the time-derivative of the right ascension of the ascending node [rad/s]
of an orbit with semi-major axis `a` [m], eccentricity `e`, and inclination `i`
[rad], using the perturbation terms specified by the symbol `pert`. The orbit
can also be specified by `orb`, which is an instance of the structure `Orbit`.

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.

If `pert` is omitted, then it defaults to `:J2`.

"""
@inline function dRAAN(a::Number, e::Number, i::Number, pert::Symbol = :J2)
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
    else
        throw(ArgumentError("The perturbation parameter $pert is not defined."))
    end
end

@inline dRAAN(orb::Orbit, pert::Symbol = :J2) = dRAAN(orb.a, orb.e, orb.i, pert)


"""
    function period(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    function period(orb::Orbit, pert::Symbol = :J2)

Compute the period [s] of an object in an orbit with semi-major axis `a` [m],
eccentricity `e`, and inclination `i` [rad], using the perturbation terms
specified by the symbol `pert`. The orbit can also be specified by `orb`, which
is an instance of the structure `Orbit`.

`pert` can be:

* `:J0`: Consider a Keplerian orbit.
* `:J2`: Consider the perturbation terms up to J2.

If `pert` is omitted, then it defaults to `:J2`.

"""
@inline function period(a::Number, e::Number, i::Number, pert::Symbol = :J2)
    n = angvel(a, e, i, pert)
    2π/n
end

@inline period(orb::Orbit, pert::Symbol = :J2) = period(orb.a, orb.e, orb.i, pert)
