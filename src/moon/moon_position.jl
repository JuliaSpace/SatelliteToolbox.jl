# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Compute the moon position.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorne, CA.
#   [2] Meeus, J (1998). Astronomical algorithms. Willmann-Bell, Inc, Richmond,
#       VA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export moon_position_i

"""
    moon_position_i(JD_TDB::Number[, model])

Compute the Moon position represented in the IAU-76/FK5 MOD (mean-equator,
mean-equinox of date) at the Julian Day `JD_TDB` (Barycentric Dynamical Time).

The `model` must be `Val(:Meeus)` or `Val(:fast)`. `Val(:Meeus)` uses the
algorithm in [2, p. 337] that provides an accuracy of 10" in the longitude and
4" in the latitude (the reference does not mention the timespan). `Val(:fast)`
uses the algorithm in [1, p. 288] that is 10x faster than `Val(:Meeus)` but can
lead to errors of 0.3° in longitude and 0.2° in latitude.

"""
moon_position_i(JD_TDB::Number) = moon_position_i(JD_TDB, Val(:Meeus))

function moon_position_i(JD_TDB::Number, ::Val{:Meeus})
    # Number of Julian centuries from J2000 epoch.
    T_TDB = (JD_TDB - JD_J2000)/36525.0

    # Moon's mean latitude reffered to the mean equinox of data [deg].
    L´ = @evalpoly(T_TDB, 218.3164477, 481_267.88123421, -0.0015786, 1/538841, -1/65_194_000)

    # Mean elongation of the Moon [deg].
    D = @evalpoly(T_TDB, 297.8501921, 445_267.1114034, -0.0018819, 1/545868, -1/113_065_000)

    # Sun's mean anomaly [deg].
    M = @evalpoly(T_TDB, 357.5291092, 35_999.0502909, -0.0001536, 1/24_490_000)

    # Moon's mean anomaly [deg].
    M´ = @evalpoly(T_TDB, 134.9633964, 477_198.8675055, 0.0087414, 1/69_699, -1/14_712_000)

    # Moon's argument of latitude (mean distance of the Moon from its ascending
    # node) [deg].
    F = @evalpoly(T_TDB, 93.2720950, 483_202.0175233, -0.0036539, -1/3_526_000, 1/863_310_000)

    # Obliquity of the ecliptic [deg].
    ϵ = @evalpoly(T_TDB, 23.439_291, -0.013_004_2, -1.64e-7, +5.04e-7)

    # Additional arguments required for the algorithm [deg].
    A₁ = @evalpoly(T_TDB, 119.75, 131.849)
    A₂ = @evalpoly(T_TDB,  53.09, 479_264.290)
    A₃ = @evalpoly(T_TDB, 313.45, 481_266.484)

    # Convert everything to [rad] and limit the angles between [0, 2π].
    L´ = mod2pi(deg2rad(L´))
    D  = mod2pi(deg2rad(D))
    M  = mod2pi(deg2rad(M))
    M´ = mod2pi(deg2rad(M´))
    F  = mod2pi(deg2rad(F))
    ϵ  = mod2pi(deg2rad(ϵ))
    A₁ = mod2pi(deg2rad(A₁))
    A₂ = mod2pi(deg2rad(A₂))
    A₃ = mod2pi(deg2rad(A₃))

    # Term used to corret the arguments of the angle M that depends on the
    # Earth's orbit eccentricity around the Sun.
    E  = @evalpoly(T_TDB, 1, -0.002516, 0.0000074)
    E² = E * E

    # Compute the sum of the terms in the tables 47.A and 47.B [2].

    # Auxiliary variables to simplify the expressions.
    tab = _tab_47a

    # Longitude and distance of the Moon.
    num_terms = size(tab)[2]

    Σl = 0.0
    Σr = 0.0

    @inbounds for k = 1:num_terms
        aD  = tab[1, k]
        aM  = tab[2, k]
        aM´ = tab[3, k]
        aF  = tab[4, k]

        arg = mod2pi(aD * D + aM * M + aM´ * M´ + aF * F)

        # Check if we need to apply the correction `E`.
        E_corr = if ((aM == 1) || (aM == -1))
            E
        elseif ((aM == 2) || (aM == -2))
            E²
        else
            one(E)
        end

        # Compute the quantities.
        sin_arg, cos_arg = sincos(arg)

        Σl += tab[5, k] * E_corr * sin_arg
        Σr += tab[6, k] * E_corr * cos_arg
    end

    # Auxiliary variables to simplify the expressions.
    tab = _tab_47b

    # Latitude of the Moon.
    num_terms = size(tab)[2]

    Σb = 0.0

    @inbounds for k = 1:num_terms
        aD  = tab[1, k]
        aM  = tab[2, k]
        aM´ = tab[3, k]
        aF  = tab[4, k]

        arg = mod2pi(aD * D + aM * M + aM´ * M´ + aF * F)

        # Check if we need to apply the correction `E`.
        E_corr = if ((aM == 1) || (aM == -1))
            E
        elseif ((aM == 2) || (aM == -2))
            E²
        else
            one(E)
        end

        Σb += tab[5, k] * E_corr * sin(arg)
    end

    # Apply the corrections to the terms.
    Σl +=  3958sin(A₁) + 1962sin(L´ - F) + 318sin(A₂)
    Σb += -2235sin(L´) + 382sin(A₃) + 175sin(A₁ - F) + 175sin(A₁ + F) +
           127sin(L´ - M´) - 115sin(L´ + M´)

    # Convert to [rad].
    Σl = mod2pi(deg2rad(Σl / 1_000_000))
    Σb = mod2pi(deg2rad(Σb / 1_000_000))

    # Compute the Moon coordinates [rad] and [m].
    λ = mod2pi(L´ + Σl)
    β = Σb
    Δ = 385_000.56e3 + Σr

    # Compute the Moon vector in MOD [m].
    #
    # Notice that λ and β provide us the geocentric latitude and longitude of
    # the Moon w.r.t. the mean equator of date in the ecliptic plane. Hence, we
    # need also to rotate the mean ecliptic to obtain the vector in the MOD.
    sin_λ, cos_λ = sincos(λ)
    sin_β, cos_β = sincos(β)
    sin_ϵ, cos_ϵ = sincos(ϵ)

    r_moon_i = SVector(
        Δ * cos_β * cos_λ,
        Δ * (cos_ϵ * cos_β * sin_λ - sin_ϵ * sin_β),
        Δ * (sin_ϵ * cos_β * sin_λ + cos_ϵ * sin_β)
    )

    return r_moon_i
end

function moon_position_i(JD_TDB::Number, ::Val{:fast})
    # Number of Julian centuries from J2000 epoch.
    T_TDB = (JD_TDB - JD_J2000)/36525.0

    # Auxiliary computation to improve performance.
    sin1, cos1 = sincos(deg2rad(134.9 + 477_198.85T_TDB))
    sin2, cos2 = sincos(deg2rad(259.2 - 413_335.38T_TDB))
    sin3, cos3 = sincos(deg2rad(235.7 + 890_534.23T_TDB))
    sin4, cos4 = sincos(deg2rad(269.9 + 954_397.70T_TDB))
    sin5       = sin(deg2rad(357.5 +  35_999.05T_TDB))
    sin6       = sin(deg2rad(186.6 + 966_404.05T_TDB))

    # Ecliptic latitude of the Moon [deg].
    λₑ = 218.32 + 481_267.8813T_TDB +
         6.29sin1 - 1.27sin2 + 0.66sin3 + 0.21sin4 - 0.19sin5 - 0.11sin6

    # Ecliptic longitude of the Moon [deg].
    ϕₑ = 5.13sin(deg2rad( 93.3 + 483_202.03T_TDB)) +
         0.28sin(deg2rad(228.2 + 960_400.87T_TDB)) -
         0.28sin(deg2rad(318.3 +   6_003.18T_TDB)) -
         0.17sin(deg2rad(217.6 - 407_332.20T_TDB))

    # Parallax [deg].
    P = 0.9508 + 0.0518cos1 + 0.0095cos2 + 0.0078cos3 + 0.0028cos4

    # Obliquity of the ecliptic [deg].
    ϵ = @evalpoly(T_TDB, 23.439_291, -0.013_004_2, -1.64e-7, +5.04e-7)

    # Convert to radians and limit to the interval [0,2π].
    λₑ = mod2pi(deg2rad(λₑ))
    ϕₑ = mod2pi(deg2rad(ϕₑ))
    P  = mod2pi(deg2rad( P))
    ϵ  = mod2pi(deg2rad( ϵ))

    # Compute the distance from Earth to the Moon [m].
    r = R0 / sin(P)

    # Auxiliary variables.
    sin_λ, cos_λ = sincos(λₑ)
    sin_ϕ, cos_ϕ = sincos(ϕₑ)
    sin_ϵ, cos_ϵ = sincos(ϵ)

    # Compute the Moon vector represented in MOD (IAU-76/KF5 mean-equator,
    # mean-equinox of date).
    r_mod = SVector(
        r * (cos_ϕ * cos_λ),
        r * (cos_ϵ * cos_ϕ * sin_λ - sin_ϵ * sin_ϕ),
        r * (sin_ϵ * cos_ϕ * sin_λ + cos_ϵ * sin_ϕ)
    )

    return r_mod
end
