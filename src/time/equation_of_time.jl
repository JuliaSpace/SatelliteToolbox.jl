## Description #############################################################################
#
#   Compute the equation of time.
#
## References ##############################################################################
#
# [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
#     Microcosm Press, Hawthorne, CA.
#
############################################################################################

export equation_of_time

"""
    equation_of_time(t::Union{Number, DateTime}) -> T

Compute the difference between the Sun apparent local time and the Sun mean local time
[rad], which is called Equation of Time, at the time `t` [UT1], which can be represented by
a Julian Day or `DateTime`. The algorithm was adapted from **[1, p. 178, 277-279]**.

!!! note

    The output type `T` is the float type of the Julian Day. If `t` is a `DateTime`, `T` is
    `Float64`.

# References

- **[1]** Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications. 4th ed.
    Microcosm Press, Hawthorne, CA.
"""
equation_of_time(t::DateTime) = equation_of_time(datetime2julian(t))

function equation_of_time(jd::Number)
    T = float(typeof(jd))

    # Number of Julian centuries from J2000 epoch. The subtraction is performed before the
    # conversion to `T` to avoid losing precision if `jd` has a lower precision than
    # `JD_J2000`.
    t_ut1 = T((jd - JD_J2000) / 36525)

    # Mean longitude of the Sun [deg].
    λ_m = mod(T(280.460) + T(36000.771) * t_ut1, 360)

    # Mean anomaly of the Sun [rad].
    #
    # Here, we should use T_TDB (Barycentric Dynamical Time). However, it is sufficient to
    # use t_ut1 because this is a low precision computation [1].
    Ms = mod(T(357.5291092) + T(35999.05034) * t_ut1, 360) |> deg2rad

    # Auxiliary variables.
    sin_Ms  = sin(Ms)
    sin_2Ms = sin(2Ms)

    # Ecliptic longitude of the Sun [rad].
    λ_ecliptic = mod(λ_m + T(1.914666471) * sin_Ms + T(0.019994643) * sin_2Ms, 360)
    λ_ecliptic = deg2rad(λ_ecliptic)

    # Compute the equation of time [deg].
    eot = -T(1.914666471) * sin_Ms -
           T(0.019994643) * sin_2Ms +
           T(2.466) * sin(2λ_ecliptic) -
           T(0.0053) * sin(4λ_ecliptic)

    return deg2rad(eot)
end
