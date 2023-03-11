# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Auxiliary functions for Sun synchronous orbit computations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export sun_sync_orbit_from_ang_vel
export sun_sync_orbit_semi_major_axis
export sun_sync_orbit_inclination

"""
    sun_sync_orbit_from_ang_vel(n::T1, e::T2 = 0; kwargs...) where {T1 <: Number, T2 <: Number}

Compute the a Sun synchronous orbit given the angular velocity `n` [rad/s] and
the orbit eccentricity `e` [ ]. If the latter is omitted, the orbit is
considered circular, i.e., `e = 0`.

This function returns the orbit semi-major axis [km] and the orbit inclination
[rad] that provides a Sun synchronous orbit. It also returns a `Bool` indicating
whether the numerical algorithm has converged.

!!! note
    Internaly, this function uses the precision obtained by promoting `T1` and
    `T2` to a float-pointing number.

# Keywords

- `max_iterations::Number`: Maximum number of iterations in the Newton-Raphson
    method. (**Default** = 3)
- `tol::Union{Nothing, Number}`: Tolerance to verify if the numerical method has
    converged. If it is `nothing`, `eps(T)` will be used, where `T` is the
    internal type for the computations. (**Default** = 1e-18)

# Extended help

## Remarks

This algorithm uses the Newton-Raphson method to compute the orbit semi-major
axis and inclination that leads to a Sun synchronous orbit with the angular
velocity `n`. The theory here considers only the perturbation terms up to J2.
"""
function sun_sync_orbit_from_ang_vel(
    n::T1,
    e::T2 = 0;
    max_iterations::Number = 30,
    tol::Union{Nothing, Number} = nothing,
) where {T1 <: Number, T2 <: Number}

    T = float(promote_type(T1, T2))

    # Check if the arguments are valid.
    n <= 0        && throw(ArgumentError("The angular velocity must be greater than 0."))
    !(0 <= e < 1) && throw(ArgumentError("The eccentricity must be within the interval [0, 1)."))

    # Obtain the tolerance.
    if isnothing(tol)
        tol = eps(T)
    end

    # Auxiliary variables.
    sqrt_m0   = √T(m0)
    aux_e     = 1 - T(e)^2
    sqrt_1_e2 = √aux_e

    # Auxiliary constant to compute the functions.
    k₁ = 3T(R0)^2 * T(J2) * sqrt_m0 / (4aux_e^2)

    # Solve for zeros of f1 and f2 using Newton-Raphson method
    # ==========================================================================
    #
    # The function `f1` is the residue related to the time derivative of RAAN,
    # and the function `f2` is the residue related to the angular velocity.

    # Initial guess based on the unperturbed model.
    a_k = (T(m0) / T(n)^2)^T(1 / 3)
    i_k = acos(-T(ne) * a_k^T(7 / 2) / (2k₁))

    # By setting the initial values of `f1` and `f2` to `10tol`, we assure that
    # the loop will be executed at least one time.
    f₁ = 10T(tol)
    f₂ = 10T(tol)

    # Loop
    it = 0
    converged = true

    while (abs(f₁) > tol) || (abs(f₂) > tol)
        a_k_₁ = a_k
        i_k_₁ = i_k

        # Auxiliary variables.
        sin_i, cos_i = sincos(i_k_₁)
        sin_i²       = sin_i^2
        cos_i²       = cos_i^2
        a_pm1d5      = 1 / (a_k_₁ * √a_k_₁) # ..................... a_k_₁^(-1.5)
        a_pm2d5      = a_pm1d5 / a_k_₁      # ..................... a_k_₁^(-2.5)
        a_pm3d5      = a_pm2d5 / a_k_₁      # ..................... a_k_₁^(-3.5)
        a_pm4d5      = a_pm3d5 / a_k_₁      # ..................... a_k_₁^(-4.5)
        aux_1        = k₁ * sqrt_1_e2 * (3cos_i² - 1)
        aux_2        = k₁ * (5cos_i² - 1)

        # Compute the Jacobian.
        df1da = -7k₁ * a_pm4d5 * cos_i
        df1di = -2k₁ * a_pm3d5 * sin_i
        df2da = - T(3 / 2) * sqrt_m0 * a_pm2d5 -
                  T(7 / 2) * (aux_1 * a_pm4d5 + aux_2 * a_pm4d5)
        df2di = -(6k₁ * sqrt_1_e2 * a_pm3d5 + 10k₁ * a_pm4d5) * cos_i * sin_i

        J_k_₁ = @SMatrix T[df1da df1di
                           df2da df2di]

        # Compute the functions.
        f₁ = T(ne) + 2k₁ * a_pm3d5 * cos_i
        f₂ = sqrt_m0 * a_pm1d5 + aux_1 * a_pm3d5 + aux_2 * a_pm3d5 - n

        # Compute the new estimate.
        vf    = @SVector T[f₁, f₂]
        v_k_₁ = @SVector T[a_k_₁, i_k_₁]
        v_k   = v_k_₁ - inv(J_k_₁) * vf
        a_k   = v_k[1]
        i_k   = v_k[2]

        # Check if the maximum number of iterations has been reached.
        it += 1

        # If the maximum number of iterations allowed has been reached, then
        # indicate that the solution did not converged and exit loop.
        if (it >= max_iterations)
            converged = false
            break
        end
    end

    !converged && @warn("""
        The algorithm to compute the Sun synchornous orbit has not converged!
        Residues: f₁ = $f₁, f₂ = $f₂"""
    )

    # Return.
    return a_k, i_k, converged
end

"""
    sun_sync_orbit_semi_major_axis(i::T1, e::T2 = 0) where {T1 <: Number, T2 <: Number}

Compute the Sun synchronous orbit semi-major axis [m] given the orbit inclination
`i` [rad] and the orbit eccentricity `e` [ ]. If the latter is omitted, the
orbit is considered circular, i.e., `e = 0`.

!!! note
    Internaly, this function uses the precision obtained by promoting `T1` and
    `T2` to a float-pointing number.

# Remarks

The theory here considers only the perturbation terms up to J2.
"""
function sun_sync_orbit_semi_major_axis(i::T1, e::T2 = 0) where {T1 <: Number, T2 <: Number}
    T = float(promote_type(T1, T2))

    # Check the inputs.
    !(0 <= e < 1) && throw(ArgumentError("The eccentricity must be within the interval [0, 1)."))

    # Auxiliary variables.
    sqrt_m0 = √T(m0)

    # Auxiliary constant.
    k₁ = 3T(R0)^2 * T(J2) * sqrt_m0 / (4 * (1 - T(e)^2)^2)

    # Compute the semi-major axis.
    a = (-2cos(i) * k₁ / T(ne))^T(2 / 7)

    # Check if the orbit is valid.
    if (a * (1 - e) <  R0)
        throw(ErrorException("It was not possible to find a valid Sun-synchronous orbit!"))
    end

    return a
end

"""
    sun_sync_orbit_inclination(a::T1, e::T2 = 0) where {T1 <: Number, T2 <: Number}

Compute the Sun synchronous orbit inclination [rad] given the orbit semi-major
axis `a` [m] and the orbit eccentricity `e` [ ]. If the latter is omitted, the
orbit is considered circular, i.e., `e = 0`.

!!! note
    Internaly, this function uses the precision obtained by promoting `T1` and
    `T2` to a float-pointing number.

# Remarks

The theory here considers only the perturbation terms up to J2.
"""
function sun_sync_orbit_inclination(a::T1, e::T2 = 0) where {T1 <: Number, T2 <: Number}
    T = float(promote_type(T1, T2))

    # Check if the arguments are valid.
    (a * (1 - e) <= R0) && throw(ArgumentError("The perigee must be larger than the Earth's radius."))
    !(0 <= e < 1) && throw(ArgumentError("The eccentricity must be within the interval [0, 1)."))

    # Auxiliary variables.
    sqrt_m0 = √T(m0)

    # Auxiliary constant.
    k₁ = 3T(R0)^2 * T(J2) * sqrt_m0 / (4 * (1 - T(e)^2)^2)

    # Compute the inclination.
    cos_i = -T(ne) * T(a)^3 * √T(a) / (2k₁)

    # Check if the cosine is valid.
    if ((cos_i < -1) || (cos_i > 1))
        throw(ErrorException("It was not possible to find a Sun-synchronous orbit!"))
    end

    return acos(cos_i)
end
