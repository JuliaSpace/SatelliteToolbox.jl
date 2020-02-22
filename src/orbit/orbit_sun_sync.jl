#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Many auxiliary functions for sun synchronous orbit computations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export compute_ss_orbit_by_ang_vel
export compute_ss_orbit_by_inclination
export compute_ss_orbit_by_semi_major_axis

"""
    compute_ss_orbit_by_ang_vel(n::Number, e::Number)

Compute the Sun-synchronous orbit given the angular velocity `n` [rad/s] and the
eccentricity `e`.

# Returns

* The semi-major axis [m].
* The inclination [rad].
* The residues of the two functions.
* A boolean variable that indicates if the numerical algorithm converged.

"""
function compute_ss_orbit_by_ang_vel(n::Number, e::Number)
    # Check if the arguments are valid.
    (n <= 0)          && throw(ArgumentError("The angular velocity must be greater than 0."))
    !( 0. <= e < 1. ) && throw(ArgumentError("The eccentricity must be within the interval 0 <= e < 1."))

    # Tolerance for the Newton-Raphson method.
    tol = 1e-18

    # Maximum number of iterations for the Newton-Raphson method.
    maxIt = 30

    # Auxiliary variables.
    sqrt_m0 = sqrt(m0)
    aux_e   = 1-e^2
    sqrt_e2 = sqrt(aux_e)

    # Auxiliary constant to compute the functions.
    K1 = 3R0^2*J2*sqrt_m0/(4aux_e^2)

    # Solve for zeros of f1 and f2 using Newton-Raphson method
    # ==========================================================================
    #
    # The function `f1` is the residue related to the time derivative of RAAN,
    # and the function `f2` is the residue related to the angular velocity.

    # Initial guess based on the unperturbed model.
    a_k = (m0/n^2)^(1/3)
    i_k = acos( -ne*a_k^(3.5)/(2K1) )

    # By setting the initial values of `f1` and `f2` to `10tol`, we assure that
    # the loop will be executed at least one time.
    f1 = 10tol
    f2 = 10tol

    # Loop
    it = 0
    converged = true

    while (abs(f1) > tol) || (abs(f2) > tol)
        a_k_1 = a_k
        i_k_1 = i_k

        # Auxiliary variables.
        sin_i, cos_i = sincos(i_k_1)
        sin_i²       = sin_i^2
        cos_i²       = cos_i^2
        a_pm1d5      = a_k_1^(-1.5)
        a_pm2d5      = a_pm1d5*a_k_1^(-1) # --> a_k_1^(-2.5)
        a_pm3d5      = a_pm2d5*a_k_1^(-1) # --> a_k_1^(-3.5)
        a_pm4d5      = a_pm3d5*a_k_1^(-1) # --> a_k_1^(-4.5)
        aux_1        = K1*sqrt_e2*(3cos_i²-1)
        aux_2        = K1*(5cos_i²-1)

        # Compute the Jacobian.
        df1da = -7K1*a_pm4d5*cos_i
        df1di = -2K1*a_pm3d5*sin_i
        df2da = -1.5*sqrt_m0*a_pm2d5 - 3.5aux_1*a_pm4d5 - 3.5aux_2*a_pm4d5
        df2di = -6K1*sqrt_e2*a_pm3d5*cos_i*sin_i - 10K1*a_pm4d5*cos_i*sin_i

        J_k_1 = SMatrix{2,2}(df1da, df1di,
                             df2da, df2di)'

        # Compute the functions.
        f1 = ne + 2K1*a_pm3d5*cos_i
        f2 = sqrt_m0*a_pm1d5 + aux_1*a_pm3d5 + aux_2*a_pm3d5 - n

        # Compute the new estimate.
        vf    = SVector{2}(f1,f2)
        v_k_1 = SVector{2}(a_k_1, i_k_1)
        v_k   = v_k_1 - inv(J_k_1)*vf
        a_k   = v_k[1]
        i_k   = v_k[2]

        # Check if the maximum number of iterations has been reached.
        it += 1

        # If the maximum number of iterations allowed has been reached, then
        # indicate that the solution did not converged and exit loop.
        if (it >= maxIt)
            converged = false
            break
        end
    end

    # Return.
    a_k, i_k, f1, f2, converged
end

"""
    compute_ss_orbit_by_inclination(i::Number, e::Number)

Compute the Sun-synchronous orbit given the inclination `i` [rad] and the
eccentricity `e`.

# Returns

The semi-major axis of the Sun-synchronous orbit [m].

"""
function compute_ss_orbit_by_inclination(i::Number, e::Number)
    !( 0. <= e < 1. ) && throw(ArgumentError("The eccentricity must be within the interval 0 <= e < 1."))

    # Auxiliary variables.
    sqrt_m0 = sqrt(m0)

    # Auxiliary constant.
    K1 = 3.0*R0^2*J2*sqrt_m0/(4.0*(1-e^2)^2)

    # Compute the semi-major axis.
    a = (-2*cos(i)*K1/ne)^(2.0/7.0)

    # Check if the orbit is valid.
    ( a*(1-e) <  R0 ) && throw(ErrorException("It was not possible to find a valid Sun-synchronous orbit with the inclination given."))

    # Return.
    a
end

"""
    compute_ss_orbit_by_semi_major_axis(a::Number, e::Number)

Compute the Sun-synchronous orbit given the semi-major axis `a` [m] and the
eccentricity `e`.

# Returns

The inclination of the Sun-synchronous orbit [rad].

"""
function compute_ss_orbit_by_semi_major_axis(a::Number, e::Number)
    # Check if the arguments are valid.
    (a*(1-e) <= R0)   && throw(ArgumentError("The perigee must be larger than the Earth radius."))
    !( 0. <= e < 1. ) && throw(ArgumentError("The eccentricity must be within the interval 0 <= e < 1."))

    # Auxiliary variables.
    sqrt_m0 = sqrt(m0)

    # Auxiliary constant.
    K1 = 3R0^2*J2*sqrt_m0/(4(1-e^2)^2)

    # Compute the inclination.
    cos_i = -ne*a^3*sqrt(a)/(2K1)

    # Check if -1 <= cos_i <= 1.
    ( (cos_i < -1) || (cos_i > 1) ) && throw(ErrorException("It was not possible to find a Sun-synchronous orbit with the semi-major axis given."))

    # Return.
    acos(cos_i)
end
