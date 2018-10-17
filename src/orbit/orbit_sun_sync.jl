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
    function compute_ss_orbit_by_ang_vel(n::Number, e::Number)

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
    if (n <= 0)
        throw(ArgumentError("The angular velocity must be greater than 0."))
    end

    if !( 0. <= e < 1. )
        throw(ArgumentError("The eccentricity must be within the interval 0 <= e < 1."))
    end

    # Tolerance for the Newton-Raphson method.
    tol = 1e-18

    # Maximum number of iterations for the Newton-Raphson method.
    maxIt = 30

    # Auxiliary variables.
    sqrt_m0::Float64 = sqrt(m0)
    sqrt_e2::Float64 = sqrt(1-e^2)

    # Auxiliary constant to compute the functions.
    K1::Float64 = 3.0*R0^2*J2*sqrt_m0/(4.0*(1-e^2)^2)

    # Declare the functions that must solved for 0.
    f1(a, i) = ne + 2.0*K1*a^(-3.5)*cos(i)
    f2(a, i) = sqrt_m0*a^(-1.5) + K1*sqrt_e2*(3.0*cos(i)^2-1.0)*a^(-3.5) + K1*(5.0*cos(i)^2-1.0)*a^(-3.5) - n

    # Declare the derivative of the functions.
    df1da(a, i) = -7.0*K1*a^(-4.5)*cos(i)
    df1di(a, i) = -2.0*K1*a^(-3.5)*sin(i)
    df2da(a, i) = -1.5*sqrt_m0*a^(-2.5)-3.5*K1*sqrt_e2*(3.0*cos(i)^2-1.0)*a^(-4.5)-3.5*K1*(5.0*cos(i)^2-1.0)*a^(-4.5)
    df2di(a, i) = -6.0*K1*sqrt_e2*a^(-3.5)*cos(i)*sin(i)-10.0*K1*a^(-3.5)*cos(i)*sin(i)

    # Jacobian.
    Jf1f2(a, i) = [df1da(a,i) df1di(a,i);
                   df2da(a,i) df2di(a,i);]

    ############################################################################
    # Solve for zeros of f1 and f2 using Newton-Raphson method.
    ############################################################################

    # Maximum number of

    # Initial guess based on the unperturbed model.
    a_k::Float64 = (m0/n^2)^(1/3)
    i_k::Float64 = acos( -ne*a_k^(3.5)/(2*K1) )

    # Loop
    it = 0;
    converged = true

    while (abs(f1(a_k, i_k)) > tol) || (abs(f2(a_k, i_k)) > tol)
        a_k_1 = a_k
        i_k_1 = i_k

        # Compute the Jacobian.
        J_k_1::Matrix{Float64} = Jf1f2(a_k_1, i_k_1)

        # Compute the new estimate.
        (a_k, i_k) = [a_k_1; i_k_1] - inv(J_k_1)*[f1(a_k_1, i_k_1); f2(a_k_1, i_k_1)];

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
    a_k, i_k, f1(a_k, i_k), f2(a_k, i_k), converged
end

"""
    function compute_ss_orbit_by_inclination(i::Number, e::Number)

Compute the Sun-synchronous orbit given the inclination `i` [rad] and the
eccentricity `e`.

# Returns

The semi-major axis of the Sun-synchronous orbit [m].

"""
function compute_ss_orbit_by_inclination(i::Number, e::Number)
    if !( 0. <= e < 1. )
        throw(ArgumentError("The eccentricity must be within the interval 0 <= e < 1."))
    end

    # Auxiliary variables.
    sqrt_m0 = sqrt(m0)

    # Auxiliary constant.
    K1 = 3.0*R0^2*J2*sqrt_m0/(4.0*(1-e^2)^2)

    # Compute the semi-major axis.
    a = (-2*cos(i)*K1/ne)^(2.0/7.0)

    # Check if the orbit is valid.
    if ( a*(1.0-e) <  R0 )
        throw(ErrorException("It was not possible to find a valid Sun-synchronous orbit with the inclination given."))
    end

    # Return.
    a
end

"""
    function compute_ss_orbit_by_semi_major_axis(a::Number, e::Number)

Compute the Sun-synchronous orbit given the semi-major axis `a` and the
eccentricity `e`.

# Args

* `a`: Semi-major axis [m].
* `e`: Eccentricity.

# Returns

The inclination of the Sun-synchronous orbit [rad].

"""
function compute_ss_orbit_by_semi_major_axis(a::Number, e::Number)
    # Check if the arguments are valid.
    if (a*(1.0-e) <= R0)
        throw(ArgumentError("The perigee must be larger than the Earth radius."))
    end

    if !( 0. <= e < 1. )
        throw(ArgumentError("The eccentricity must be within the interval 0 <= e < 1."))
    end

    # Auxiliary variables.
    sqrt_m0 = sqrt(m0)

    # Auxiliary constant.
    K1 = 3.0*R0^2*J2*sqrt_m0/(4.0*(1-e^2)^2)

    # Compute the inclination.
    cos_i = -ne*a^3*sqrt(a)/(2*K1)

    # Check if -1 <= cos_i <= 1.
    if ( (cos_i < -1) || (cos_i > 1) )
        throw(ErrorException("It was not possible to find a Sun-synchronous orbit with the semi-major axis given."))
    end

    # Return.
    acos(cos_i)
end
