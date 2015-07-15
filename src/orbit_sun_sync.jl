#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#    Many auxiliary functions for sun-synchronous orbit computations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2015-07-15: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br> The function
#    `minimum_swath_grss` no longer return the minimum Half-FOV. This
#    functionality was moved to the function `minimum_half_fov`.
#
# 2014-12-18: Ronan Arraes Jardim Chagas <ronan.chagas@inpe.br>
#    Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

#==#
#
# @brief Compute the distance between adjacent ground tracks in a given latitude
# for a ground repeating, sun synchronous orbit.
#
# @param[in] T Orbit period [s].
# @param[in] i Inclination [rad].
# @param[in] To Orbit cycle [days].
# @param[in] lat Latitude [rad].
#
# @return The distance between adjacent ground tracks [m].
#
# @remarks The function does not check if the orbit is a GRSS orbit.
#
#==#

function adjacent_track_distance_grss(T::Real, i::Real, To::Real, lat::Real)
    # Angle between two adjacent traces.
    theta = (T/To)*pi/43200.0

    # Angle of swath.
    beta = asin(sin(theta)*sin(i))

    # Minimum swath.
    beta*R0*cos(lat)
end

#==#
#
# @brief Compute the distance between adjacent ground tracks in a given latitude
# for a ground repeating, sun synchronous orbit.
#
# @param[in] a Semi-major axis [m].
# @param[in] e Eccentricity [s].
# @param[in] i Inclination [rad].
# @param[in] To Orbit cycle [days].
# @param[in] lat Latitude [rad].
#
# @return The distance between adjacent ground tracks [m].
#
# @remarks The function does not check if the orbit is a GRSS orbit.
#
#==#

function adjacent_track_distance_grss(a::Real, e::Real, i::Real, To::Real,
                                      lat::Real)
    # Orbit period.
    T = t_J2(a,e,i)

    # Compute the distance between adjacent ground tracks.
    adjacent_track_distance_grss(T, i, To, lat)
end

#==#
#
# @brief Compute the angle between two adjacent ground tracks in a given
# latitude measured on the satellite position for a ground repeating, sun
# synchronous orbit.
#
# @param[in] h Orbit altitude in the Equator [m].
# @param[in] T Orbit period [s].
# @param[in] i Inclination [rad].
# @param[in] To Orbit cycle [days].
# @param[in] lat Latitude [rad].
#
# @return The angle between adjacent ground tracks [rad].
#
# @remarks The function does not check if the orbit is a GRSS orbit.
#
#==#

function adjacent_track_angle_grss(h::Real, T::Real, i::Real, To::Real,
                                   lat::Real)
    # Angle between two adjacent traces.
    theta = (T/To)*pi/43200.0

    # Angle of swath.
    beta = asin(sin(theta)*sin(i))

    # Get the Earth radius in the plane perpendicular to the Equatorial plane.
    Rp = R0*cos(lat)
    
    # Compute the angle between the two ground tracks.
    b = sqrt(Rp^2+(Rp+h)^2-2*Rp*(Rp+h)*cos(beta/2.0))
    asin(Rp/b*sin(beta/2.0))
end

#==#
#
# @brief Compute the angle between two adjacent ground tracks in a given
# latitude measured on the satellite position for a ground repeating, sun
# synchronous orbit.
#
# @param[in] h Orbit altitude in the Equator [m].
# @param[in] a Semi-major axis [m].
# @param[in] e Eccentricity [s].
# @param[in] i Inclination [rad].
# @param[in] To Orbit cycle [days].
# @param[in] lat Latitude [rad].
#
# @return The angle between adjacent ground tracks [rad].
#
# @remarks The function does not check if the orbit is a GRSS orbit.
#
#==#

function adjacent_track_angle_grss(h::Real, a::Real, e::Real, i::Real,
                                   To::Integer, lat::Real)
    # Period (J2).
    T = t_J2(a,e,i)

    # Compute the minimum half FOV.
    adjacent_track_angle_grss(h, T, i, To, lat)
end

#==#
# 
# @brief Compute the sun-synchronous orbit given the angular velocity and the
# eccentricity.
#
# @param[in] n Angular velocity [rad/s].
# @param[in] e Eccentricity.
#
# @return The semi-major axis [m], the inclination [rad], the residues of the
# two functions and a boolean variable that indicates if the method converged.
#
#==#

function compute_ss_orbit_by_ang_vel(n::Real, e::Real)
    # Check if the arguments are valid.
    if (n <= 0)
        throw(ArgumentError("The angular velocity must be greater than 0."))
    end

    if !( 0. <= e < 1. )
        throw(ArgumentError("The eccentricity must be within the interval 0 <= e < 1."))
    end

    # Tolerance for the Newton-Raphson method.
    const tol = 1e-18

    # Maximum number of iterations for the Newton-Raphson method.
    const maxIt = 30

    # Auxiliary variables.
    sqrt_m0 = sqrt(m0)
    sqrt_e2 = sqrt(1-e^2)

    # Auxiliary constant to compute the functions.
    K1 = 3.0*R0^2*J2*sqrt_m0/(4.0*(1-e^2)^2)

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
    
    ################################################################################
    # Solve for zeros of f1 and f2 using Newton-Raphson method.
    ################################################################################

    # Maximum number of 
    
    # Initial guess based on the unperturbed model.
    a_k = (m0/n^2)^(1/3)
    i_k = acos( -ne*a_k^(3.5)/(2*K1) )

    # Loop
    it = 0;
    converged = true
    
    while (abs(f1(a_k, i_k)) > tol) || (abs(f2(a_k, i_k)) > tol)
        a_k_1 = a_k
        i_k_1 = i_k

        # Compute the Jacobian.
        J_k_1 = Jf1f2(a_k_1, i_k_1)

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

#==#
# 
# @brief Compute the sun-synchronous orbit given the inclination and the
# eccentricity
#
# @param[in] i Inclination [rad].
# @param[in] e Eccentricity.
#
# @return The semi-major axis [m].
#
#==#

function compute_ss_orbit_by_inclination(i::Real, e::Real)
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
        throw(ErrorException("It was not possible to find a valid sun-synchronous orbit with the inclination given."))
    end

    # Return.
    a
end

#==#
# 
# @brief Compute the sun-synchronous orbit given the number of revolutions per
# day and the eccentricity.
#
# @param[in] numRevPD Number of revolutions per day.
# @param[in] e Eccentricity.
#
# @return The semi-major axis [m], the inclination [rad], the residues of the
# two functions and a boolean variable that indicates if the method converged.
#
#==#

function compute_ss_orbit_by_num_rev_per_day(numRevPD::Real, e::Real)
    compute_ss_orbit_by_ang_vel(numRevPD*2*pi/86400.0, e)
end

#==#
# 
# @brief Compute the sun-synchronous orbit given the semi-major axis and the
# eccentricity
#
# @param[in] a Semi-major axis [m].
# @param[in] e Eccentricity.
#
# @return The inclination [rad].
#
#==#

function compute_ss_orbit_by_semi_major_axis(a::Real, e::Real)
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
        throw(ErrorException("It was not possible to find a sun-synchronous orbit with the semi-major axis given."))
    end
    
    # Return.
    acos(cos_i)
end

#==#
# 
# @brief Compute a list of repeating sun-synchronous orbits.
#
# @param[in] minRep Minimum repetition time of the orbit [days].
# @param[in] maxRep Maximum repetition time of the orbit [days].
# @param[in] minAlt Minimum altitude of the orbits on the list [m].
# @param[in] maxAlt Maximum altitude of the orbits on the list [m].
# @param[in] e Eccentricity.
#
# @return A matrix containing the orbits found.
#
# @remarks
#
# 1) If minAlt or maxAlt is < 0.0, then the altitude will not be checked when a
# orbit is added to the list.
#
# 2) The output matrix has the following format:
#
# Semi-major axis [m] | Altitude [m] | Period [s] | Int | Num | Den
# --------------------|--------------|------------|-----|-----|----
#
# in which the period is Int + Num/Den.
#
#==#

function list_ss_orbits_by_rep_period(minRep::Int,       maxRep::Int,
                                      minAlt::Real=-1.0, maxAlt::Real=-1.0,
                                      e::Real=0.0)
    numRevPD = 0
    # Check if the arguments are valid.
    if (minRep <= 0)
        throw(ArgumentError("The minimum repetition time must be greater than 0."))
    end

    if (maxRep <= 0)
        throw(ArgumentError("The maximum repetition time must be greater than 0."))
    end

    if (minRep > maxRep)
        throw(ArgumentError("The minimum repetition time must be smaller or equal to the maximum repetition time."))
    end
    
    if !( 0. <= e < 1. )
        throw(ArgumentError("The eccentricity must be within the interval 0 <= e < 1."))
    end
    
    # Matrix to store the available orbits.
    ss_orbits = Array(Float64, 0, 7)

    # Integer part of the number of orbits per day.
    intNumOrb = [13., 14., 15., 16., 17.]
    
    # Loop for the possible repetition times.
    for den = minRep:maxRep
        for num = 0:den-1
            # Check if the fraction num/den is irreducible.
            if ( gcd(num, den) == 1.0 )
                # Loop through the integer parts.
                for ino in intNumOrb
                    addOrbit = false

                    # Number of revolutions per day.
                    numRevPD = ino+float64(num)/float64(den)

                    (a, i, f1r, f2r, converged) =
                        compute_ss_orbit_by_num_rev_per_day(numRevPD, e)

                    # Check if the orbit is valid.
                    (!is_orbit_valid(a, e)) && continue
                        
                    # Check if the altitude interval must be verified.
                    if (minAlt > 0) && (maxAlt > 0)
                        if (minAlt < a-R0) && (a-R0 < maxAlt) && (converged)
                            addOrbit = true
                        end
                    else
                        addOrbit = converged
                    end

                    # Check if the orbit must be added to the list.
                    if (addOrbit)
                        # Compute the period of the orbit considering
                        # perturbations (J2).
                        period = t_J2(a, e, i)

                        ss_orbits =
                            vcat(ss_orbits, [a a-R0 i period ino num den])
                    end
                end
            end
        end
    end

    # Return the list of orbits.    
    ss_orbits
end

#==#
# 
# @brief Sort the list of sun-synchronous orbits by height.
#
# @param[in] ss_orbits List of sun-synchronous orbits (@see
# list_ss_orbits_by_rep_period).
#
# @return A matrix containing a list of sun-synchronous orbits sorted by height.
#
#==#

sort_list_ss_orbits_by_height(ss_orbits::Array{Float64, 2}) =
    sortrows(ss_orbits, by=x->x[1])
