# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions to convert osculating elements to mean elements.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A; Crawford, P (2008). SGP4 orbit determination. AIAA/AAS
#       Astrodynamics Specialist Conference, Honoulu, HI.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export rv_to_mean_elements_sgp4, rv_to_tle

"""
    rv_to_mean_elements_sgp4(vJD::AbstractVector{T}, vr::AbstractVector{Tv}, vv::AbstractVector{Tv}, W = I; mean_elements_epoch::Symbol = :end, max_it::Int = 25, sgp4_gc = sgp4_gc_wgs84, atol::Number = 2e-4, rtol::Number = 2e-4) where

Compute the mean elements for SGP4 based on the position `vr` and velocity
vectors `vr` represented in TEME reference frame. The epoch of those
measurements [Julian Day] must be in `vJD`.

The matrix `W` defined the weights for the least-square algorithm.

# Keywords

* `mean_elements_epoch`: If it is  `:end`, the epoch of the mean elements will
                         be equal to the last value in `vJD`. Otherwise, if it
                         is `:begin`, the epoch will be selected as the first
                         value in `vJD`.
* `max_it`: The maximum allowed number of iterations.
* `sgp4_gc`: SPG4 constants (see `SGP4_GravCte`).
* `atol`: The tolerance for the absolute value of the residue. If, at any
          iteration, the residue is lower than `atol`, then the iterations stop.
* `rtol`: The tolerance for the relative difference between the residues. If, at
          any iteration, the relative difference between the residues in two
          consecutive iterations is lower than `rtol`, then the iterations stop.

# Returns

* The epoch of the elements [Julian Day].
* The mean elements for SGP4 algorithm:
    * Semi-major axis [m];
    * Eccentricity [ ];
    * Inclination [rad];
    * Right ascension of the ascending node [rad];
    * Argument of perigee [rad];
    * True anomaly [rad].
* The covariance matrix of the mean elements estimation.

"""
function rv_to_mean_elements_sgp4(vJD::AbstractVector{T},
                                  vr::AbstractVector{Tv},
                                  vv::AbstractVector{Tv},
                                  W = I;
                                  mean_elements_epoch::Symbol = :end,
                                  max_it::Int = 25,
                                  sgp4_gc = sgp4_gc_wgs84,
                                  atol::Number = 2e-4,
                                  rtol::Number = 2e-4) where
    {T,Tv<:AbstractVector}

    # Number of measurements.
    num_meas = length(vr)

    # Check if the orbit epoch must be the first or the last element.
    if mean_elements_epoch == :end
        r₁ = last(vr)
        v₁ = last(vv)
        epoch = last(vJD)
    else
        r₁ = first(vr)
        v₁ = first(vv)
        epoch = first(vJD)
    end

    # Initial guess of the mean elements.
    #
    # NOTE: x₁ is the previous estimate and x₂ is the current estimate.
    x₁ = SVector{6,T}(r₁[1], r₁[2], r₁[3], v₁[1], v₁[2], v₁[3])
    x₂ = x₁

    # Number of states in the input vector.
    num_states = length(x₁)

    # Covariance matrix.
    P = SMatrix{num_states, num_states, T}(I)

    # Variable to store the last residue.
    σ_i₋₁ = T(NaN)

    # Variable to store how many iterations the residue increased. This is used
    # to account for divergence.
    Δd = 0

    # Header
    @printf("          %10s %20s %20s\n", "Iter.", "Residue", "Res. variation")

    # Loop until the maximum allowed iteration.
    @inbounds for it = 1:max_it
        x₁ = x₂

        # Variables to store the summations to compute the least square fitting
        # algorithm.
        ΣAᵀWA = @SMatrix zeros(num_states, num_states)
        ΣAᵀWb = @SVector zeros(num_states)

        # Variable to store the residue in this iteration.
        σ_i = T(0)

        for k = 1:num_meas
            # Obtain the measured ephemerides.
            y = vcat(vr[k], vv[k])

            # Obtain the computed ephemerides considering the current estimate
            # of the mean elements.
            Δt = (vJD[k] - epoch)*86400
            r̂, v̂, ~ = _sgp4_sv(Δt, sgp4_gc, epoch, x₁...)
            ŷ = vcat(r̂, v̂)

            # Compute the error.
            b = y - ŷ

            # Compute the Jacobian.
            A = _sgp4_jacobian(Δt, epoch, x₁, ŷ)

            # Accumulation.
            ΣAᵀWA += A'*W*A
            ΣAᵀWb += A'*W*b
            σ_i   += sum(W*(b.^2))
        end

        # Normalize the residue.
        σ_i /= num_meas

        # Update the estimate.
        P  = pinv(ΣAᵀWA)
        δx = P*ΣAᵀWb

        # Limit the correction to avoid divergence.
        for i = 1:num_states
            abs(δx[i] / x₁[i]) > 0.01 && (δx[i] = 0.01 * abs(x₁[i]) * sign(δx[i]))
        end

        x₂ = x₁ + δx

        # Compute the residue variation.
        σ_p = (σ_i - σ_i₋₁)/σ_i₋₁

        if it != 1
            @printf("PROGRESS: %10d %20g %20g %%\n", it, σ_i, 100σ_p)
        else
            @printf("PROGRESS: %10d %20g %20s\n", it, σ_i, "---")
        end

        # Check if the residue is increasing.
        if σ_i < σ_i₋₁
            Δd = 0
        else
            Δd += 1
        end

        # If the residue increased by three iterations and the residue is higher
        # than 5e8, then we abort because the iterations are diverging.
        ( (Δd ≥ 3) && (σ_i > 5e11) ) && error("The iterations diverged!")

        # Check if the condition to stop has been reached.
        ((abs(σ_p) < rtol) || (σ_i < atol) || (it ≥ max_it)) && break

        σ_i₋₁ = σ_i
    end

    # Convert the state vector to mean elements.
    orb = rv_to_kepler((@view x₂[1:6])...)

    # Assemble the output vector.
    xo = @SVector [orb.a, orb.e, orb.i, orb.Ω, orb.ω, orb.f]

    # Return the mean elements for SGP4 and the covariance matrix.
    return epoch, xo, P
end

"""
    rv_to_tle(args...; name::String = "UNDEFINED", sat_num::Int = 9999, classification::Char = 'U', int_designator = "999999", elem_set_number::Int = 0, rev_num, kwargs...)

Convert a set of position and velocity vectors represented in TEME reference
frame to a TLE. The arguments `args` and keywords `kwargs` are the same as those
described in the function `rv_to_mean_elements_sgp4`.

Additionally, the user can specify some parameters of the generated TLE.

"""
function rv_to_tle(args...;
                   name::String = "UNDEFINED",
                   sat_num::Int = 9999,
                   classification::Char = 'U',
                   int_designator = "999999",
                   elem_set_number::Int = 0,
                   rev_num::Int = 0,
                   sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84,
                   kwargs...)

    # Convert the position and velocity vectors to mean elements.
    JD, x, P = rv_to_mean_elements_sgp4(args...; kwargs...)

    # Compute the data as required by the TLE format.
    dt  = JDtoDate(DateTime, JD)
    dt₀ = JDtoDate(DateTime, DatetoJD(year(dt), 1, 1, 0, 0, 0))

    dt_year    = year(dt)
    epoch_year = dt_year < 1980 ? dt_year - 1900 : dt_year - 2000
    epoch_day  = (dt - dt₀).value/1000/86400 + 1

    # Obtain the Keplerian elements with the right units.
    a₀ = x[1]/(1000sgp4_gc.R0)
    e₀ = x[2]
    i₀ = rad2deg(x[3])
    Ω₀ = rad2deg(x[4])
    ω₀ = rad2deg(x[5])
    M₀ = rad2deg(f_to_M(e₀, x[6]))

    # Obtain the mean motion [rad/min].
    n₀ = sgp4_gc.XKE / sqrt(a₀ * a₀ * a₀)

    # Construct the TLE.
    tle = TLE(name,
              sat_num,
              classification,
              int_designator,
              epoch_year,
              epoch_day,
              JD,
              0.0,
              0.0,
              0.0, # x[7],
              elem_set_number,
              0,
              i₀,
              Ω₀,
              e₀,
              ω₀,
              M₀,
              720n₀/π,
              rev_num,
              0)

    # Return the TLE string.
    return tle
end


################################################################################
#                              Private functions
################################################################################

# Compute the SGP4 algorithm considering all variables in a state vector.
function _sgp4_sv(Δt::Number,
                  sgp4_gc::SGP4_GravCte,
                  epoch::Number,
                  rx_TEME::Number,
                  ry_TEME::Number,
                  rz_TEME::Number,
                  vx_TEME::Number,
                  vy_TEME::Number,
                  vz_TEME::Number,
                  bstar::Number = 0)

    r_TEME = @SVector [rx_TEME, ry_TEME, rz_TEME]
    v_TEME = @SVector [vx_TEME, vy_TEME, vz_TEME]

    orb_TEME = rv_to_kepler(r_TEME, v_TEME, epoch)

    # Obtain the required mean elements to initialize the SGP4.
    a₀ = orb_TEME.a/(1000sgp4_gc.R0) # .................... Semi-major axis [ER]
    e₀ = orb_TEME.e                  # ........................ Eccentricity [ ]
    i₀ = orb_TEME.i                  # ....................... Inclination [rad]
    Ω₀ = orb_TEME.Ω                  # .............................. RAAN [rad]
    ω₀ = orb_TEME.ω                  # ................... Arg. of perigee [rad]
    M₀ = f_to_M(e₀, orb_TEME.f)      # ...................... Mean anomaly [rad]

    # Obtain the mean motion [rad/min].
    n₀ = sgp4_gc.XKE / sqrt(a₀ * a₀ * a₀)

    r, v, sgp4d = sgp4(Δt/60, sgp4_gc, epoch, n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar)

    # Return the elements using SI units.
    return 1000r, 1000v, sgp4d
end

function _sgp4_jacobian(Δt::T,
                        epoch::T,
                        x₁::SVector{NS, T},
                        y₁::SVector{NO, T};
                        pert::T = 1e-3,
                        perttol::T = 1e-5,
                        sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84) where {NS, NO, T}

    num_states = NS
    dim_obs    = NO
    M          = Matrix{T}(undef, dim_obs, num_states)

    # Auxiliary variables.
    x₂ = copy(x₁)

    @inbounds for j = 1:num_states
        # State that will be perturbed.
        α = x₂[j]

        # Obtain the perturbation, taking care to avoid small values.
        ϵ = T(0)
        pert_i = pert

        for _ = 1:5
            ϵ = α * pert_i
            abs(ϵ) > perttol && break
            pert_i *= 1.4
        end

        α += ϵ

        # Avoid division by zero in cases that α is very small. In this
        # situation, we force |α| = perttol.
        abs(α) < perttol && (α = sign(α) * perttol)

        # Modify the perturbed state.
        x₂ = setindex(x₂, α, j)

        # Obtain the Jacobian by finite differentiation.
        r, v, ~ = _sgp4_sv(Δt, sgp4_gc, epoch, x₂...)
        y₂      = vcat(r,v)

        M[:,j] .= (y₂ .- y₁)./ϵ

        # Restore the value of the perturbed state for the next iteration.
        x₂ = setindex(x₂, x₁[j], j)
    end

    return M
end
