# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions to convert osculating elements to mean elements using the J2
#   osculating algorithm.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A; Crawford, P (2008). SGP4 orbit determination. AIAA/AAS
#       Astrodynamics Specialist Conference, Honoulu, HI.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export rv_to_mean_elements_j2osc

"""
    rv_to_mean_elements_j2osc(vJD::AbstractVector{T}, vr::AbstractVector{Tv}, vv::AbstractVector{Tv}, W = I; mean_elements_epoch::Symbol = :end, max_it::Int = 50, j2_gc = sgp4_gc_wgs84, atol::Number = 2e-4, rtol::Number = 2e-4) where {T,Tv<:AbstractVector}

Compute the mean elements for SGP4 based on the position `vr` and velocity
vectors `vr` represented in TEME reference frame. The epoch of those
measurements [Julian Day] must be in `vJD`.

The matrix `W` defined the weights for the least-square algorithm.

# Keywords

* `estimate_bstar`: If `true`, then the BSTAR parameters of the TLE will be
                    estimated.
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
    * True anomaly [rad];
    * BSTAR (0 if `estimate_bstar` is `false`).
* The covariance matrix of the mean elements estimation.

"""
function rv_to_mean_elements_j2osc(vJD::AbstractVector{T},
                                   vr::AbstractVector{Tv},
                                   vv::AbstractVector{Tv},
                                   W = I;
                                   mean_elements_epoch::Symbol = :end,
                                   max_it::Int = 50,
                                   j2_gc = j2_gc_egm08,
                                   atol::Number = 2e-4,
                                   rtol::Number = 2e-4) where {T, Tv<:AbstractVector}

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

        @views for k = 1:num_meas
            # Obtain the measured ephemerides.
            y = vcat(vr[k], vv[k])

            # Obtain the computed ephemerides considering the current estimate
            # of the mean elements.
            Δt = (vJD[k] - epoch)*86400

            r̂, v̂ = _j2osc_sv(Δt, j2_gc, epoch,
                             x₁[1], x₁[2], x₁[3],
                             x₁[4], x₁[5], x₁[6])

            ŷ = vcat(r̂, v̂)

            # Compute the error.
            b = y - ŷ

            # Compute the Jacobian.
            A = _j2osc_jacobian(Δt, epoch, x₁, ŷ)

            # Accumulation.
            ΣAᵀWA += A'*W*A
            ΣAᵀWb += A'*W*b
            σ_i   += sum(W*(b.^2))
        end

        # Normalize the residue.
        σ_i /= num_meas

        # Update the estimate.
        P  = SMatrix{num_states, num_states, T}(pinv(ΣAᵀWA))
        δx = P*ΣAᵀWb

        # Limit the correction to avoid divergence.
        for i = 1:num_states
            abs(δx[i] / x₁[i]) > 0.01 &&
                (δx = setindex(δx, 0.01 * abs(x₁[i]) * sign(δx[i]), i))
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
    orb = rv_to_kepler(x₂[1], x₂[2], x₂[3], x₂[4], x₂[5], x₂[6], epoch)

    # Return the mean elements for J2 osculating algorithm and the covariance
    # matrix.
    return orb, P
end

################################################################################
#                              Private functions
################################################################################

# Compute the J2 osculating algorithm considering all variables in a state
# vector.
function _j2osc_sv(Δt::Number,
                   j2_gc::J2_GravCte,
                   epoch::Number,
                   rx_i::Number,
                   ry_i::Number,
                   rz_i::Number,
                   vx_i::Number,
                   vy_i::Number,
                   vz_i::Number)

    r_i = @SVector [rx_i, ry_i, rz_i]
    v_i = @SVector [vx_i, vy_i, vz_i]

    # Obtain the required mean elements to initialize the J2 osculating
    # propagator.
    orb_i = rv_to_kepler(r_i, v_i, epoch)

    j2oscd = j2osc_init(orb_i.t, orb_i.a, orb_i.e, orb_i.i, orb_i.Ω, orb_i.ω,
                        orb_i.f, 0, 0; j2_gc = j2_gc)
    r, v = j2osc!(j2oscd, Δt)

    # Return the elements using SI units.
    return r, v
end

function _j2osc_jacobian(Δt::T,
                         epoch::T,
                         x₁::SVector{NS, T},
                         y₁::SVector{NO, T};
                         pert::T = 1e-3,
                         perttol::T = 1e-5,
                         j2_gc::J2_GravCte = j2_gc_egm08) where {NS, NO, T}

    num_states = NS
    dim_obs    = NO
    M          = zeros(T, dim_obs, num_states)

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
        r, v = _j2osc_sv(Δt, j2_gc, epoch, x₂...)
        y₂   = vcat(r, v)

        M[:,j] .= (y₂ .- y₁)./ϵ

        # Restore the value of the perturbed state for the next iteration.
        x₂ = setindex(x₂, x₁[j], j)
    end

    return M
end
