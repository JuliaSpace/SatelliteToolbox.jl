# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Functions to convert osculating elements to mean elements using SGP4.
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

export rv_to_mean_elements_sgp4, rv_to_tle

"""
    rv_to_mean_elements_sgp4(vjd::AbstractVector{T}, vr::AbstractVector{Tv}, vv::AbstractVector{Tv}, W = I; estimate_bstar::Bool = true, mean_elements_epoch::Symbol = :end, sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84, print_debug::Bool = true, max_iterations::Int = 50, atol::Number = 2e-4, rtol::Number = 2e-4) where {T, Tv<:AbstractVector}

Compute the mean elements for the orbit propagator SGP4 based on the position
vectors `vr` and velocity vectors `vr`, both represented in TEME reference
frame. The epoch of those measurements, represented as Julian days, must be in
`vjd`.

The matrix `W` defined the weights for the least-square algorithm.

# Keywords

- `estimate_bstar::Bool`: If `true`, then the BSTAR parameters of the TLE will
    be estimated.
- `mean_elements_epoch::Symbol`: If it is  `:end`, the epoch of the mean
    elements will be equal to the last value in `vjd`. Otherwise, if it is
    `:begin`, the epoch will be selected as the first value in `vjd`.
- `sgp4_gc::SGP4_GravCte`: SPG4 gravitational constants (see `SGP4_GravCte`).
- `print_debug::Bool`: If `true`, then debug information will be printed to the
    `stdout`. (**Default** = `true`)
- `max_iterations::Int`: The maximum allowed number of iterations.
- `atol::Number`: The tolerance for the absolute value of the residue. If, at
    any iteration, the residue is lower than `atol`, then the iterations stop.
- `rtol::Number`: The tolerance for the relative difference between the
    residues. If, at any iteration, the relative difference between the residues
    in two consecutive iterations is lower than `rtol`, then the iterations stop.

# Returns

- The epoch of the elements [Julian Day].
- The mean elements for the SGP4 orbit propagator:
    - Semi-major axis [m];
    - Eccentricity [ ];
    - Inclination [rad];
    - Right ascension of the ascending node [rad];
    - Argument of perigee [rad];
    - True anomaly [rad];
    - BSTAR (0 if `estimate_bstar` is `false`).
- The covariance matrix of the mean elements estimation.
"""
function rv_to_mean_elements_sgp4(
    vjd::AbstractVector{T},
    vr::AbstractVector{Tv},
    vv::AbstractVector{Tv},
    W = I;
    estimate_bstar::Bool = true,
    mean_elements_epoch::Symbol = :end,
    sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84,
    print_debug::Bool = true,
    max_iterations::Int = 50,
    atol::Number = 2e-4,
    rtol::Number = 2e-4
) where {T, Tv<:AbstractVector}
    # Number of measurements.
    num_meas = length(vr)

    # Check if the orbit epoch must be the first or the last element.
    if mean_elements_epoch == :end
        r₁    = last(vr)
        v₁    = last(vv)
        epoch = last(vjd)
    else
        r₁    = first(vr)
        v₁    = first(vv)
        epoch = first(vjd)
    end

    # Initial guess of the mean elements.
    #
    # NOTE: x₁ is the previous estimate and x₂ is the current estimate.
    x₁ = estimate_bstar ?
        SVector{7, T}(r₁[1], r₁[2], r₁[3], v₁[1], v₁[2], v₁[3], T(0.00001)) :
        SVector{7, T}(r₁[1], r₁[2], r₁[3], v₁[1], v₁[2], v₁[3], T(0))
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
    print_debug && @printf("          %10s %20s %20s\n", "Iter.", "Residue", "Res. variation")

    # Loop until the maximum allowed iteration.
    @inbounds for it in 1:max_iterations
        x₁ = x₂

        # Variables to store the summations to compute the least square fitting
        # algorithm.
        ΣAᵀWA = @SMatrix zeros(num_states, num_states)
        ΣAᵀWb = @SVector zeros(num_states)

        # Variable to store the residue in this iteration.
        σ_i = T(0)

        @views for k in 1:num_meas
            # Obtain the measured ephemerides.
            y = vcat(vr[k], vv[k])

            # Obtain the computed ephemerides considering the current estimate
            # of the mean elements.
            Δt = (vjd[k] - epoch) * 86400

            r̂, v̂ = estimate_bstar ?
                _sgp4_sv(
                    Δt,
                    sgp4_gc,
                    epoch,
                    x₁[1],
                    x₁[2],
                    x₁[3],
                    x₁[4],
                    x₁[5],
                    x₁[6],
                    x₁[7]
                ) :
                _sgp4_sv(
                    Δt,
                    sgp4_gc,
                    epoch,
                    x₁[1],
                    x₁[2],
                    x₁[3],
                    x₁[4],
                    x₁[5],
                    x₁[6]
                )

            ŷ = vcat(r̂, v̂)

            # Compute the error.
            b = y - ŷ

            # Compute the Jacobian.
            A = _sgp4_jacobian(Δt, epoch, x₁, ŷ; estimate_bstar = estimate_bstar)

            # Accumulation.
            ΣAᵀWA += A' * W * A
            ΣAᵀWb += A' * W * b
            σ_i   += sum(W * (b.^2))
        end

        # Normalize the residue.
        σ_i /= num_meas

        # Update the estimate.
        P  = SMatrix{num_states, num_states, T}(pinv(ΣAᵀWA))
        δx = P * ΣAᵀWb

        # Limit the correction to avoid divergence, but it should not be applied
        # to B*.
        for i in 1:6
            threshold = 0.1
            if abs(δx[i] / x₁[i]) > threshold
                δx = setindex(δx, threshold * abs(x₁[i]) * sign(δx[i]), i)
            end
        end

        x₂ = x₁ + δx

        # Compute the residue variation.
        σ_p = (σ_i - σ_i₋₁) / σ_i₋₁

        if print_debug
            if it != 1
                @printf("PROGRESS: %10d %20g %20g %%\n", it, σ_i, 100σ_p)
            else
                @printf("PROGRESS: %10d %20g %20s\n", it, σ_i, "---")
            end
        end

        # Check if the residue is increasing.
        if σ_i < σ_i₋₁
            Δd = 0
        else
            Δd += 1
        end

        # If the residue increased by three iterations and the residue is higher
        # than 5e11, then we abort because the iterations are diverging.
        ((Δd ≥ 3) && (σ_i > 5e11)) && error("The iterations diverged!")

        # Check if the condition to stop has been reached.
        ((abs(σ_p) < rtol) || (σ_i < atol) || (it ≥ max_iterations)) && break

        σ_i₋₁ = σ_i
    end

    # Convert the state vector to mean elements.
    orb = rv_to_kepler(x₂[1], x₂[2], x₂[3], x₂[4], x₂[5], x₂[6])
    bstar = x₂[7]

    # Assemble the output vector.
    xo = @SVector [orb.a, orb.e, orb.i, orb.Ω, orb.ω, orb.f, bstar]

    # Return the mean elements for SGP4 and the covariance matrix.
    return epoch, xo, P
end

"""
    rv_to_tle(args...; name::String = "UNDEFINED", sat_num::Int = 9999, classification::Char = 'U', int_designator = "999999", elem_set_number::Int = 0, rev_num, kwargs...)

Convert a set of position and velocity vectors represented in TEME reference
frame to a TLE. The arguments `args` and keywords `kwargs` are the same as those
described in the function `rv_to_mean_elements_sgp4`.

Additionally, the user can specify some parameters of the generated TLE.

This function returns the TLE and the covariance of the estimated elements
(state vector).
"""
function rv_to_tle(
    args...;
    name::String = "UNDEFINED",
    sat_num::Int = 9999,
    classification::Char = 'U',
    int_designator = "999999",
    elem_set_number::Int = 0,
    rev_num::Int = 0,
    sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84,
    kwargs...
)
    # Convert the position and velocity vectors to mean elements.
    JD, x, P = rv_to_mean_elements_sgp4(args...; kwargs...)

    # Compute the data as required by the TLE format.
    dt  = jd_to_date(DateTime, JD)
    dt₀ = jd_to_date(DateTime, date_to_jd(year(dt), 1, 1, 0, 0, 0))

    dt_year    = year(dt)
    epoch_year = dt_year < 1980 ? dt_year - 1900 : dt_year - 2000
    epoch_day  = (dt - dt₀).value / 1000 / 86400 + 1

    # Obtain the Keplerian elements with the right units.
    a₀ = x[1] / (1000sgp4_gc.R0)
    e₀ = x[2]
    i₀ = rad2deg(x[3])
    Ω₀ = rad2deg(x[4])
    ω₀ = rad2deg(x[5])
    M₀ = rad2deg(f_to_M(e₀, x[6]))

    # Obtain the mean motion [rad/min].
    n₀ = sgp4_gc.XKE / sqrt(a₀ * a₀ * a₀)

    # Construct the TLE.
    tle = TLE(
        name,
        sat_num,
        classification,
        int_designator,
        epoch_year,
        epoch_day,
        JD,
        0.0,
        0.0,
        x[7],
        elem_set_number,
        0,
        i₀,
        Ω₀,
        e₀,
        ω₀,
        M₀,
        720n₀ / π,
        rev_num,
        0
    )

    # Return the TLE string.
    return tle, P
end


################################################################################
#                              Private functions
################################################################################

# Compute the SGP4 algorithm considering all variables in a state vector.
function _sgp4_sv(
    Δt::Number,
    sgp4_gc::SGP4_GravCte,
    epoch::Number,
    rx_TEME::Number,
    ry_TEME::Number,
    rz_TEME::Number,
    vx_TEME::Number,
    vy_TEME::Number,
    vz_TEME::Number,
    bstar::Number = 0
)
    r_TEME = @SVector [rx_TEME, ry_TEME, rz_TEME]
    v_TEME = @SVector [vx_TEME, vy_TEME, vz_TEME]

    orb_TEME = rv_to_kepler(r_TEME, v_TEME, epoch)

    # Obtain the required mean elements to initialize the SGP4.
    a₀ = orb_TEME.a / (1000sgp4_gc.R0) # .................. Semi-major axis [ER]
    e₀ = orb_TEME.e                    # ...................... Eccentricity [ ]
    i₀ = orb_TEME.i                    # ..................... Inclination [rad]
    Ω₀ = orb_TEME.Ω                    # ............................ RAAN [rad]
    ω₀ = orb_TEME.ω                    # ................. Arg. of perigee [rad]
    M₀ = f_to_M(e₀, orb_TEME.f)        # .................... Mean anomaly [rad]

    # Obtain the mean motion [rad/min].
    n₀ = sgp4_gc.XKE / sqrt(a₀ * a₀ * a₀)

    r, v, ~ = sgp4(Δt / 60, sgp4_gc, epoch, n₀, e₀, i₀, Ω₀, ω₀, M₀, bstar)

    # Return the elements using SI units.
    return 1000r, 1000v
end

function _sgp4_jacobian(
    Δt::T,
    epoch::T,
    x₁::SVector{NS, T},
    y₁::SVector{NO, T};
    estimate_bstar::Bool = true,
    pert::T = 1e-3,
    perttol::T = 1e-5,
    sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84
) where {NS, NO, T}
    num_states = NS
    dim_obs    = NO
    M          = zeros(T, dim_obs, num_states)

    # Auxiliary variables.
    x₂ = copy(x₁)

    # If B* does not need to be estimated, then does not compute the last
    # column.
    J = num_states - !estimate_bstar

    @inbounds for j in 1:J
        # State that will be perturbed.
        α = x₂[j]

        # Obtain the perturbation, taking care to avoid small values.
        ϵ = T(0)
        pert_i = pert

        for _ in 1:5
            ϵ = α * pert_i
            abs(ϵ) > perttol && break
            pert_i *= 1.4
        end

        α += ϵ

        # Avoid division by zero in cases that α is very small. In this
        # situation, we force |α| = perttol.
        if abs(α) < perttol
            α = sign(α) * perttol
        end

        # Modify the perturbed state.
        x₂ = setindex(x₂, α, j)

        # Obtain the Jacobian by finite differentiation.
        r, v = _sgp4_sv(Δt, sgp4_gc, epoch, x₂...)
        y₂   = vcat(r,v)

        M[:, j] .= (y₂ .- y₁) ./ ϵ

        # Restore the value of the perturbed state for the next iteration.
        x₂ = setindex(x₂, x₁[j], j)
    end

    return M
end
