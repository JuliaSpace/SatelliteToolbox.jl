#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Deprecation warnings.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

#                     Introduced in SatelliteToolbox v0.4
# ==============================================================================

@deprecate parse_gfc(filename) create_gravity_model_coefs(parse_icgem(filename))
@deprecate satellite_orbit_compute_f(a, e, i, M) M_to_f(e, M)
@deprecate satellite_orbit_compute_f(a, e, i, M, tol) M_to_f(e, M, tol)

#                     Introduced in SatelliteToolbox v0.5
# ==============================================================================

@deprecate legendre_fully_normalized(ϕ::Number, n_max::Number, ph_term::Bool)         legendre_fully_normalized(ϕ, n_max, -1, ph_term)
@deprecate legendre_schmidt_quasi_normalized(ϕ::Number, n_max::Number, ph_term::Bool) legendre_schmidt_quasi_normalized(ϕ, n_max, -1, ph_term)
@deprecate legendre_conventional(ϕ::Number, n_max::Number, ph_term::Bool)             legendre_conventional(ϕ, n_max, -1, ph_term)
@deprecate legendre(::Type{Val{:conv}}, ϕ::Number, n_max::Number, ph_term::Bool)      legendre(Val{:conv}, ϕ, n_max, -1, ph_term)
@deprecate legendre(::Type{Val{:schmidt}}, ϕ::Number, n_max::Number, ph_term::Bool)   legendre(Val{:schmidt}, ϕ, n_max, -1, ph_term)
@deprecate legendre(::Type{Val{:full}}, ϕ::Number, n_max::Number, ph_term::Bool)      legendre(Val{:full}, ϕ, n_max, -1, ph_term)
@deprecate legendre(ϕ::Number, n_max::Number, ph_term::Bool)                          legendre(ϕ, n_max, -1, ph_term)

@deprecate dlegendre_fully_normalized(ϕ::Number, n_max::Number, ph_term::Bool)         dlegendre_fully_normalized(ϕ, n_max, -1, ph_term)
@deprecate dlegendre_schmidt_quasi_normalized(ϕ::Number, n_max::Number, ph_term::Bool) dlegendre_schmidt_quasi_normalized(ϕ, n_max, -1, ph_term)
@deprecate dlegendre_conventional(ϕ::Number, n_max::Number, ph_term::Bool)             dlegendre_conventional(ϕ, n_max, -1, ph_term)
@deprecate dlegendre(::Type{Val{:conv}}, ϕ::Number, n_max::Number, ph_term::Bool)      dlegendre(Val{:conv}, ϕ, n_max, -1, ph_term)
@deprecate dlegendre(::Type{Val{:schmidt}}, ϕ::Number, n_max::Number, ph_term::Bool)   dlegendre(Val{:schmidt}, ϕ, n_max, -1, ph_term)
@deprecate dlegendre(::Type{Val{:full}}, ϕ::Number, n_max::Number, ph_term::Bool)      dlegendre(Val{:full}, ϕ, n_max, -1, ph_term)
@deprecate dlegendre(ϕ::Number, n_max::Number, ph_term::Bool)                          dlegendre(ϕ, n_max, -1, ph_term)

#                     Introduced in SatelliteToolbox v0.6
# ==============================================================================

@deprecate satellite_beta_angle(JD0, a, e, i, RAAN, numDays)                beta_angle(JD0, a, e, i, RAAN, numDays)
@deprecate satellite_check_station(r_e, rs_e, minElev)                      ground_station_visible(r_e, rs_e, minElev)
@deprecate satellite_check_station(r_e, lat_s, lon_s, h_s, minElev)         ground_station_visible(r_e, lat_s, lon_s, h_s, minElev)
@deprecate satellite_eclipse_time(JD0, a, e, i, w, RAAN, numDays)           eclipse_time_summary(JD0, a, e, i, RAAN, w, numDays)
@deprecate satellite_eclipse_time(JD0, a, e, i, w, RAAN, numDays, relative) eclipse_time_summary(JD0, a, e, i, RAAN, w, numDays, relative)

# Deprecation warnings related to Val{T} => Val(T)
# ==============================================================================

# IGRF
@deprecate igrf12(date, r, λ, Ω, ::Type{Val{T}}; show_warns = true) where T igrf12(date, r, λ, Ω, Val(T); show_warns = show_warns)

# Legendre
@deprecate dlegendre!(::Type{Val{T}}, dP, ϕ, P, ph_term = false) where T            dlegendre!(Val(T), P, ϕ, dP, ph_term)
@deprecate dlegendre(::Type{Val{T}}, ϕ, n_max, m_max = -1, ph_term = false) where T dlegendre(Val(T), ϕ, n_max, m_max, ph_term)
@deprecate legendre!(::Type{Val{T}}, P, ϕ, ph_term = false) where T                 legendre!(Val(T), P, ϕ, ph_term)
@deprecate legendre(::Type{Val{T}}, ϕ, n_max, m_max = -1, ph_term = false) where T  legendre(Val(T), ϕ, n_max, m_max, ph_term)

# Space indices
for sym in (:F10, :F10obs, :F10adj, :Kp, :Ap, :Kp_vect, :Ap_vect, :S10, :S81a,
            :M10, :M81a, :Y10, :Y81a, :DstΔTc)
    qsym = Meta.quot(sym)
    @eval @deprecate get_space_index(::Type{Val{$qsym}}, JD) get_space_index(Val($qsym), JD)
end

for sym in (:F10M, :F10Mobs, :F10Madj)
    qsym = Meta.quot(sym)
    @eval @deprecate get_space_index(::Type{Val{$qsym}}, JD; window = 81) get_space_index(Val($qsym), JD; window = window)
end
