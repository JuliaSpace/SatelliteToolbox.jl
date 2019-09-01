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

@deprecate satellite_check_station(r_e, rs_e, minElev)              ground_station_visible(r_e, rs_e, minElev)
@deprecate satellite_check_station(r_e, lat_s, lon_s, h_s, minElev) ground_station_visible(r_e, lat_s, lon_s, h_s, minElev)
