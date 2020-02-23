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
# ------------------------------------------------------------------------------

@deprecate igrf12(date, r, λ, Ω, ::Type{Val{T}}; show_warns = true) where T igrf12(date, r, λ, Ω, Val(T); show_warns = show_warns)

# Legendre
# ------------------------------------------------------------------------------

@deprecate dlegendre!(::Type{Val{T}}, dP, ϕ, P, ph_term = false) where T            dlegendre!(Val(T), P, ϕ, dP, ph_term)
@deprecate dlegendre(::Type{Val{T}}, ϕ, n_max, m_max = -1, ph_term = false) where T dlegendre(Val(T), ϕ, n_max, m_max, ph_term)
@deprecate legendre!(::Type{Val{T}}, P, ϕ, ph_term = false) where T                 legendre!(Val(T), P, ϕ, ph_term)
@deprecate legendre(::Type{Val{T}}, ϕ, n_max, m_max = -1, ph_term = false) where T  legendre(Val(T), ϕ, n_max, m_max, ph_term)

# Orbit propagation
# ------------------------------------------------------------------------------

# Individual parameters

@deprecate init_orbit_propagator(::Type{Val{:J2}}, epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, dn_o2 = 0, ddn_o6 = 0, j2_gc = j2_gc_egm08) where T #=
=#  init_orbit_propagator(Val(:J2), epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, dn_o2, ddn_o6, j2_gc)

@deprecate init_orbit_propagator(::Type{Val{:J4}}, epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, dn_o2 = 0, ddn_o6 = 0, j4_gc = j4_gc_egm08) where T #=
=#  init_orbit_propagator(Val(:J4), epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, dn_o2, ddn_o6, j4_gc)

@deprecate init_orbit_propagator(::Type{Val{:twobody}}, epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, μ = m0) where T #=
=#  init_orbit_propagator(Val(:twobody), epoch, a_0, e_0, i_0, Ω_0, ω_0, f_0, μ)

# Orbit structure

@deprecate init_orbit_propagator(::Type{Val{:J2}}, orb_0::Orbit, dn_o2::Number = 0, ddn_o6::Number = 0, j2_gc::J2_GravCte = j2_gc_egm08) #=
=#  init_orbit_propagator(Val(:J2), orb_0, dn_o2, ddn_o6, j2_gc)

@deprecate init_orbit_propagator(::Type{Val{:J4}}, orb_0::Orbit, dn_o2::Number = 0, ddn_o6::Number = 0, j4_gc::J4_GravCte = j4_gc_egm08) #=
=#  init_orbit_propagator(Val(:J4), orb_0, dn_o2, ddn_o6, j4_gc)

@deprecate init_orbit_propagator(::Type{Val{:twobody}}, orb_0::Orbit, μ::Number = m0) #=
=#  init_orbit_propagator(Val(:twobody), orb_0, μ)

# TLE

@deprecate init_orbit_propagator(::Type{Val{:J2}}, tle::TLE, j2_gc::J2_GravCte = j2_gc_egm08) #=
=#  init_orbit_propagator(Val(:J2), tle, j2_gc)

@deprecate init_orbit_propagator(::Type{Val{:J4}}, tle::TLE, j4_gc::J4_GravCte = j4_gc_egm08) #=
=#  init_orbit_propagator(Val(:J4), tle, j4_gc)

@deprecate init_orbit_propagator(::Type{Val{:twobody}}, tle::TLE, μ = m0) #=
=#  init_orbit_propagator(Val(:twobody), tle, μ)

@deprecate init_orbit_propagator(::Type{Val{:sgp4}}, tle::TLE, sgp4_gc::SGP4_GravCte = sgp4_gc_wgs84) #=
=#  init_orbit_propagator(Val(:sgp4), tle, sgp4_gc)

# Space indices
# ------------------------------------------------------------------------------

for sym in (:F10, :F10obs, :F10adj, :Kp, :Ap, :Kp_vect, :Ap_vect, :S10, :S81a,
            :M10, :M81a, :Y10, :Y81a, :DstΔTc)
    qsym = Meta.quot(sym)
    @eval @deprecate get_space_index(::Type{Val{$qsym}}, JD) get_space_index(Val($qsym), JD)
end

for sym in (:F10M, :F10Mobs, :F10Madj)
    qsym = Meta.quot(sym)
    @eval @deprecate get_space_index(::Type{Val{$qsym}}, JD; window = 81) get_space_index(Val($qsym), JD; window = window)
end

# Transformations
# ------------------------------------------------------------------------------

# ECEF to ECEF

@deprecate rECEFtoECEF(::Type{Val{To}}, ::Type{Val{Tf}}, JD, eop_data) where {To,Tf}    rECEFtoECEF(Val(To), Val(Tf), JD, eop_data)
@deprecate rECEFtoECEF(T, ::Type{Val{To}}, ::Type{Val{Tf}}, JD, eop_data) where {To,Tf} rECEFtoECEF(T, Val(To), Val(Tf), JD, eop_data)

# ECEF to ECI

@deprecate rECEFtoECI(::Type{Val{To}}, ::Type{Val{Tf}}, JD) where {To,Tf}              rECEFtoECI(Val(To), Val(Tf), JD)
@deprecate rECEFtoECI(::Type{Val{To}}, ::Type{Val{Tf}}, JD, eop_data) where {To,Tf}    rECEFtoECI(Val(To), Val(Tf), JD, eop_data)
@deprecate rECEFtoECI(T, ::Type{Val{To}}, ::Type{Val{Tf}}, JD) where {To,Tf}           rECEFtoECI(T, Val(To), Val(Tf), JD)
@deprecate rECEFtoECI(T, ::Type{Val{To}}, ::Type{Val{Tf}}, JD, eop_data) where {To,Tf} rECEFtoECI(T, Val(To), Val(Tf), JD, eop_data)

# ECI to ECEF

@deprecate rECItoECEF(::Type{Val{To}}, ::Type{Val{Tf}}, JD) where {To,Tf}              rECItoECEF(Val(To), Val(Tf), JD)
@deprecate rECItoECEF(::Type{Val{To}}, ::Type{Val{Tf}}, JD, eop_data) where {To,Tf}    rECItoECEF(Val(To), Val(Tf), JD, eop_data)
@deprecate rECItoECEF(T, ::Type{Val{To}}, ::Type{Val{Tf}}, JD) where {To,Tf}           rECItoECEF(T, Val(To), Val(Tf), JD)
@deprecate rECItoECEF(T, ::Type{Val{To}}, ::Type{Val{Tf}}, JD, eop_data) where {To,Tf} rECItoECEF(T, Val(To), Val(Tf), JD, eop_data)

# ECI to ECI

@deprecate rECItoECI(::Type{Val{To}}, ::Type{Val{Tf}}, JD) where {To,Tf}                    rECItoECI(Val(To), Val(Tf), JD)
@deprecate rECItoECI(::Type{Val{To}}, JDo, ::Type{Val{Tf}}, JDf) where {To,Tf}              rECItoECI(Val(To), JDo, Val(Tf), JDf)
@deprecate rECItoECI(::Type{Val{To}}, ::Type{Val{Tf}}, JD, eop_data) where {To,Tf}          rECItoECI(Val(To), Val(Tf), JD, eop_data)
@deprecate rECItoECI(::Type{Val{To}}, JDo, ::Type{Val{Tf}}, JDf, eop_data) where {To,Tf}    rECItoECI(Val(To), JDo, Val(Tf), JDf, eop_data)
@deprecate rECItoECI(T::T_ROT, ::Type{Val{To}}, ::Type{Val{Tf}}, JD) where {To,Tf}                 rECItoECI(T, Val(To), Val(Tf), JD)
@deprecate rECItoECI(T::T_ROT, ::Type{Val{To}}, JDo, ::Type{Val{Tf}}, JDf) where {To,Tf}           rECItoECI(T, Val(To), JDo, Val(Tf), JDf)
@deprecate rECItoECI(T::T_ROT, ::Type{Val{To}}, ::Type{Val{Tf}}, JD, eop_data) where {To,Tf}       rECItoECI(T, Val(To), Val(Tf), JD, eop_data)
@deprecate rECItoECI(T::T_ROT, ::Type{Val{To}}, JDo, ::Type{Val{Tf}}, JDf, eop_data) where {To,Tf} rECItoECI(T, Val(To), JDo, Val(Tf), JDf, eop_data)

# State vector, ECEF to ECI

@deprecate svECEFtoECI(sv, ::Type{Val{To}}, ::Type{Val{Tf}}, JD)           where {To,Tf} svECEFtoECI(sv, Val(To), Val(Tf), JD)
@deprecate svECEFtoECI(sv, ::Type{Val{To}}, ::Type{Val{Tf}}, JD, eop_data) where {To,Tf} svECEFtoECI(sv, Val(To), Val(Tf), JD, eop_data)

# State vector, ECI to ECEF

@deprecate svECItoECEF(sv, ::Type{Val{To}}, ::Type{Val{Tf}}, JD)           where {To,Tf} svECItoECEF(sv, Val(To), Val(Tf), JD)
@deprecate svECItoECEF(sv, ::Type{Val{To}}, ::Type{Val{Tf}}, JD, eop_data) where {To,Tf} svECItoECEF(sv, Val(To), Val(Tf), JD, eop_data)
