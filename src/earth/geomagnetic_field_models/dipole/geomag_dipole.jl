# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Dipole model for the Earth geomagnetic field.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] http://helios.fmi.fi/~juusolal/geomagnetism/Lectures/Chapter3_dipole.pdf
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export geomag_dipole

"""
    geomag_dipole(r_e::AbstractVector{T}, year::Number = 2019) where T

Compute the geomagnetic field [nT] using the simplified dipole model at position
`r_e` (ECEF reference frame). This function considers that the latitude of the
South magnetic pole (which lies in the North hemisphere) is `pole_lat` [rad] and
the longitude is `pole_lon` [rad]. Furthermore, the dipole moment is considered
to be `m` [A.m²].

    geomag_dipole(r_e::AbstractVector, year::Number = 2019)

Compute the geomagnetic field [nT] using the simplified dipole model at position
`r_e` (ECEF reference frame). This function uses the year `year` to obtain the
position of the South magnetic pole (which lies in the North hemisphere) and the
dipole moment. If `year` is omitted, then it will be considered as 2019.

# Remarks

1. In both functions, the output vector will be represented in the ECEF
    reference frame.
2. The returned vector will have the same type `T` of the input vector.
"""
function geomag_dipole(r_e::AbstractVector{T}, year::Number = 2019) where T
    d2r = T(π / 180)

    # Obtain the pole location and the dipole moment.
    pole_lat = T(_itp_dipole_lat(year)) * d2r
    pole_lon = T(_itp_dipole_lon(year)) * d2r
    m        = T(_itp_dipole_mag(year)) * T(1e22)

    return geomag_dipole(r_e, pole_lat, pole_lon, m)
end

function geomag_dipole(
    r_e::AbstractVector{T},
    pole_lat::Number,
    pole_lon::Number,
    m::Number
) where T
    # DCM that converts the ECEF into the geomagnetic coordinates.
    Dge = angle_to_dcm(T(pole_lon), T(π / 2) - T(pole_lat), 0, :ZYX)

    # Compute the dipole momentum represented in the ECEF reference frame.
    k₀_e = T(1e-7) * m * (Dge' * SVector{3, T}(0, 0, -1))

    # Compute the distance from the Earth center of the desired point.
    r = norm(r_e)

    # Compute the unitary vector that points to the desired direction.
    eᵣ_e = SVector{3, T}(r_e) / r

    # Compute the geomagnetic field vector [nT].
    B_e = 1 / r^3 * (3 * eᵣ_e*eᵣ_e' - I) * k₀_e * T(1e9)

    return B_e
end
