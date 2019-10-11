#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Dipole model for the Earth geomagnetic field.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] http://helios.fmi.fi/~juusolal/geomagnetism/Lectures/Chapter3_dipole.pdf
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export geomag_dipole

"""
    function geomag_dipole(r_e::AbstractVector, pole_lat::Number, pole_lon::Number, m::Number)

Compute the geomagnetic field [nT] using the simplified dipole model at position
`r_e` (ECEF reference frame). This function considers that the latitude of the
South magnetic pole (which lies in the North hemisphere) is `pole_lat` [rad] and
the longitude is `pole_lon` [rad]. Furthermore, the dipole moment is considered
to be `m` [A.m²].

    function geomag_dipole(r_e::AbstractVector, year::Number = 2019)

Compute the geomagnetic field [nT] using the simplified dipole model at position
`r_e` (ECEF reference frame). This function uses the year `year` to obtain the
position of the South magnetic pole (which lies in the North hemisphere) and the
dipole moment. If `year` is omitted, then it will be considered as 2019.

# Remarks

In both functions, the output vector will be represented in the ECEF reference
frame.

"""
function geomag_dipole(r_e::AbstractVector, year::Number = 2019)
    d2r = π/180

    # Obtain the pole location and the dipole moment.
    pole_lat = _itp_dipole_lat(year)*d2r
    pole_lon = _itp_dipole_lon(year)*d2r
    m        = _itp_dipole_mag(year)*1e22

    return geomag_dipole(r_e, pole_lat, pole_lon, m)
end

function geomag_dipole(r_e::AbstractVector, pole_lat::Number, pole_lon::Number,
                       m::Number)

    # DCM that convertes the ECEF into the geomagnetic coordinates.
    Dge = angle_to_dcm(pole_lon, π/2-pole_lat, 0, :ZYX)

    # Compute the dipole momentum represented in the ECEF reference frame.
    k₀_e = 1e-7m*(Dge'*SVector{3}(0,0,-1))

    # Compute the distance from the Earth center of the desired point.
    r = norm(r_e)

    # Compute the unitary vector that points to the desired direction.
    eᵣ_e = SVector{3}(r_e/r)

    # Compute the geomagnetic field vector [nT].
    B_e = 1/r^3*( 3eᵣ_e*eᵣ_e' - I )*k₀_e*1e9

    return B_e
end
