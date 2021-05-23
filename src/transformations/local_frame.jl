# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#    Coordinate transformations related with local reference frames.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export ecef_to_ned, ned_to_ecef

"""
    ecef_to_ned(r_ecef::AbstractVector, lat::Number, lon::Number, h::Number; translate::Bool = false)

Convert a vector `r_ecef` represented in the Earth-Centered, Earth-Fixed (ECEF)
frame to the local reference frame NED (North, East, Down) at the geodetic
position `lat` [rad], `lon` [rad], and `h` [m].

If `translate` is `false`, then this function computes only the rotation between
ECEF and NED. Otherwise, it will also translate the vector considering the
distance between the Earth's center and NED origin.

# Remarks

This algorithm was based on the information in **[1]**.

# References

- **[1]**: [Transformations between ECEF and ENU
    coordinates](https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates)
"""
function ecef_to_ned(
    r_ecef::AbstractVector,
    lat::Number,
    lon::Number,
    h::Number;
    translate::Bool = false
)
    # Create the matrix that rotates the ECEF into NED.
    D_ned_ecef = angle_to_dcm(lon, -(lat + π / 2), 0, :ZYX)

    # Check if we need to translate the vector considering NED origins.
    if !translate
        Δr_ecef = r_ecef
    else
        # We need now to translate the vector. Thus, we need to obtain the ECEF
        # position of the NED origin.
        r_ned_ecef = geodetic_to_ecef(lat, lon, h)
        Δr_ecef = r_ecef - r_ned_ecef
    end

    # Now we can compute the vector in the NED.
    r_ned = D_ned_ecef * Δr_ecef

    return r_ned
end

"""
    ned_to_ecef(r_ned::AbstractVector, lat::Number, lon::Number, h::Number; translate::Bool = false)

Convert a vector `r_ned` represented in the local reference frame NED (North,
East, Down) at the geodetic position `lat` [rad], `lon` [rad], and `h` [m] to
the Earth-Centered, Earth-Fixed (ECEF) frame.

If `translate` is `false`, then this function computes only the rotation between
NED and ECEF. Otherwise, it will also translate the vector considering the
distance between the Earth's center and NED origin.

# Remarks

This algorithm was based on the information in **[1]**.

# References

- **[1]** [Transformations between ECEF and ENU
    coordinates](https://gssc.esa.int/navipedia/index.php/Transformations_between_ECEF_and_ENU_coordinates)
"""
function ned_to_ecef(
    r_ned::AbstractVector,
    lat::Number,
    lon::Number,
    h::Number;
    translate::Bool = false
)
    # Create the matrix that rotates the NED into ECEF.
    D_ecef_ned = angle_to_dcm(0, lat + π/2, -lon, :XYZ)

    # Now we can compute the vector in ECEF.
    r_ecef = D_ecef_ned * r_ned

    # Check if we need to translate the vector considering NED origins.
    if !translate
        return r_ecef
    else
        # We need now to translate the vector. Thus, we need to obtain the ECEF
        # position of the NED origin.
        r_ned_ecef = geodetic_to_ecef(lat, lon, h)
        return r_ecef + r_ned_ecef
    end
end
