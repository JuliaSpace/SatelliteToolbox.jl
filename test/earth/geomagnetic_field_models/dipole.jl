#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests related to simplified dipole model of the Earth geomagnetic field.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] http://helios.fmi.fi/~juusolal/geomagnetism/Lectures/Chapter3_dipole.pdf
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

# File: ./src/earth/geomagnetic_field_models/dipole
# =================================================

# Function: geomag_dipole
# -----------------------

@testset "Function geomag_dipole" begin

    d2r = π/180

    # Test 1
    # --------------------------------------------------------------------------
    #
    # Position aligned with the dipole moment vector.

    year     = SatelliteToolbox._dipole_pole[32,1]
    pole_lat = SatelliteToolbox._dipole_pole[32,2] * d2r
    pole_lon = SatelliteToolbox._dipole_pole[32,3] * d2r
    m        = SatelliteToolbox._dipole_mag[32,2]  * 1e22

    # Distance from the Earth center.
    r = R0 + 196e3

    # Compute using the simplified equations [1] [nT].
    k0     = 1e-7m
    B_norm = 2k0 / r^3 * 1e9

    # Rotate to ECEF.
    Deg = angle_to_dcm(-(π / 2 - pole_lat), -pole_lon, 0, :YZX)
    B_g = SVector{3}(0, 0, -B_norm)
    B_e_expected = Deg * B_g

    # Compute using the function.
    r_g = SVector{3}(0, 0, r)
    r_e = Deg * r_g

    B_e_result_f64 = geomag_dipole(r_e, year)
    B_e_result_f32 = geomag_dipole(Float32.(r_e), year)

    @test B_e_expected ≈ B_e_result_f64 atol=1e-9
    @test eltype(B_e_result_f64) === Float64

    @test B_e_expected ≈ B_e_result_f32 atol=1e-1
    @test eltype(B_e_result_f32) === Float32

    # Test 2
    # --------------------------------------------------------------------------
    #
    # Position at the magnetic Equator, which must have half the magnitude of
    # that of the test 1.

    B_g = SVector{3}(0, 0, B_norm / 2)
    B_e_expected = Deg * B_g

    r_g = SVector{3}(r, 0, 0)
    r_e = Deg * r_g

    B_e_result_f64 = geomag_dipole(r_e, year)
    B_e_result_f32 = geomag_dipole(Float32.(r_e), year)

    @test B_e_expected ≈ B_e_result_f64 atol=1e-9
    @test eltype(B_e_result_f64) === Float64

    @test B_e_expected ≈ B_e_result_f32 atol=1e-1
    @test eltype(B_e_result_f32) === Float32
end
