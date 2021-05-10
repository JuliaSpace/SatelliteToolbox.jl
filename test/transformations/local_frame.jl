# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to the transformations of local reference frames.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/transformations/local_frame.jl
# ==========================================

# Functions: ecef_to_ned and ned_to_ecef
# --------------------------------------

@testset "Functions ecef_to_ned and ned_to_ecef" begin
    r_ecef = SVector(R0, 0, -R0)
    lat = 0.5
    lon = 0.44
    h = 130e3

    r_ned = ecef_to_ned(r_ecef, lat, lon, h; translate = true)

    @test r_ned[1] ≈  -8.345950969454717e6
    @test r_ned[2] ≈  -2.7167002618976594e6
    @test r_ned[3] ≈   4.496865561149651e6

    r_ecef_conv = ned_to_ecef(r_ned, lat, lon, h; translate = true)

    @test r_ecef_conv[1] ≈ r_ecef[1]
    @test r_ecef_conv[2] ≈ r_ecef[2]
    @test r_ecef_conv[3] ≈ r_ecef[3]

    r_ned = ecef_to_ned(r_ecef, lat, lon, h)
    @test norm(r_ned) == norm(r_ecef)

    B_l = SVector(24178.985570422887, -39579.98354881559, -14897.692106044478)

    # This is a validated code that rotates NED to ECEF without translations.
    D_ecef_l = angle_to_dcm(lat, -lon, 0, :YZX)
    B_ecef_exp = D_ecef_l*SVector(-B_l[3], +B_l[2], +B_l[1])
    B_ecef = ned_to_ecef(B_l, lat, lon, h)

    @test B_ecef[1] ≈ B_ecef_exp[1]
    @test B_ecef[2] ≈ B_ecef_exp[2]
    @test B_ecef[3] ≈ B_ecef_exp[3]
end

