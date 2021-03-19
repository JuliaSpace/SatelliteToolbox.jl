# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to the nutation and EO computation of IAU-2006.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Wallace, P. T., Capitaine, N (2006). Precession-nutation procedures
#       consistent with IAU 2006 resolutions. Astronomy & Astrophysics.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/transformations/iau2006/nutation_eo.jl
# ==================================================

# Function nutation_eo_iau2006
# ----------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Numerical example at Appendix of [1]
#
# According to this example, using:
#
#   UTC    = 2006 January 15, 21h 24.375
#   JD_TT  = 2400000.5+53750.892855138888889
#   JD_UT1 = 2400000.5+53750.892104561342593
#
# one gets:
#
#   mϵ_2000 = +84378.576696215 arcsec (called ϵ_A in [1])
#   Δϵ_2000 =   +8.656841020   arcsec (called Δϵ in [1])
#   ΔΨ_2000 =   –1.071332969   arcsec (called ΔΨ in [1])
#   EO      = -277.646996035   arcsec
#
################################################################################

@testset "Function nutation_eo_iau2006" begin
    JD_TT  = 2400000.5+53750.892855138888889
    JD_UT1 = 2400000.5+53750.892104561342593

    r2a = 180/π*3600

    mϵ_2000, Δϵ_2000, ΔΨ_2000, EO = nutation_eo_iau2006(JD_TT)

    @test mϵ_2000*r2a ≈ 84378.576696215 atol = 1e-9
    @test Δϵ_2000*r2a ≈     8.656841020 atol = 1e-7
    @test ΔΨ_2000*r2a ≈    -1.071332969 atol = 1e-7
    @test EO*r2a      ≈  -277.646996035 atol = 1e-6
end
