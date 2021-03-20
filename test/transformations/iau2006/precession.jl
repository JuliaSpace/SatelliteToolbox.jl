# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to the precession computation of IAU-2006.
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

# File: ./src/transformations/iau2006/precession.jl
# =================================================

# Function precession_iau2006
# ---------------------------

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
#   Ψ_a = +304.359364139   arcsec
#   ω_a = +84381.404629617 arcsec
#   χ_a = +0.628998164     arcsec
#
################################################################################

@testset "Function precession_iau2006" begin
    JD_TT  = 2400000.5+53750.892855138888889

    r2a = 180/π*3600

    Ψ_a, ω_a, χ_a = precession_iau2006(JD_TT)

    @test Ψ_a*r2a ≈   +304.359364139 atol = 1e-9
    @test ω_a*r2a ≈ +84381.404629617 atol = 1e-9
    @test χ_a*r2a ≈     +0.628998164 atol = 1e-9
end
