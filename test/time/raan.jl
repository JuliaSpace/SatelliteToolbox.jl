# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Tests related to the conversions between the local time at ascending or
#   descending node and RAAN.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/time/raan.jl
# ==========================================================================================

# Function ltan_to_raan
# ------------------------------------------------------------------------------------------

@testset "Function ltan_to_raan" begin
    t₀ = DateTime("2021-06-19T19:35:35")
    raan = ltan_to_raan(22.5, t₀)
    @test raan ≈ 4.289024646407403

    raan = ltan_to_raan(Time("22:30:00"), t₀)
    @test raan ≈ 4.289024646407403
end

# Function ltdn_to_raan
# ------------------------------------------------------------------------------------------

@testset "Function ltdn_to_raan" begin
    t₀ = DateTime("2021-06-19T19:35:35")
    raan = ltdn_to_raan(10.5, t₀)
    @test raan ≈ 4.289024646407403

    raan = ltdn_to_raan(Time("10:30:00"), t₀)
    @test raan ≈ 4.289024646407403
end

# Function raan_to_ltan
# ------------------------------------------------------------------------------------------

@testset "Function raan_to_ltan" begin
    t₀ = DateTime("2021-06-19T19:35:35")
    ltan = raan_to_ltan(4.289024646407403, t₀)
    @test ltan ≈ 22.5
end

# Function raan_to_ltdn
# ------------------------------------------------------------------------------------------

@testset "Function raan_to_ltdn" begin
    t₀ = DateTime("2021-06-19T19:35:35")
    ltdn = raan_to_ltdn(4.289024646407403, t₀)
    @test ltdn ≈ 10.5
end
