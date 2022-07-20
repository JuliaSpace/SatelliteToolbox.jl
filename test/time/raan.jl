# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to the conversions between the local time at ascending or
#   descending node and RAAN.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/time/raan.jl
# ========================

# Function ltan_to_raan
# ---------------------

@testset "Function ltan_to_raan" begin
    jd₀ = date_to_jd(2021, 6, 19, 19, 35, 35)
    raan = ltan_to_raan(jd₀, 22.5)
    @test raan ≈ 4.289024646407403

    raan = ltan_to_raan(jd₀, Time("22:30:00"))
    @test raan ≈ 4.289024646407403
end

# Function ltdn_to_raan
# ---------------------

@testset "Function ltdn_to_raan" begin
    jd₀ = date_to_jd(2021, 6, 19, 19, 35, 35)
    raan = ltdn_to_raan(jd₀, 10.5)
    @test raan ≈ 4.289024646407403

    raan = ltdn_to_raan(jd₀, Time("10:30:00"))
    @test raan ≈ 4.289024646407403
end
