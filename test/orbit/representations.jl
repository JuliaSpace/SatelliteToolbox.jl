# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to orbit representations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
# ==============================================================================
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/orbit/representations/rv.jl
# =======================================

# Functions: kepler_to_rv and rv_to_kepler
# ----------------------------------------

################################################################################
#                                 TEST RESULTS
################################################################################
#
# Scenario 01
# ===========
#
# Example 2-5: Finding position and velocity vectors (COE2RV Test Case) [1, p.
# 119-120].
#
# Cartesian representation:
#
#     r = 6525.344    I + 6861.535    J + 6449.125    K km
#     v =    4.902276 I +    5.533124 J -    1.975709 K km
#
# Orbit elements
#
#     ╔═════════════════╦══════════════╗
#     ║    Parameter    ║    Values    ║
#     ╠═════════════════╬══════════════╣
#     ║ p               ║ 11067.790 km ║
#     ║ Eccentricity    ║ 0.83285      ║
#     ║ Inclination     ║ 87.87°       ║
#     ║ RAAN            ║ 227.89°      ║
#     ║ Arg. of Perigee ║ 53.38°       ║
#     ║ True Anomaly    ║ 92.335°      ║
#     ╚═════════════════╩══════════════╝
#
################################################################################

@testset "Function kepler_to_rv and rv_to_kepler" begin

    # Float64
    # ==========================================================================

    let
        ## kepler_to_rv
        ## ============

        p    = 11067.790 * 1000
        e    = 0.83285
        i    = 87.87  * π / 180
        RAAN = 227.89 * π / 180
        w    = 53.38  * π / 180
        f    = 92.335 * π / 180
        a    = p / (1 - e^2)

        r_i, v_i = kepler_to_rv(a, e, i, RAAN, w, f)

        @test r_i[1] / 1000 ≈ +6525.344 atol=5e-2
        @test r_i[2] / 1000 ≈ +6861.535 atol=5e-2
        @test r_i[3] / 1000 ≈ +6449.125 atol=5e-2
        @test eltype(r_i) == Float64

        @test v_i[1] / 1000 ≈ +4.902276 atol=1e-4
        @test v_i[2] / 1000 ≈ +5.533124 atol=1e-4
        @test v_i[3] / 1000 ≈ -1.975709 atol=1e-4
        @test eltype(v_i) == Float64

        ## rv_to_kepler
        ## ============

        r_i = [6525.344; 6861.535; 6449.125] * 1000
        v_i = [4.902276; 5.533124; -1.975709] * 1000

        oe = rv_to_kepler(r_i..., v_i...)

        a, e, i, RAAN, w, f = oe.a, oe.e, oe.i, oe.Ω, oe.ω, oe.f
        p = a * (1 - e^2)

        @test p / 1000       ≈ 11067.790 atol=5e-2
        @test e              ≈ 0.83285   atol=1e-5
        @test i    * 180 / π ≈ 87.87     atol=1e-2
        @test RAAN * 180 / π ≈ 227.89    atol=1e-2
        @test w    * 180 / π ≈ 53.38     atol=1e-2
        @test f    * 180 / π ≈ 92.335    atol=1e-3
        @test oe.a isa Float64
    end

    # Float32
    # ==========================================================================

    let
        ## kepler_to_rv
        ## ============

        p    = 11067.790f0 * 1000
        e    = 0.83285f0
        i    = 87.87f0  * Float32(π / 180)
        RAAN = 227.89f0 * Float32(π / 180)
        w    = 53.38f0  * Float32(π / 180)
        f    = 92.335f0 * Float32(π / 180)
        a    = p / (1 - e^2)

        r_i, v_i = kepler_to_rv(a, e, i, RAAN, w, f)

        @test r_i[1] / 1000 ≈ +6525.344 atol=5e-2
        @test r_i[2] / 1000 ≈ +6861.535 atol=5e-2
        @test r_i[3] / 1000 ≈ +6449.125 atol=5e-2
        @test eltype(r_i) == Float32

        @test v_i[1] / 1000 ≈ +4.902276 atol=1e-4
        @test v_i[2] / 1000 ≈ +5.533124 atol=1e-4
        @test v_i[3] / 1000 ≈ -1.975709 atol=1e-4
        @test eltype(v_i) == Float32

        ## rv_to_kepler
        ## ============

        r_i = [6525.344f0; 6861.535f0; 6449.125f0]  * 1000
        v_i = [4.902276f0; 5.533124f0; -1.975709f0] * 1000

        oe = rv_to_kepler(r_i..., v_i...)

        a, e, i, RAAN, w, f = oe.a, oe.e, oe.i, oe.Ω, oe.ω, oe.f
        p = a * (1 - e^2)

        @test p / 1000       ≈ 11067.790 atol=5e-2
        @test e              ≈ 0.83285   atol=1e-5
        @test i * 180 / π    ≈ 87.87     atol=1e-2
        @test RAAN * 180 / π ≈ 227.89    atol=1e-2
        @test w * 180 / π    ≈ 53.38     atol=1e-2
        @test f * 180 / π    ≈ 92.335    atol=1e-3
        @test oe.a isa Float32
    end
end

@testset "Issue #25 - Special cases in kepler_to_rv" begin
    # Circular and equatorial
    # =======================

	orb_i = KeplerianElements(0.0, 42164e3, 0, 0, 0, 0, 0 )
    orb_o = sv_to_kepler( kepler_to_sv( orb_i ) )
	orb_e = KeplerianElements(0.0, 42164e3, 0, 0, 0, 0, 0 )

    @test orb_o.a ≈ orb_e.a   atol = 1e-7
    @test orb_o.e ≈ orb_e.e   atol = 1e-7
    @test orb_o.i ≈ orb_e.i   atol = 1e-7
    @test orb_o.Ω ≈ orb_e.Ω   atol = 1e-7
    @test orb_o.ω ≈ orb_e.ω   atol = 1e-7
    @test orb_o.f ≈ orb_e.f   atol = 1e-7

	orb_i = KeplerianElements(0.0, 42164e3, 0, 0, 10*π/180, 11*π/180, 12*π/180 )
    orb_o = sv_to_kepler( kepler_to_sv( orb_i ) )
    orb_e = KeplerianElements(0.0, 42164e3, 0, 0, 0, 0, (10+11+12)*π/180 )

    @test orb_o.a ≈ orb_e.a   atol = 1e-7
    @test orb_o.e ≈ orb_e.e   atol = 1e-7
    @test orb_o.i ≈ orb_e.i   atol = 1e-7
    @test orb_o.Ω ≈ orb_e.Ω   atol = 1e-7
    @test orb_o.ω ≈ orb_e.ω   atol = 1e-7
    @test orb_o.f ≈ orb_e.f   atol = 1e-7

    # Elliptical equatorial
    # =====================

    orb_i = KeplerianElements(0.0, 42164e3, 0.05, 0, 10*π/180, 11*π/180, 12*π/180 )
    orb_o = sv_to_kepler( kepler_to_sv( orb_i ) )
    orb_e = KeplerianElements(0.0, 42164e3, 0.05, 0, 0, (10+11)*π/180, 12*π/180 )

    @test orb_o.a ≈ orb_e.a   atol = 1e-7
    @test orb_o.e ≈ orb_e.e   atol = 1e-7
    @test orb_o.i ≈ orb_e.i   atol = 1e-7
    @test orb_o.Ω ≈ orb_e.Ω   atol = 1e-7
    @test orb_o.ω ≈ orb_e.ω   atol = 1e-7
    @test orb_o.f ≈ orb_e.f   atol = 1e-7

    # Circular inclined
    # =================

    orb_i = KeplerianElements(0.0, 42164e3, 0, π/4, 10*π/180, 11*π/180, 12*π/180 )
    orb_o = sv_to_kepler( kepler_to_sv( orb_i ) )
    orb_e = KeplerianElements(0.0, 42164e3, 0, π/4, 10*π/180, 0, (11+12)*π/180 )

    @test orb_o.a ≈ orb_e.a   atol = 1e-7
    @test orb_o.e ≈ orb_e.e   atol = 1e-7
    @test orb_o.i ≈ orb_e.i   atol = 1e-7
    @test orb_o.Ω ≈ orb_e.Ω   atol = 1e-7
    @test orb_o.ω ≈ orb_e.ω   atol = 1e-7
    @test orb_o.f ≈ orb_e.f   atol = 1e-7
end

@testset "Issue #72 - Bug with elliptical and equatorial orbit" begin
    a = 6674.790266053491 * 1000
    e = 0.0055622826070485095
    i = 0.0 |> deg2rad
    Ω = 0.0 |> deg2rad
    ω = 330.27258118831503 |> deg2rad
    f = 42.80749332919855 |> deg2rad

    rr, vv = kepler_to_rv(a, e, i, Ω, ω, f)

    oe = rv_to_kepler(rr, vv, 0)

    @test oe.a ≈ a atol = 1e-7
    @test oe.e ≈ e atol = 1e-7
    @test oe.i ≈ i atol = 1e-7
    @test oe.Ω ≈ Ω atol = 1e-7
    @test oe.ω ≈ ω atol = 1e-7
    @test oe.f ≈ f atol = 1e-7
end
