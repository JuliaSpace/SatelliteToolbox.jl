# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to the EOP data.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/transformations/eop.jl
# ==================================

@testset "EOP data (IAU1980)" begin
    # Compare the fetched data to the one available locally
    # ==========================================================================

    # Fetch EOP from the internet.
    eop = get_iers_eop()

    # Parse the EOP from the local file.
    eop_local = read_iers_eop("./eop_IAU1980.txt")

    # The data in `eop_IAU1980.txt` is final. Hence, it must never change.
    JD = date_to_jd(2004, 4, 20, 19, 19, 19)
    @test eop.x(JD)           == eop_local.x(JD)
    @test eop.y(JD)           == eop_local.y(JD)
    @test eop.UT1_UTC(JD)     == eop_local.UT1_UTC(JD)
    @test eop.LOD(JD)         == eop_local.LOD(JD)
    @test eop.dPsi(JD)        == eop_local.dPsi(JD)
    @test eop.dEps(JD)        == eop_local.dEps(JD)
    @test eop.x_err(JD)       == eop_local.x_err(JD)
    @test eop.y_err(JD)       == eop_local.y_err(JD)
    @test eop.UT1_UTC_err(JD) == eop_local.UT1_UTC_err(JD)
    @test eop.LOD_err(JD)     == eop_local.LOD_err(JD)
    @test eop.dPsi_err(JD)    == eop_local.dPsi_err(JD)
    @test eop.dEps_err(JD)    == eop_local.dEps_err(JD)
    @test eop.x(JD)           == eop_local.x(JD)
    @test eop.y(JD)           == eop_local.y(JD)
    @test eop.UT1_UTC(JD)     == eop_local.UT1_UTC(JD)
    @test eop.LOD(JD)         == eop_local.LOD(JD)
    @test eop.dPsi(JD)        == eop_local.dPsi(JD)
    @test eop.dEps(JD)        == eop_local.dEps(JD)

    # Show
    # ==========================================================================

    buf = IOBuffer()
    io = IOContext(buf)
    show(io, MIME"text/plain"(), eop_local)
    result = String(take!(buf))

    expected = """
        EOPData_IAU1980:
             Data │ Timespan
         ─────────┼──────────────────────────────────────────────
                x │ 2004-04-01T00:00:00 -- 2004-04-30T00:00:00
                y │ 2004-04-01T00:00:00 -- 2004-04-30T00:00:00
          UT1-UTC │ 2004-04-01T00:00:00 -- 2004-04-30T00:00:00
              LOD │ 2004-04-01T00:00:00 -- 2004-04-30T00:00:00
               dψ │ 2004-04-01T00:00:00 -- 2004-04-30T00:00:00
               dϵ │ 2004-04-01T00:00:00 -- 2004-04-30T00:00:00"""

    @test result == expected

    buf = IOBuffer()
    io = IOContext(buf, :color => true)
    show(io, MIME"text/plain"(), eop_local)
    result = String(take!(buf))

    expected = """
        EOPData_IAU1980:
        \e[1m     Data \e[90m│ \e[0m\e[1mTimespan\e[0m
        \e[90m ─────────┼──────────────────────────────────────────────\e[0m
        \e[1m        x \e[90m│ \e[0m2004-04-01T00:00:00 -- 2004-04-30T00:00:00
        \e[1m        y \e[90m│ \e[0m2004-04-01T00:00:00 -- 2004-04-30T00:00:00
        \e[1m  UT1-UTC \e[90m│ \e[0m2004-04-01T00:00:00 -- 2004-04-30T00:00:00
        \e[1m      LOD \e[90m│ \e[0m2004-04-01T00:00:00 -- 2004-04-30T00:00:00
        \e[1m       dψ \e[90m│ \e[0m2004-04-01T00:00:00 -- 2004-04-30T00:00:00
        \e[1m       dϵ \e[90m│ \e[0m2004-04-01T00:00:00 -- 2004-04-30T00:00:00"""

    @test result == expected
end

@testset "EOP data (IAU2000A)" begin
    # Compare the fetched data to the one available locally
    # ==========================================================================

    # Fetch EOP from the internet.
    eop = get_iers_eop(Val(:IAU2000A))

    # Parse the EOP from the local file.
    eop_local = read_iers_eop("./eop_IAU2000A.txt", Val(:IAU2000A))

    # The data in `eop_IAU1980.txt` is final. Hence, it must never change.
    JD = date_to_jd(2004, 4, 20, 19, 19, 19)
    @test eop.x(JD)           == eop_local.x(JD)
    @test eop.y(JD)           == eop_local.y(JD)
    @test eop.UT1_UTC(JD)     == eop_local.UT1_UTC(JD)
    @test eop.LOD(JD)         == eop_local.LOD(JD)
    @test eop.dX(JD)          == eop_local.dX(JD)
    @test eop.dY(JD)          == eop_local.dY(JD)
    @test eop.x_err(JD)       == eop_local.x_err(JD)
    @test eop.y_err(JD)       == eop_local.y_err(JD)
    @test eop.UT1_UTC_err(JD) == eop_local.UT1_UTC_err(JD)
    @test eop.LOD_err(JD)     == eop_local.LOD_err(JD)
    @test eop.dX_err(JD)      == eop_local.dX_err(JD)
    @test eop.dY_err(JD)      == eop_local.dY_err(JD)
    @test eop.x(JD)           == eop_local.x(JD)
    @test eop.y(JD)           == eop_local.y(JD)
    @test eop.UT1_UTC(JD)     == eop_local.UT1_UTC(JD)
    @test eop.LOD(JD)         == eop_local.LOD(JD)
    @test eop.dX(JD)          == eop_local.dX(JD)
    @test eop.dY(JD)          == eop_local.dY(JD)

    # Show
    # ==========================================================================

    buf = IOBuffer()
    io = IOContext(buf)
    show(io, MIME"text/plain"(), eop_local)
    result = String(take!(buf))

    expected = """
        EOPData_IAU2000A:
             Data │ Timespan
         ─────────┼──────────────────────────────────────────────
                x │ 2004-04-01T00:00:00 -- 2004-04-30T00:00:00
                y │ 2004-04-01T00:00:00 -- 2004-04-30T00:00:00
          UT1-UTC │ 2004-04-01T00:00:00 -- 2004-04-30T00:00:00
              LOD │ 2004-04-01T00:00:00 -- 2004-04-30T00:00:00
               dX │ 2004-04-01T00:00:00 -- 2004-04-30T00:00:00
               dY │ 2004-04-01T00:00:00 -- 2004-04-30T00:00:00"""

    @test result == expected

    buf = IOBuffer()
    io = IOContext(buf, :color => true)
    show(io, MIME"text/plain"(), eop_local)
    result = String(take!(buf))

    expected = """
        EOPData_IAU2000A:
        \e[1m     Data \e[90m│ \e[0m\e[1mTimespan\e[0m
        \e[90m ─────────┼──────────────────────────────────────────────\e[0m
        \e[1m        x \e[90m│ \e[0m2004-04-01T00:00:00 -- 2004-04-30T00:00:00
        \e[1m        y \e[90m│ \e[0m2004-04-01T00:00:00 -- 2004-04-30T00:00:00
        \e[1m  UT1-UTC \e[90m│ \e[0m2004-04-01T00:00:00 -- 2004-04-30T00:00:00
        \e[1m      LOD \e[90m│ \e[0m2004-04-01T00:00:00 -- 2004-04-30T00:00:00
        \e[1m       dX \e[90m│ \e[0m2004-04-01T00:00:00 -- 2004-04-30T00:00:00
        \e[1m       dY \e[90m│ \e[0m2004-04-01T00:00:00 -- 2004-04-30T00:00:00"""

    @test result == expected
end
