# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Tests related to the transformations between Geodetic and Geocentric frames.
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

# File: ./src/transformations/geodetic_geocentric.jl
# ==================================================

# Function: ecef_to_geodetic
# --------------------------

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 3-3: Converting ECEF to Lat Lon [1, p. 173].
#
# According to this example, using:
#
#   r_ecef = 6524.834 i + 6862.875 j + 6448.296 k [km]
#
# one gets:
#
#   Geodetic Latitude  = 34.352496°
#            Longitude = 46.4464°
#            Altitude  = 5085.22 km
#
# Scenario 02
# ===========
#
# At the poles, we have a singularity. We know that, if:
#
#   r_ecef = 0 i + 0 j + Z k
#
# then
#
#   Geodetic Latitude  = 90° for Z > 0 and -90° for Z < 0
#            Longitude =  0°
#            Altitude  = Z - b_wgs84
#
################################################################################

@testset "Function ecef_to_geodetic" begin
    # Scenario 01
    # ===========

    r = [6524.834e3, 6862.875e3, 6448.296e3]

    ϕ_gd, λ_gd, h = ecef_to_geodetic(r)

    @test rad2deg(ϕ_gd) ≈ 34.352496 atol = 1e-6
    @test rad2deg(λ_gd) ≈ 46.4464   atol = 1e-4
    @test h/1000        ≈ 5085.22   atol = 1e-2

    # Scenario 02
    # ===========

    aux = rand(0:1000)

    Z = R0 + aux
    ϕ_gd, λ_gd, h = ecef_to_geodetic([0;0;Z])

    @test rad2deg(ϕ_gd) ≈ 90
    @test rad2deg(λ_gd) ≈ 0
    @test h             ≈ Z - SatelliteToolbox.b_wgs84

    Z = -R0 + aux
    ϕ_gd, λ_gd, h = ecef_to_geodetic([0;0;Z])

    @test rad2deg(ϕ_gd) ≈ -90
    @test rad2deg(λ_gd) ≈ 0
    @test h             ≈ -Z - SatelliteToolbox.b_wgs84
end

################################################################################
#                                 Test Results
################################################################################
#
# Using the FORTRAN code available in:
#
#   https://www.astro.uni.torun.pl/~kb/Papers/geod/Geod-BG.htm
#
# we obtained the following matrix of data:
#
# | Geocentric lat [rad] | Distance [m] | Geodetic lat [rad]   | Altitude [m]        |
# |----------------------|--------------|----------------------|---------------------|
# |  19.86               | R0 + 1987    |  1.0134554245512695  | 17352.756962650223  |
# |  π/2                 | R0 + 190686  |  1.5707963267948966  | 212070.68575482070  |
# |  31.35 * π/180       | R0 + 752000  |  0.54808101129276687 | 757773.37237201189  |
# | -19.86               | R0 + 1987    | -1.0134554245512695  | 17352.756962650223  |
# | -π/2                 | R0 + 190686  | -1.5707963267948966  | 212070.68575482070  |
# | -31.35 * π/180       | R0 + 752000  | -0.54808101129276687 | 757773.37237201189  |
# |  0                   | R0           |  0.0000000000000000  | 0.0000000000000000  |
# |  15                  | 10000        |  1.3567978765139961  | -6353137.8454352869 |
# | -15                  | 10000        | -1.3567978765139961  | -6353137.8454352869 |
#
################################################################################
#
# FORTRAN code used to obtain the results
# ==============================================================================
#
#      subroutine GEOD(r,z,fi,h)
#          implicit real*8(a-h,o-z)
#          data a,frec /6378137.d0,298.257223563d0/
#          b=dsign(a-a/frec,z)
#          if(r.eq.0d0) return
#          E=((z+b)*b/a-a)/r
#          F=((z-b)*b/a+a)/r
#          P=(E*F+1.)*4d0/3.
#          Q=(E*E-F*F)*2.
#          D=P*P*P+Q*Q
#          if(D.ge.0d0) then
#              s=dsqrt(D)+Q
#              s=dsign(dexp(dlog(dabs(s))/3d0),s)
#              v=P/s-s
#              v=-(Q+Q+v*v*v)/(3*P)
#          else
#              v=2.*dsqrt(-P)*dcos(dacos(Q/P/dsqrt(-P))/3.)
#          endif
#          G=.5*(E+dsqrt(E*E+v))
#          t=dsqrt(G*G+(F-v*G)/(G+G-E))-G
#          fi=datan((1.-t*t)*a/(2*b*t))
#          h=(r-a*t)*dcos(fi)+(z-b)*dsin(fi)
#          end
#
#      program TEST
#          implicit real*8(a-h,o-z)
#          R0 = 6378137.d0
#          PI = 4 * atan(1.d0)
#
#          dlat = 19.86d0
#          r    = R0 + 1987.d0
#          call GEOD(r * dcos(dlat), r * dsin(dlat), fi, h)
#          PRINT *, fi, h
#
#          dlat = 90.00d0 * PI / 180.d0
#          r    = R0 + 190686.d0
#          call GEOD(r * dcos(dlat), r * dsin(dlat), fi, h)
#          PRINT *, fi, h
#
#          dlat = 31.25d0 * PI / 180.d0
#          r    = R0 + 752000.d0
#          call GEOD(r * dcos(dlat), r * dsin(dlat), fi, h)
#          PRINT *, fi, h
#
#          dlat = -19.86d0
#          r    = R0 + 1987.d0
#          call GEOD(r * dcos(dlat), r * dsin(dlat), fi, h)
#          PRINT *, fi, h
#
#          dlat = -90.00d0 * PI / 180.d0
#          r    = R0 + 190686.d0
#          call GEOD(r * dcos(dlat), r * dsin(dlat), fi, h)
#          PRINT *, fi, h
#
#          dlat = -31.25d0 * PI / 180.d0
#          r    = R0 + 752000.d0
#          call GEOD(r * dcos(dlat), r * dsin(dlat), fi, h)
#          PRINT *, fi, h
#
#          dlat = 0.d0
#          r    = R0
#          call GEOD(r * dcos(dlat), r * dsin(dlat), fi, h)
#          PRINT *, fi, h
#
#          dlat = 15.d0 * PI / 180.d0
#          r    = 10000.d0
#          call GEOD(r * dcos(dlat), r * dsin(dlat), fi, h)
#          PRINT *, fi, h
#
#          dlat = -15.d0 * PI / 180.d0
#          r    = 10000.d0
#          call GEOD(r * dcos(dlat), r * dsin(dlat), fi, h)
#          PRINT *, fi, h
#      end program
#
################################################################################


# Function: geocentric_to_geodetic
# --------------------------------

@testset "Function geocentric_to_geodetic" begin
    # North hemisphere
    # ==========================================================================

    ϕ_gd, h = geocentric_to_geodetic(19.86, R0 + 1987)
    @test ϕ_gd ≈ 1.0134554245512695
    @test h    ≈ 17352.756962650223

    ϕ_gd, h = geocentric_to_geodetic(π/2, R0 + 190686)
    @test ϕ_gd ≈ 1.5707963267948966
    @test h    ≈ 212070.68575482070

    ϕ_gd, h = geocentric_to_geodetic(deg2rad(31.25), R0 + 752000)
    @test ϕ_gd ≈ 0.54808101129276687
    @test h    ≈ 757773.37237201189

    # South hemisphere
    # ==========================================================================

    ϕ_gd, h = geocentric_to_geodetic(-19.86, R0 + 1987)
    @test ϕ_gd ≈ -1.0134554245512695
    @test h    ≈ 17352.756962650223

    ϕ_gd, h = geocentric_to_geodetic(-π/2, R0 + 190686)
    @test ϕ_gd ≈ -1.5707963267948966
    @test h    ≈ 212070.68575482070

    ϕ_gd, h = geocentric_to_geodetic(-deg2rad(31.25), R0 + 752000)
    @test ϕ_gd ≈ -0.54808101129276687
    @test h    ≈ 757773.37237201189

    # Special cases
    # ==========================================================================

    ϕ_gd, h = geocentric_to_geodetic(0, R0)
    @test ϕ_gd ≈ 0.0
    @test h    ≈ 0.0

    # D < 0
    ϕ_gd, h = geocentric_to_geodetic(deg2rad(15), 10e3)
    @test ϕ_gd ≈ 1.3567978765139961
    @test h    ≈ -6353137.8454352869

    ϕ_gd, h = geocentric_to_geodetic(-deg2rad(15), 10e3)
    @test ϕ_gd ≈ -1.3567978765139961
    @test h    ≈ -6353137.8454352869
end
