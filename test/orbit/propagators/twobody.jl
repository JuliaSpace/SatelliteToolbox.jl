#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Tests of the Two Body orbit propagator.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-04-08: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

################################################################################
#                                 Test Results
################################################################################
#
# Scenario 01
# ===========
#
# Example 2-4. Solving Kepler's Problem [1, p. 94-95].
#
#   Initial position:
#
#       r0_ijk = + 1131.340 i - 2282.343 j + 6672.423 k [km]
#       v0_ijk = - 5.64305  i + 4.30333  j + 2.42879  k [km/s]
#
#    After 40 min., the solution of Kepler's problem leads to the following
#    position:
#
#       rf_ijk = - 4219.7527 i + 4363.0292 j - 3958.7666 k [km]
#       vf_ijk = - 3.689866  i - 1.916735  j - 6.112511  k [km/s]
#
#   Using the SatelliteToolbox Two Body orbit propagator, the following result
#   was obtained:
#
#       julia> orb = rv_to_kepler(1131340, -2282343, 6672423, -5643.05, 4303.33, 2428.79)
#       julia> orbp = init_orbit_propagator(Val{:twobody}, orb)
#       julia> (o,r,v) = step!(orbp,40*60)
#       julia> @printf("rf_ijk = %+f i %+f j %+f k", r[1]/1000, r[2]/1000, r[3]/1000)
#           rf_ijk = -4219.752738 i +4363.029177 j -3958.766617 k
#       julia> @printf("vf_ijk = %+f i %+f j %+f k", v[1]/1000, v[2]/1000, v[3]/1000)
#           vf_ijk = +3.689866 i -1.916735 j -6.112511 k
#
#   Where it is possible to see that the solution is exactly the same as
#   presented in [1].
#
################################################################################

@testset "Two Body orbit propagator" begin
    orb = rv_to_kepler(1131340., -2282343., 6672423., -5643.05, 4303.33, 2428.79)
    orbp = init_orbit_propagator(Val{:twobody}, orb)
    (o,r,v) = step!(orbp,40*60)

    # Testing position.
    @test r[1]/1000 ≈ -4219.7527 atol=1e-3
    @test r[2]/1000 ≈ +4363.0292 atol=1e-3
    @test r[3]/1000 ≈ -3958.7666 atol=1e-3

    # Testing velocity.
    @test v[1]/1000 ≈ +3.689866 atol=1e-6
    @test v[2]/1000 ≈ -1.916735 atol=1e-6
    @test v[3]/1000 ≈ -6.112511 atol=1e-6
end
