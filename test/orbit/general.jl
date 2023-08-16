# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==========================================================================================
#
#   Tests related to general orbit functions.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# File: ./src/orbit/general.jl
# ==========================================================================================

# Function: orbital_angular_velocity
# ------------------------------------------------------------------------------------------

############################################################################################
#                                       Test Results
############################################################################################
#
# We will test the conversion by initializing the related propagator with the same inputs to
# the function and comparing the angular velocity.
#
############################################################################################

@testset "Function orbital_angular_velocity" begin
    for T in (Float64, Float32)
        a   = T(7130.984e3)
        e   = T(0.001111)
        i   = T(98.4106) |> deg2rad
        orb = KeplerianElements(0, a, e, i, 0, 0, 0)

        # J₀
        # ==================================================================================

        angvel = orbital_angular_velocity(a, e, i; perturbation = :J0)
        @test eltype(angvel) == T
        @test angvel ≈ √(GM_EARTH / a^3)

        angvel = orbital_angular_velocity(orb; perturbation = :J0)
        @test eltype(angvel) == T
        @test angvel ≈ √(GM_EARTH / a^3)

        # J₂
        # ==================================================================================

        orbp = Propagators.init(Val(:J2), orb)

        angvel = orbital_angular_velocity(
            a,
            e,
            i;
            J2 = 0.0010826261738522227,
            m0 = 3.986004415e14,
        )
        @test eltype(angvel) == T
        @test angvel ≈ orbp.j2d.n̄ + orbp.j2d.δω

        angvel = orbital_angular_velocity(orb)
        @test eltype(angvel) == T
        @test angvel ≈ orbp.j2d.n̄ + orbp.j2d.δω

        # J₄
        # ==================================================================================

        orbp = Propagators.init(Val(:J4), orb)

        angvel = orbital_angular_velocity(
            a,
            e,
            i;
            perturbation = :J4,
            J2 = 0.0010826261738522227,
            J4 = -1.6198975999169731e-6,
            m0 = 3.986004415e14,
        )
        @test eltype(angvel) == T
        @test angvel ≈ orbp.j4d.n̄ + orbp.j4d.δω

        angvel = orbital_angular_velocity(
            orb;
            perturbation = :J4,
            J2 = 0.0010826261738522227,
            J4 = -1.6198975999169731e-6,
            m0 = 3.986004415e14,
        )
        @test eltype(angvel) == T
        @test angvel ≈ orbp.j4d.n̄ + orbp.j4d.δω
    end
end

@testset "Function orbital_angular_velocity [ERRORS]" begin
    @test_throws ArgumentError orbital_angular_velocity(7130.9e3, 0, 0; perturbation = :J3)
end

# Function: orbital_angular_velocity_to_semimajor_axis
# ------------------------------------------------------------------------------------------

############################################################################################
#                                       Test Results
############################################################################################
#
# We will test the conversion by initializing the related propagator with the computed
# semi-major axis and comparing the angular velocity.
#
############################################################################################

@testset "Function orbital_angular_velocity_to_semimajor_axis" begin
    for T in (Float64, Float32)
        angvel = T(0.06) |> deg2rad
        e      = T(0.1)
        i      = T(98.4106) |> deg2rad

        # J₀
        # ==================================================================================

        â, conv = orbital_angular_velocity_to_semimajor_axis(angvel, e, i; perturbation = :J0)

        @test eltype(â) == T
        @test â ≈ (GM_EARTH / angvel^2)^(1 // 3)
        @test conv == true

        # J₂
        # ==================================================================================

        â, conv = orbital_angular_velocity_to_semimajor_axis(
            angvel,
            e,
            i;
            m0 = 3.986004415e14
        )
        orb  = KeplerianElements(0, â, e, i, 0, 0, 0)
        orbp = Propagators.init(Val(:J2), orb)

        @test eltype(â) == T
        @test (orbp.j2d.n̄ + orbp.j2d.δω) ≈ angvel
        @test conv == true

        # J₄
        # ==================================================================================

        â, conv = orbital_angular_velocity_to_semimajor_axis(
            angvel,
            e,
            i;
            perturbation = :J4,
            m0 = 3.986004415e14
        )
        orb  = KeplerianElements(0, â, e, i, 0, 0, 0)
        orbp = Propagators.init(Val(:J4), orb)

        @test eltype(â) == T
        @test (orbp.j4d.n̄ + orbp.j4d.δω) ≈ angvel
        @test conv == true
    end
end

@testset "Function orbital_angular_velocity_to_semimajor_axis [ERRORS]" begin
    @test_throws ArgumentError orbital_angular_velocity_to_semimajor_axis(0.001, 0, 0; perturbation = :J3)
end

# Function: orbital_period
# ------------------------------------------------------------------------------------------

############################################################################################
#                                       Test Results
############################################################################################
#
# We will test the conversion by initializing the related propagator with the same inputs to
# the function and comparing the period.
#
############################################################################################

@testset "Function orbital_period" begin
    for T in (Float64, Float32)
        a   = T(7130.984e3)
        e   = T(0.001111)
        i   = T(98.4106) |> deg2rad
        orb = KeplerianElements(0, a, e, i, 0, 0, 0)

        # J₀
        # ==================================================================================

        P = orbital_period(a, e, i; perturbation = :J0)
        @test eltype(P) == T
        @test P ≈ 2π / √(GM_EARTH / a^3)

        P = orbital_period(orb; perturbation = :J0)
        @test eltype(P) == T
        @test P ≈ 2π / √(GM_EARTH / a^3)

        # J₂
        # ==================================================================================

        orbp = Propagators.init(Val(:J2), orb)

        P = orbital_period(a, e, i; m0 = 3.986004415e14, J2 = 0.0010826261738522227)
        @test eltype(P) == T
        @test P ≈ T(2π) / (orbp.j2d.n̄ + orbp.j2d.δω)

        P = orbital_period(orb)
        @test eltype(P) == T
        @test P ≈ T(2π) / (orbp.j2d.n̄ + orbp.j2d.δω)

        # J₄
        # ==================================================================================

        orbp = Propagators.init(Val(:J4), orb)

        P = orbital_period(
            a,
            e,
            i;
            perturbation = :J4,
            J2 = 0.0010826261738522227,
            J4 = -1.6198975999169731e-6,
            m0 = 3.986004415e14,
        )
        @test eltype(P) == T
        @test P ≈ T(2π) / (orbp.j4d.n̄ + orbp.j4d.δω)

        P = orbital_period(
            orb;
            perturbation = :J4,
            J2 = 0.0010826261738522227,
            J4 = -1.6198975999169731e-6,
            m0 = 3.986004415e14,
        )
        @test eltype(P) == T
        @test P ≈ T(2π) / (orbp.j4d.n̄ + orbp.j4d.δω)
    end
end

# Function: raan_time_derivative
# ------------------------------------------------------------------------------------------

############################################################################################
#                                       Test Results
############################################################################################
#
# We will test the conversion by initializing the related propagator with the same inputs to
# the function and comparing the RAAN time-derivative.
#
############################################################################################

@testset "Function raan_time_derivative" begin
    for T in (Float64, Float32)
        a   = T(7130.984e3)
        e   = T(0.001111)
        i   = T(98.4106) |> deg2rad
        orb = KeplerianElements(0, a, e, i, 0, 0, 0)

        # J₀
        # ==================================================================================

        ∂Ω = raan_time_derivative(a, e, i; perturbation = :J0)
        @test eltype(∂Ω) == T
        @test ∂Ω == T(0)

        ∂Ω = raan_time_derivative(orb; perturbation = :J0)
        @test eltype(∂Ω) == T
        @test ∂Ω == T(0)

        # J₂
        # ==================================================================================

        orbp = Propagators.init(Val(:J2), orb)

        ∂Ω = raan_time_derivative(a, e, i; m0 = 3.986004415e14, J2 = 0.0010826261738522227)
        @test eltype(∂Ω) == T
        @test ∂Ω ≈ orbp.j2d.δΩ

        ∂Ω = raan_time_derivative(orb)
        @test eltype(∂Ω) == T
        @test ∂Ω ≈ orbp.j2d.δΩ

        # J₄
        # ==================================================================================

        orbp = Propagators.init(Val(:J4), orb)

        ∂Ω = raan_time_derivative(
            a,
            e,
            i;
            perturbation = :J4,
            J2 = 0.0010826261738522227,
            J4 = -1.6198975999169731e-6,
            m0 = 3.986004415e14,
        )
        @test eltype(∂Ω) == T
        @test ∂Ω ≈ orbp.j4d.δΩ

        ∂Ω = raan_time_derivative(
            orb;
            perturbation = :J4,
            J2 = 0.0010826261738522227,
            J4 = -1.6198975999169731e-6,
            m0 = 3.986004415e14,
        )
        @test eltype(∂Ω) == T
        @test ∂Ω ≈ orbp.j4d.δΩ
    end
end

@testset "Function raan_time_derivative [ERRORS]" begin
    @test_throws ArgumentError raan_time_derivative(7000e3, 0, 0; perturbation = :J3)
end
