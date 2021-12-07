# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Miscellaneous functions for the functions of IAU-2006 theory.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Compute the sum as used in many computations of IAU-2006 theory. The first
# parameter is a tuple with the coefficient matrices, the second is the Julian
# century [TT], and the rest are the fundamental arguments.
@inline function _iau2006_sum(
    coefs::Tuple,
    t_tt::Number,
    M_s::Number,
    M_m::Number,
    u_Mm::Number,
    D_s::Number,
    Ω_m::Number,
    λ_M☿::Number,
    λ_M♀::Number,
    λ_Me::Number,
    λ_M♂::Number,
    λ_M♃::Number,
    λ_M♄::Number,
    λ_M⛢::Number,
    λ_M♆::Number,
    p_λ::Number
)
    # Result of the sum.
    r = 0.0

    # Auxiliary variable to compute the powers of t_tt.
    t_tt_power = one(t_tt)

    # Number of sums.
    num_sums = length(coefs)

    @inbounds for i in 1:num_sums
        ci = coefs[i]

        # Result of this sum.
        rp = zero(r)

        # Number of coefficients in this sum.
        num_coefs = size(ci, 2)

        # Notice that the matrices were transposed when created to improve the
        # performance due to memory alignment.
        for j = 1:num_coefs
            As     = ci[ 2,j]
            Ac     = ci[ 3,j]
            ap     = ci[ 4,j] * M_m  + ci[ 5,j] * M_s  + ci[ 6,j] * u_Mm +
                     ci[ 7,j] * D_s  + ci[ 8,j] * Ω_m  + ci[ 9,j] * λ_M☿ +
                     ci[10,j] * λ_M♀ + ci[11,j] * λ_Me + ci[12,j] * λ_M♂ +
                     ci[13,j] * λ_M♃ + ci[14,j] * λ_M♄ + ci[15,j] * λ_M⛢ +
                     ci[16,j] * λ_M♆ + ci[17,j] * p_λ
            sj, cj = sincos(ap)
            rp    += (As * sj + Ac * cj) / 1e6
        end

        # Accumulate in the output variable.
        r += rp * t_tt_power

        # Update the t_tt power for the next pass.
        t_tt_power *= t_tt
    end

    return r
end
