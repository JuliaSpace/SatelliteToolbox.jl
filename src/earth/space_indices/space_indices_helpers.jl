# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   This file contains helpers related to space indices.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

for sym in [:F10, :F10obs, :F10adj, :F10M, :F10Mobs, :F10Madj, :Kp, :Ap,
            :Kp_vect, :Ap_vect, :S10, :S81a, :M10, :M81a, :Y10, :Y81a, :DstΔTc]

    qsym = Meta.quot(sym)

    @eval begin
        export $sym
        $sym() = Val($qsym)
    end
end
