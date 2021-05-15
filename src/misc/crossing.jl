# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Find crossing of a function.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

"""
    find_crossing(f::Function, t₀::Number, t₁::Number, s₀, s₁; Δ = 1e-3, max = 100)

Return the crossing time `tc` in which the function `f(t)` goes from the state
`s₀` to the state `s₁`. It is assumed that `f(t₀) = s₀` and `f(t₁) = s₁`.

If the computed interval is smalled than `Δ`, or if the number of iterations is
higher than `max`, then the algorithm stops.

# Examples

```julia-repl
julia> SatelliteToolbox.find_crossing(
    t -> (sin(t) > 0),
    -0.3,
    0.3,
    false,
    true;
    Δ = 1e-10
)
6.984919309616089e-11
```
"""
function find_crossing(
    f::Function,
    t₀::Number,
    t₁::Number,
    s₀,
    s₁,
    vargs...;
    Δ = 1e-3,
    max = 100
)
    it = 0

    T = typeof((t₁ + t₀) / 2)

    while it <= max
        # Call the function at the middle of the interval.
        ti = (t₁ + t₀) / 2
        si = f(ti, vargs...)

        # Compute the new interval.
        if si == s₀
            t₀ = ti
        elseif si == s₁
            t₁ = ti
        else
            error("The function `f` returned an unexpected state.")
        end

        # If the interval is small enough, then just return.
        (t₁ - t₀) < Δ && break

        it += 1
    end

    return T(t₁)
end
