# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   The exponential atmospheric model.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export expatmosphere

"""
    function expatmosphere(h::Number)

Compute the atmospheric density [kg/m³] at the altitude `h` \\[m] (above the
ellipsoid) using the exponential atmospheric model:

                    ┌            ┐
                    │    h - h₀  │
    ρ(h) = ρ₀ ⋅ exp │ - ──────── │ ,
                    │      H     │
                    └            ┘

in which `ρ₀`, `h₀`, and `H` are parameters obtained from tables that depend
only on `h`.

"""
function expatmosphere(h::Number)
    # Check the bounds.
    (h < 0) && @error "h must be positive."

    # Transform `h` to km.
    h /= 1000

    # Get the values for the exponential model.
    id = (h >= 1000) ? 28 : findfirst( (x)->x > 0 , _expatmosphere_h₀ .- h ) - 1

    h₀ = _expatmosphere_h₀[id]
    ρ₀ = _expatmosphere_ρ₀[id]
    H  = _expatmosphere_H[id]

    # Compute the density.
    ρ₀ * exp( - ( h - h₀ )/H )
end

################################################################################
#                              Private Variables
################################################################################

"""
Base altitude for the exponential atmospheric model [km].

"""
const _expatmosphere_h₀ = [0, 25, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120,
                           130, 140, 150, 180, 200, 250, 300, 350, 400, 450,
                           500, 600, 700, 800, 900, 1000]

"""
Nominal density for the exponential atmospheric model [kg/m³].

"""
const _expatmosphere_ρ₀ = [1.225, 3.899e-2, 1.774e-2, 3.972e-3, 1.057e-3,
                           3.206e-4, 8.770e-5, 1.905e-5, 3.396e-6, 5.297e-7,
                           9.661e-8, 2.438e-8, 8.484e-9, 3.845e-9, 2.070e-9,
                           5.464e-10, 2.789e-10, 7.248e-11, 2.418e-11,
                           9.518e-12, 3.725e-12, 1.585e-12, 6.967e-13,
                           1.454e-13, 3.614e-14, 1.170e-14, 5.245e-15,
                           3.019e-15]

"""
Scale height for the exponential atmospheric model [km].

"""
const _expatmosphere_H = [7.249, 6.349, 6.682, 7.554, 8.382, 7.714, 6.549,
                          5.799, 5.382, 5.877, 7.263, 9.473, 12.636, 16.149,
                          22.523, 29.740, 37.105, 45.546, 53.628, 53.298,
                          58.515, 60.828, 63.822, 71.835, 88.667, 124.64,
                          181.05, 268.00]
