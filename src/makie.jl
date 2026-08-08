## Description #############################################################################
#
# Functions for creating Makie theme for the SatelliteToolbox.jl ecosystem.
#
############################################################################################

"""
    makie_palette(n::Int; dark::Bool = true) → Vector{Colorant}

Return the first `n` colors from the categorical palette. Set `dark = false` to get the
light-background variant. Throws `ArgumentError` if `n > 6`.
"""
function makie_palette(::Any; kwargs...)
    return error(
        "The function `makie_palette` is provided by a package extension. Load Makie.jl " *
        "(e.g., `using CairoMakie`) to use it. If Makie.jl is already loaded, check the " *
        "arguments: the valid call is `makie_palette(n::Int; kwargs...)`."
    )
end

"""
    makie_theme(; kwargs...) -> Theme
    makie_theme(T; kwargs...) -> Theme

Return the SatelliteToolbox.jl Makie theme, ready for `set_theme!` or `with_theme`. For the
light version, use `T = Val(:light)`, for the dark version, use `T = Val(:dark)`.

Recommended figure size for 16:9 slides:

    Figure(size = (1280, 720))

Recommended export resolution:

    save("plot.png", fig, px_per_unit = 2)  # 2560 × 1440 px

# Keywords

- `fontscale::Real`: Factor to uniformly scale every font size (tick labels, axis labels,
    titles, legends, etc.) by that factor. For example, `fontscale = 1.25` makes all text
    25% larger, which is useful when a figure is shrunk into a small slide area.
    (**Default** = `1`)
- `mono_ticklabels::Bool`: If `true`, render the axis and colorbar tick labels in IBM Plex
    Mono.
    (**Default** = `false`)
"""
function makie_theme(; kwargs...)
    return makie_theme(Val(:light); kwargs...)
end

function makie_theme(variant::Symbol; kwargs...)
    return makie_theme(Val(variant); kwargs...)
end

function makie_theme(::Any; kwargs...)
    return error(
        "The function `makie_theme` is provided by a package extension. Load Makie.jl " *
        "(e.g., `using CairoMakie`) to use it. If Makie.jl is already loaded, check the " *
        "arguments: the valid calls are `makie_theme()`, `makie_theme(:dark)`, and " *
        "`makie_theme(:light)`."
    )
end
