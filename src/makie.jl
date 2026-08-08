## Description #############################################################################
#
# Functions to create Makie themes for the SatelliteToolbox.jl ecosystem.
#
############################################################################################

export makie_palette, makie_theme

"""
    makie_palette(n::Int; kwargs...) -> Vector{Colorant}

Return the first `n` colors of the 6-color categorical palette used by the
SatelliteToolbox.jl Makie theme.

!!! note

    This function is defined in a package extension. Hence, Makie.jl must be loaded (e.g.,
    `using CairoMakie`) before calling it, otherwise it throws an `ErrorException`.

See also: [`makie_theme`](@ref)

# Keywords

- `dark::Bool`: If `true`, return the palette designed for dark backgrounds. Otherwise,
    return the palette designed for light backgrounds, matching the default theme returned
    by [`makie_theme`](@ref).
    (**Default**: `false`)

# Extended help

## Throws

- `ArgumentError`: If `n` is negative or greater than the number of colors in the palette.
- `ErrorException`: If Makie.jl is not loaded.
"""
function makie_palette(::Any; kwargs...)
    return error(
        "The function `makie_palette` is provided by a package extension. Load Makie.jl " *
        "(e.g., `using CairoMakie`) to use it. If Makie.jl is already loaded, check the " *
        "arguments: the valid call is `makie_palette(n::Int; kwargs...)`."
    )
end

"""
    makie_theme(; kwargs...) -> Makie.Theme
    makie_theme(variant::Symbol; kwargs...) -> Makie.Theme
    makie_theme(variant::Val; kwargs...) -> Makie.Theme

Return the SatelliteToolbox.jl Makie theme, ready to be applied with `set_theme!` or
`with_theme`.

`variant` selects the color scheme: `:dark` (or `Val(:dark)`) for dark backgrounds, and
`:light` (or `Val(:light)`) for light backgrounds. If `variant` is omitted, the light theme
is returned.

!!! note

    This function is defined in a package extension. Hence, Makie.jl must be loaded (e.g.,
    `using CairoMakie`) before calling it, otherwise it throws an `ErrorException`.

See also: [`makie_palette`](@ref)

# Keywords

- `fontscale::Real`: Factor to uniformly scale every font size (tick labels, axis labels,
    titles, legends, etc.). For example, `fontscale = 1.25` makes all text 25% larger,
    which is useful when a figure is shrunk into a small slide area.
    (**Default**: `1`)
- `mono_ticklabels::Bool`: If `true`, render the axis and colorbar tick labels using the
    monospaced font IBM Plex Mono. The fixed-width digits keep numeric ticks vertically
    aligned, which is useful for plots dominated by numbers.
    (**Default**: `false`)

# Extended help

The recommended figure size for 16:9 slides is:

    Figure(size = (1280, 720))

The recommended export resolution is:

    save("plot.png", fig; px_per_unit = 2)  # 2560 × 1440 px

## Throws

- `ArgumentError`: If the theme variant is not `:dark` or `:light`.
- `ErrorException`: If Makie.jl is not loaded.
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
