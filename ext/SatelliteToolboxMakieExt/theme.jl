## Description #############################################################################
#
# Functions to create Makie theme for the SatelliteToolbox ecosystem.
#
############################################################################################

############################################################################################
#                                     Public Functions                                     #
############################################################################################

function makie_palette(n::Int; dark::Bool = false)
    colors = dark ? CATEGORICAL_DARK : CATEGORICAL_LIGHT

    (0 <= n <= length(colors)) || throw(
        ArgumentError("The categorical palette has $(length(colors)) colors; requested: $n.")
    )

    return colors[1:n]
end

function makie_theme(::Val{:dark}; fontscale::Real = 1, mono_ticklabels::Bool = false)
    return _build_theme(;
        bg              = NAVY_PRIMARY,
        card            = NAVY_CARD,
        separator       = SEPARATOR_DARK,
        border          = BORDER_DARK,
        text_primary    = TEXT_PRIMARY_DARK,
        text_secondary  = TEXT_SECONDARY_DARK,
        text_tertiary   = TEXT_TERTIARY_DARK,
        amber           = AMBER_DARK,
        cyan            = CYAN_DARK,
        magenta         = MAGENTA_DARK,
        categorical     = CATEGORICAL_DARK,
        fontscale       = fontscale,
        mono_ticklabels = mono_ticklabels,
    )
end

function makie_theme(::Val{:light}; fontscale::Real = 1, mono_ticklabels::Bool = false)
    return _build_theme(;
        bg              = SURFACE,
        card            = SURFACE_CARD,
        separator       = SEPARATOR_LIGHT,
        border          = BORDER_LIGHT,
        text_primary    = TEXT_PRIMARY_LIGHT,
        text_secondary  = TEXT_SECONDARY_LIGHT,
        text_tertiary   = TEXT_TERTIARY_LIGHT,
        amber           = AMBER_LIGHT,
        cyan            = CYAN_LIGHT,
        magenta         = MAGENTA_LIGHT,
        categorical     = CATEGORICAL_LIGHT,
        fontscale       = fontscale,
        mono_ticklabels = mono_ticklabels,
    )
end

function makie_theme(::Val{T}; kwargs...) where T
    throw(ArgumentError(
        "Unknown Makie theme variant `:$T`. The available options are `:dark` and `:light`."
    ))
end

############################################################################################
#                                    Private Functions                                     #
############################################################################################

# Directory holding the bundled `.ttf` font files, relative to this source file.
const _FONT_DIR = normpath(joinpath(@__DIR__, "..", "..", "assets", "fonts"))

# Return the absolute path to a bundled font file, so Makie loads it directly instead of
# looking the family up in the host system's installed fonts.
_font(file::AbstractString) = joinpath(_FONT_DIR, file)

"""
    _build_theme(; kwargs...) -> Theme

Assemble a Makie `Theme` from the provided color roles and font scale. Both `dark_theme`
and `light_theme` call this with their respective palettes, so all styling common to the
two variants lives here.

# Keywords

- `bg::Colorant`: Figure background color.
- `card::Colorant`: Slightly offset surface color for the axis background, legends, and
    tooltips.
- `separator::Colorant`: Color of the major and minor grid lines.
- `border::Colorant`: Color of box outlines such as the legend frame.
- `text_primary::Colorant`: Primary text color for body text, axis labels, and titles.
- `text_secondary::Colorant`: Secondary text color for tick labels, spines, and metadata.
- `text_tertiary::Colorant`: Tertiary text color for minor ticks.
- `amber::Colorant`: Primary accent color.
- `cyan::Colorant`: Secondary accent color used for legend and colorbar titles.
- `magenta::Colorant`: Tertiary accent color.
- `categorical::AbstractVector{<:Colorant}`: Categorical color cycle for data series.
- `fontscale::Real`: Factor by which every font size is uniformly multiplied.
    (**Default**: `1`)
- `mono_ticklabels::Bool`: If `true`, render the axis and colorbar tick labels in IBM Plex
    Mono instead of the regular body font.
    (**Default**: `false`)
"""
function _build_theme(;
    bg::Colorant,
    card::Colorant,
    separator::Colorant,
    border::Colorant,
    text_primary::Colorant,
    text_secondary::Colorant,
    text_tertiary::Colorant,
    amber::Colorant,
    cyan::Colorant,
    magenta::Colorant,
    categorical::AbstractVector{<:Colorant},
    fontscale::Real = 1,
    mono_ticklabels::Bool = false,
)
    # Helper to scale a base font size by `fontscale`, keeping the result a plain `Float64`.
    fs(base) = base * fontscale

    # Tick labels use the monospace font when requested (fixed-width digits align numeric
    # ticks); otherwise `:regular` keeps them in the body font inherited from `fonts`.
    ticklabelfont = mono_ticklabels ? _font("IBMPlexMono-Regular.ttf") : :regular

    return Theme(;

        # == Figure ========================================================================
        backgroundcolor = bg,
        textcolor       = text_primary,

        # Fonts are referenced by absolute path to the `.ttf` files bundled with this package
        # (see `assets/fonts/`), so text renders identically regardless of which fonts are
        # installed on the host system. IBM Plex Sans matches the INPE presentation template
        # and ships all four styles, so each maps to its dedicated weight.
        fonts = Attributes(;
            regular     = _font("IBMPlexSans-Regular.ttf"),     # ................ body text
            bold        = _font("IBMPlexSans-Bold.ttf"),        # ......... titles / display
            italic      = _font("IBMPlexSans-Italic.ttf"),
            bold_italic = _font("IBMPlexSans-BoldItalic.ttf"),
        ),

        # == Axis ==========================================================================
        Axis = Attributes(;
            backgroundcolor = card,  # slightly offset from figure bg,
                                     # keeps the plot area from blending into it.

            # Spines match the tick color so they read as subtle chrome rather than an accent.
            # Makie has no single `spinecolor`; each side must be set individually, otherwise
            # the spines fall back to the default black (wrong on the dark background).
            leftspinecolor     = text_secondary,
            rightspinecolor    = text_secondary,
            bottomspinecolor   = text_secondary,
            topspinecolor      = text_secondary,
            leftspinevisible   = true,
            rightspinevisible  = true,
            bottomspinevisible = true,
            topspinevisible    = true,
            spinewidth         = 1.5,

            # Grid.
            xgridcolor        = separator,
            ygridcolor        = separator,
            xgridwidth        = 1.2,
            ygridwidth        = 1.2,
            xgridstyle        = :dash,
            ygridstyle        = :dash,
            xminorgridcolor   = separator,
            yminorgridcolor   = separator,
            xminorgridwidth   = 0.8,
            yminorgridwidth   = 0.8,
            xminorgridstyle   = :dot,
            yminorgridstyle   = :dot,
            xminorgridvisible = true,
            yminorgridvisible = true,
            xminorticks       = IntervalsBetween(2),
            yminorticks       = IntervalsBetween(2),

            # Ticks.
            xtickcolor         = text_secondary,
            ytickcolor         = text_secondary,
            xtickwidth         = 1.0,
            ytickwidth         = 1.0,
            xminortickcolor    = text_tertiary,
            yminortickcolor    = text_tertiary,
            xminortickwidth    = 0.75,
            yminortickwidth    = 0.75,
            xminorticksize     = 3.0,
            yminorticksize     = 3.0,
            xminorticksvisible = true,
            yminorticksvisible = true,
            xticklabelcolor    = text_secondary,
            yticklabelcolor    = text_secondary,
            xticklabelsize     = fs(18),
            yticklabelsize     = fs(18),
            xticklabelfont     = ticklabelfont,
            yticklabelfont     = ticklabelfont,

            # Axis labels
            xlabelcolor   = text_primary,
            ylabelcolor   = text_primary,
            xlabelsize    = fs(21),
            ylabelsize    = fs(21),
            xlabelpadding = 10,
            ylabelpadding = 10,

            # Axis title
            titlecolor = text_primary,
            titlesize  = fs(24),
            titlefont  = :bold,
            titlegap   = 10,
        ),

        # == Legend ========================================================================
        Legend = Attributes(;
            backgroundcolor = card,
            framecolor      = border,
            framewidth      = 0.8,
            framevisible    = true,
            labelcolor      = text_primary,
            labelsize       = fs(18),
            titlecolor      = cyan,
            titlesize       = fs(20),
            titlefont       = :bold,
            padding         = (12f0, 12f0, 10f0, 10f0),
            patchsize       = (24f0, 12f0),
            rowgap          = 6,
            titlegap        = 8,
        ),

        # == Colorbar ======================================================================
        Colorbar = Attributes(;
            tickcolor      = text_secondary,
            ticklabelcolor = text_secondary,
            ticklabelsize  = fs(16),
            ticklabelfont  = ticklabelfont,
            labelcolor     = text_primary,
            labelsize      = fs(18),
            spinecolor     = cyan, # secondary accent, consistent with axes.
        ),

        # == Text ==========================================================================
        Text = Attributes(;
            color    = text_primary,
            fontsize = fs(18),
        ),

        # == Lines =========================================================================
        Lines = Attributes(;
            linewidth = 2.0,
            cycle     = Cycle([:color], covary = true),
        ),

        # == Scatter =======================================================================
        Scatter = Attributes(;
            markersize   = 8,
            strokewidth  = 0.5,
            strokecolor  = bg,
            cycle        = Cycle([:color, :marker], covary = true),
        ),

        # == BarPlot =======================================================================
        BarPlot = Attributes(;
            gap         = 0.15,
            strokewidth = 0,
            cycle       = Cycle([:color], covary = true),
        ),

        # == Heatmap =======================================================================
        Heatmap = Attributes(;
            colormap = :viridis,
        ),

        # == Color and Marker Cycle ========================================================
        palette = Attributes(;
            color     = categorical,
            marker    = [:circle, :rect, :diamond, :utriangle, :cross, :star5],
            linestyle = [:solid, :dash, :dot, :dashdot],
        ),
    )
end
