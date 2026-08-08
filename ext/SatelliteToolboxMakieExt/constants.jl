## Description #############################################################################
#
# Definition of constants.
#
############################################################################################

# == Dark Theme Colors =====================================================================

const NAVY_PRIMARY   = colorant"#0A1929" # ................................ slide background
const NAVY_DEEP      = colorant"#04101F" # ................. darker half of split background
const NAVY_CARD      = colorant"#0D2137" # ..... slightly raised surface (legends, tooltips)
const SEPARATOR_DARK = colorant"#162940" # ................. grid lines, horizontal dividers
const BORDER_DARK    = colorant"#1E3A5F" # ....................... axis spines, box outlines

const TEXT_PRIMARY_DARK   = colorant"#F1F5F9" # ..................... body text, axis labels
const TEXT_SECONDARY_DARK = colorant"#94A3B8" # ............ captions, tick labels, metadata
const TEXT_TERTIARY_DARK  = colorant"#64748B" # ..................... minor tick labels only

const AMBER_DARK   = colorant"#F59E0B" # . primary (structural chrome, critical inline text)
const CYAN_DARK    = colorant"#38BDF8" # ......... secondary (data series, section headings)
const MAGENTA_DARK = colorant"#F472B6" # ..... tertiary (events, anomalies, second emphasis)

# Order optimized for distinguishability on dark backgrounds. Safe for deuteranopia and
# protanopia; avoid pairing green + coral as a critical duo.
const CATEGORICAL_DARK = [
    colorant"#38BDF8",   # 1  cyan
    colorant"#FBBF24",   # 2  amber
    colorant"#34D399",   # 3  green
    colorant"#F472B6",   # 4  magenta
    colorant"#A78BFA",   # 5  violet
    colorant"#F87171",   # 6  coral
]

# == Light Theme Colors ====================================================================

const SURFACE         = colorant"#FFFFFF" # ............................... slide background
const SURFACE_ALT     = colorant"#EEF4FB" # ............... lighter half of split background
const SURFACE_CARD    = colorant"#F7FAFD" # .. slightly recessed surface (legends, tooltips)
const SEPARATOR_LIGHT = colorant"#DDE8F5" # ................ grid lines, horizontal dividers
const BORDER_LIGHT    = colorant"#CBD5E1" # ...................... axis spines, box outlines

const TEXT_PRIMARY_LIGHT   = colorant"#0A1929" # .................... body text, axis labels
const TEXT_SECONDARY_LIGHT = colorant"#334155" # ........... captions, tick labels, metadata
const TEXT_TERTIARY_LIGHT  = colorant"#64748B" # .................... minor tick labels only

const AMBER_LIGHT   = colorant"#D97706" # ......... primary (darkened for contrast on white)
const CYAN_LIGHT    = colorant"#0284C7" # ....... secondary (darkened for contrast on white)
const MAGENTA_LIGHT = colorant"#DB2777" # ........ tertiary (darkened for contrast on white)

# Same hues as CATEGORICAL_DARK, darkened ~2 stops for contrast on white.
const CATEGORICAL_LIGHT = [
    colorant"#0284C7",   # 1  cyan
    colorant"#D97706",   # 2  amber
    colorant"#059669",   # 3  green
    colorant"#DB2777",   # 4  magenta
    colorant"#7C3AED",   # 5  violet
    colorant"#DC2626",   # 6  coral
]
