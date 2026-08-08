module SatelliteToolboxMakieExt

using Makie
using Makie: Colorant
using SatelliteToolbox
import SatelliteToolbox: makie_theme, makie_palette

############################################################################################
#                                        Constants                                         #
############################################################################################

include("./constants.jl")

############################################################################################
#                                         Includes                                         #
############################################################################################

include("./theme.jl")

end # module SatelliteToolboxMakieExt

