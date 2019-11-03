module SatelliteToolboxTLE

using Crayons, Dates, Parameters, Printf

import Base: show

################################################################################
#                                  Constants
################################################################################

# Escape sequences related to the crayons.
const _d = Crayon(reset = true)
const _b = crayon"bold"
const _g = crayon"bold green"
const _u = crayon"bold blue"
const _y = crayon"bold yellow"

################################################################################
#                                    Types
################################################################################

include("types.jl")

################################################################################
#                                   Includes
################################################################################

include("main.jl")
include("private.jl")

end
