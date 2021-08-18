#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Definition of the Ellipsoid struct.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html
#       Accessed on 2017-08-07.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export Ellipsoid

"""
    Ellipsoid{T}

Ellipsoid of rotation to be used for geocentric, geodetic and ecef transformations.

!!! note
    The constructor only accepts the fields `a` and `f`, with the other fields pre-computed 
    automatically from those two

# Fields
- `a` : Semi-major axis [m].
- `f` : Flattening of the ellipsoid.
- `b` : Semi-minor axis [m]
- `e²` : Eccentricity squared
- `el²` : Second Eccentricity squared
"""
struct Ellipsoid{T}
    ## Main Variables
    a::T # Semi-major axis in [m]
    f::T # Flattening of the ellipsoid
    ## Auxiliary variables, pre-computed just for convenience
    b::T # Semi-minor axis in [m]
    e²::T # Eccentricity squared
    el²::T # Second eccentricity squared
    
    ## Constructor
    function Ellipsoid(a,f)
        @assert f < 1 "The flattening should be lower than 1"
        b = a * (1 - f)
        e² = f * (2 - f)
        el² = e² / (1 - e²)
        new{typeof(el²)}(a,f,b,e²,el²)
    end
end
 
"""
    Ellipsoid(a,f)
    Ellipsoid{T}(a,f)

Generate an ellipsoid of rotation ([`Ellipsoid`](@ref)) as a function of the semi-major axis in [m] and the flattening
"""
Ellipsoid{T}(a,f) where T = Ellipsoid(T(a),T(f))

# Define the default ellipsoid based on WGS-84 [2]
const wgs84_ellipsoid = Ellipsoid(6378137.0, 1/298.257223563)
