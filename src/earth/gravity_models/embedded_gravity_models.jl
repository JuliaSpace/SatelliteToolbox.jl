#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   This function reads the embedded gravity models.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] http://icgem.gfz-potsdam.de/home
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export load_gravity_model
export EGM96, JGM2, JGM3

################################################################################
#                               Public Functions
################################################################################

"""
    function load_gravity_model(T)

Load an embedded gravity model coefficients `T`. The current supported values
for `T` are:

| `T`       | Model Name                     | Maximum Degree |
|:---------:|:-------------------------------|:---------------|
| `EGM96()` | Earth Gravitational Model 1996 | 360            |
| `JGM2()`  | Joint Earth Gravity Model 2    | 70             |
| `JGM3()`  | Joint Earth Gravity Model 3    | 70             |
|-----------|--------------------------------|----------------|

For other models, you can downlad the `gfc` file at

    http://icgem.gfz-potsdam.de/home

and load it using the function `parse_gfc`.

# Args

* `T`: Embedded model to be loaded.

# Returns

An instance of the structure `GravityModel_Coefs` with the parsed values.

"""
function load_gravity_model(::Type{Val{:egm96}})
    dir      = @__DIR__
    filename = dir * "/coefficients/EGM96.gfc"

    return parse_gfc(filename)
end

function load_gravity_model(::Type{Val{:jgm2}})
    dir      = @__DIR__
    filename = dir * "/coefficients/JGM2.gfc"

    return parse_gfc(filename)
end

function load_gravity_model(::Type{Val{:jgm3}})
    dir      = @__DIR__
    filename = dir * "/coefficients/JGM3.gfc"

    return parse_gfc(filename)
end

################################################################################
#                             Auxiliary Functions
################################################################################

@inline EGM96() = Val{:egm96}
@inline JGM2()  = Val{:jgm2}
@inline JGM3()  = Val{:jgm3}
