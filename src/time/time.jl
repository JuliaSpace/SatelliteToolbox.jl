#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# INPE - Instituto Nacional de Pesquisas Espaciais
# ETE  - Engenharia e Tecnologia Espacial
# DSE  - Divisão de Sistemas Espaciais
#
# Author: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions related Data and Time conversion.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
#   [2] https://support.microsoft.com/en-us/help/214019/method-to-determine-whether-a-year-is-a-leap-year
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-04-09: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export is_leap_year

################################################################################
#                                  Functions
################################################################################

"""
### function is_leap_year(year::Int)

Check if the year `year` is a leap year.

##### Args

* year: Year, must be > 0.

##### Returns

* **TRUE**: `year` is a leap year.
* **FALSE**: `year` is not a leap year.

##### Remarks

This algorithm was based on [3].

"""

function is_leap_year(year::Int)
    # Check if `year` is positive. This algorithm does not handle negative
    # years.
    (year < 0) && throw(ArgumentError("The year must be positive."))

    if (year % 4) == 0
        if (year % 100) == 0
            if (year % 400) == 0
                return true
            else
                return false
            end
        else
            return true
        end
    else
        return false
    end
end

