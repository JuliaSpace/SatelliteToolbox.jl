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
#   Functions related with the IERS (International Earth Orientation and
#   Reference Systems Service) EOP (Earth Orientation Parameters) data.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Vallado, D. A (2013). Fundamentals of Astrodynamics and Applications.
#       Microcosm Press, Hawthorn, CA, USA.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Changelog
#
# 2018-04-10: Ronan Arraes Jardim Chagas <ronan.arraes@inpe.br>
#   Initial version.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

using HTTP

export get_iers_eop, get_eop_iau_1980, get_eop_iau_2000A

"""
### function get_iers_eop(data_type::Symbol = :IAU2000A)

Download and parse the IERS EOP C04 data. The data type is specified by
`data_type` symbol. Supported values are:

* `IAU1980`: Get IERS EOP C04 IAU1980 data.
* `IAU2000A`: Get IERS EOP C04 IAU2000A data.

##### Args

* data_type: (OPTIONAL) Symbol that defines the type of data that will be
             downloaded (**DEFAULT** = `:IAU2000A`).

##### Returns

A matrix with the IERS EOP C4 data.

"""

function get_iers_eop(data_type::Symbol = :IAU2000A)
    if data_type == :IAU1980
        return get_iers_eop_iau_1980()
    elseif data_type == :IAU2000A
        return get_iers_eop_iau_2000A()
    else
        throw(ArgumentError("Unknow EOP type. It must be :IAU1980 or :IAU2000A."))
    end
end

"""
### function get_iers_eop_iau_1980(url::String = "https://datacenter.iers.org/eop/-/somos/5Rgv/latest/213")

Get the IERS EOP C04 IAU1980 data from the URL `url`.

##### Args

* url: (OPTIONAL) Url in which the coefficients will be obtained (**DEFAULT** =
       https://datacenter.iers.org/eop/-/somos/5Rgv/latest/213)

##### Returns

A matrix with the IERS EOP C4 IAU1980 data.

"""

function get_iers_eop_iau_1980(url::String = "https://datacenter.iers.org/eop/-/somos/5Rgv/latest/213")
    r = HTTP.request("GET", url)

    # Parse the data removing the header.
    eop = readdlm(r.body; skipblanks=true, skipstart=14)

    # Check if the dimension seems correct.
    if size(eop,2) != 16
        error("The fetched data does not have the correct dimension. Please, check the URL.")
    end

    eop
end

"""
### function get_iers_eop_iau_2000A(url::String = "https://datacenter.iers.org/eop/-/somos/5Rgv/latest/214")

Get the IERS EOP C04 IAU2000A data from the URL `url`.

##### Args

* url: (OPTIONAL) Url in which the coefficients will be obtained (**DEFAULT** =
       https://datacenter.iers.org/eop/-/somos/5Rgv/latest/214)

##### Returns

A matrix with the IERS EOP C4 IAU2000A data.

"""

function get_iers_eop_iau_2000A(url::String = "https://datacenter.iers.org/eop/-/somos/5Rgv/latest/214")
    r = HTTP.request("GET", url)

    # Parse the data removing the header.
    eop = readdlm(r.body; skipblanks=true, skipstart=14)

    # Check if the dimension seems correct.
    if size(eop,2) != 16
        error("The fetched data does not have the correct dimension. Please, check the URL.")
    end

    eop
end
