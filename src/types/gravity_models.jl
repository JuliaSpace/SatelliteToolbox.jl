# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   Structures related to the gravity models.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export ICGEM, GravityModel_Coefs

"""
    ICGEM

Structure to store the information contained in ICGEM files.

"""
@with_kw struct ICGEM
    # Mandatory keywords.
    product_type::String
    modelname::String
    gravity_constant::Float64
    radius::Float64
    max_degree::Int
    errors::Symbol

    # Optional keywords.
    tide_system::Symbol
    norm::Symbol

    # Coefficients.
    Clm::Matrix{Float64}           = Matrix{Float64}(undef,0,0)
    Slm::Matrix{Float64}           = Matrix{Float64}(undef,0,0)
    sigmaC::Matrix{Float64}        = Matrix{Float64}(undef,0,0)
    sigmaC_formal::Matrix{Float64} = Matrix{Float64}(undef,0,0)
    sigmaS::Matrix{Float64}        = Matrix{Float64}(undef,0,0)
    sigmaS_formal::Matrix{Float64} = Matrix{Float64}(undef,0,0)
end

"""
    GravityModel_Coefs{T}

Structure to store the information about a gravity model.

"""
@with_kw struct GravityModel_Coefs{T}
    name::String
    μ::T
    R0::T
    legendre_norm::Symbol
    n_max::Int
    C::Matrix{T}
    S::Matrix{T}
end
