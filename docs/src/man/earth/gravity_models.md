Gravity models
==============

```@meta
CurrentModule = SatelliteToolbox
DocTestSetup = quote
    using SatelliteToolbox
end
```

SatelliteToolbox.jl supports gravity models in which the coefficients are
provided in an [ICGEM](http://icgem.gfz-potsdam.de/home) file.

Before using a gravity model, we must load the coefficients from the input file.
This process can be accomplished by the function `parse_icgem`:

```julia
julia> egm96_model = parse_icgem("EGM96.gfc")
ICGEM
  product_type: String "gravity_field"
  modelname: String "EGM96"
  gravity_constant: Float64 3.986004415e14
  radius: Float64 6.3781363e6
  max_degree: Int64 360
  errors: Symbol formal
  tide_system: Symbol tide_free
  norm: Symbol fully_normalized
  Clm: Array{Float64}((361, 361)) [1.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; -6.34777471268e-11 5.36287680346e-11 … 6.27016396298e-11 0.0; 7.95017157688e-12 1.10754853834e-12 … 1.83971631467e-11 -4.47516389678e-25]
  Slm: Array{Float64}((361, 361)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 2.58517700048e-11 … -6.18987012099e-11 0.0; 0.0 -6.27959348935e-11 … -3.10123632209e-11 -8.30224945525e-11]
  sigmaC: Array{Float64}((361, 361)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 6.4425138e-11 3.8751548e-11 … 6.3040162e-11 0.0; 5.0033977e-11 5.0033977e-11 … 5.0033977e-11 5.0033977e-11]
  sigmaC_formal: Array{Float64}((0, 0)) Matrix{Float64}(undef, 0, 0)
  sigmaS: Array{Float64}((361, 361)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 3.8983376e-11 … 6.2932021e-11 0.0; 5.0033977e-11 5.0033977e-11 … 5.0033977e-11 5.0033977e-11]
  sigmaS_formal: Array{Float64}((0, 0)) Matrix{Float64}(undef, 0, 0)
```

It loads all the information in the ICGEM file and stores in an `ICGEM`
structure.

SatelliteToolbox.jl currently supports three functions using gravity models:

1. Computing the gravitational potential at a specific location;
2. Computing the gravitational potential derivatives (spherical coordinates) at
   a specific location; and
3. Computing the gravitational acceleration at a specific location.

For all those tasks, we initially need to process the information in the `ICGEM`
structure into an object of type `GravityModel_Coefs`. This action is
accomplished by the function `create_gravity_model_coefs`:

```julia
julia> egm96_coefs = create_gravity_model_coefs(egm96_model)
GravityModel_Coefs{Float64}
  name: String "EGM96"
  μ: Float64 3.986004415e14
  R0: Float64 6.3781363e6
  legendre_norm: Symbol full
  n_max: Int64 360
  C: Array{Float64}((361, 361)) [1.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; -6.34777471268e-11 5.36287680346e-11 … 6.27016396298e-11 0.0; 7.95017157688e-12 1.10754853834e-12 … 1.83971631467e-11 -4.47516389678e-25]
  S: Array{Float64}((361, 361)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 2.58517700048e-11 … -6.18987012099e-11 0.0; 0.0 -6.27959348935e-11 … -3.10123632209e-11 -8.30224945525e-11]
```

Finally, the object `egm96_coefs` is used to perform the computations.

Notice that SatelliteToolbox.jl provides three built-in gravity models. Their
`GravityModel_Coefs` objects can be loaded using the function
`load_gravity_model(T)`, where `T` can be:

- `EGM96()` for the Earth Gravitational Model 1996 (EGM96);
- `JGM2()` for the Joint Earth Gravity Model 2 (JGM2); or
- `JGM3()` for the Joint Earth Gravity Model 3 (JGM3).

```jldoctest
julia> egm96_coefs = load_gravity_model(EGM96())
GravityModel_Coefs{Float64}
  name: String "EGM96"
  μ: Float64 3.986004415e14
  R0: Float64 6.3781363e6
  legendre_norm: Symbol full
  n_max: Int64 360
  C: Array{Float64}((361, 361)) [1.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; -6.34777471268e-11 5.36287680346e-11 … 6.27016396298e-11 0.0; 7.95017157688e-12 1.10754853834e-12 … 1.83971631467e-11 -4.47516389678e-25]
  S: Array{Float64}((361, 361)) [0.0 0.0 … 0.0 0.0; 0.0 0.0 … 0.0 0.0; … ; 0.0 2.58517700048e-11 … -6.18987012099e-11 0.0; 0.0 -6.27959348935e-11 … -3.10123632209e-11 -8.30224945525e-11]
```

## Computing the gravitational potential

The gravitational potential [J/kg] at a specific location can be computed using
the function `compute_U`. The function signature is:

```julia
compute_U(gm_coefs::GravityModel_Coefs{T}, r::AbstractVector, n_max::Number = -1, m_max::Number = -1) where T<:Number
```

where:

- `gm_coefs` is the object that holds the gravity model coefficients (see
  [`GravityModel_Coefs`](@ref)); and
- `r` is the position vector of the desired location represented in ITRF [m].

The arguments `n_max` and `m_max` can be used to limit the maximum degree and
order the algorithm will use when computing the spherical harmonics. All the
available information within the model is considered if they are negative.

```jldoctest
julia> egm96_coefs = load_gravity_model(EGM96());

julia> compute_U(egm96_coefs, [6378e3, 0, 0])
6.2530209756063506e7
```

## Computing the gravitational potential derivatives

The gravitational potential derivatives w.r.t. the spherical coordinates
(`∂U/∂r`, `∂U/∂ϕ`, `∂U/∂λ`) at specific location can be computed using
the function `compute_dU`. The function signature is:

```julia
compute_dU(gm_coefs::GravityModel_Coefs{T}, r::AbstractVector, n_max::Number = -1, m_max::Number = -1) where T<:Number
```

where:

- `gm_coefs` is the object that holds the gravity model coefficients (see
  [`GravityModel_Coefs`](@ref)); and
- `r` is the position vector of the desired location represented in ITRF [m].

The arguments `n_max` and `m_max` can be used to limit the maximum degree and
order the algorithm will use when computing the spherical harmonics. All the
available information within the model is considered if they are negative.

!!! info
    In this case, `ϕ` is the geocentric latitude and `λ` is the geocentric
    longitude.

```julia
julia> egm96_coefs = load_gravity_model(EGM96());

julia> compute_dU(egm96_coefs, [6378e3, 0, 0])
(-9.81470672242183, 50.712111877176156, -116.12345916122457)
```

## Computing the gravitational acceleration

The gravitational acceleration vector [m/s²] represented in the ITRF at specific
location can be computed using the function `compute_g`. The function signature
is:

```julia
compute_g(gm_coefs::GravityModel_Coefs{T}, r::AbstractVector, n_max::Number = -1, m_max::Number = -1) where T<:Number
```

- `gm_coefs` is the object that holds the gravity model coefficients (see
  [`GravityModel_Coefs`](@ref)); and
- `r` is the position vector of the desired location represented in ITRF [m].

The arguments `n_max` and `m_max` can be used to limit the maximum degree and
order the algorithm will use when computing the spherical harmonics. All the
available information within the model is considered if they are negative.

!!! info
    Notice that this function computes the **gravitational acceleration**.
    Hence, the acceleration due to Earth's rotation rate **is not** included.

```jldoctest
julia> egm96_coefs = load_gravity_model(EGM96());

julia> compute_g(egm96_coefs, [6378e3, 0, 0])
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -9.81470672242183
 -1.8206876632365096e-5
  7.95109938494452e-6
```
