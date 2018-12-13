#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions to parse the ICGEM files.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] Barthelmes, F., and Förste, C (2011). The ICGEM-format. GFZ Postdam,
#       Departmant 1 "Geodesy and Remote Sensing".
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export parse_icgem

################################################################################
#                                Private Macros
################################################################################

"""
    macro _keyword_found(keyword, keywords_found, current_line)

Check if the `keyword` exists in the list `keywords_found`. If `true`, then
throw an error indicating that the problem occurred on the `current_line`.

"""
macro _keyword_found(keyword, keywords_found, current_line)
    quote
        local lk  = $(esc(keyword))
        local lkf = $(esc(keywords_found))
        local lcl = $(esc(current_line))

        (findfirst(x -> x == lk, lkf) != nothing) &&
        error("Invalid ICGEM format (line $lcl): The keyword \"$lk\" already exists!")
    end
end

"""
    macro _parse_float(input)

Parse the `input` to float substituting all `D`s and `d`s  to `e`, so that we
can convert numbers in FORTRAN format.

"""
macro _parse_float(input)
    quote
        local li = $(esc(input))

        data = replace(li, r"[D,d]" => "e")
        parse(Float64, data)
    end
end

################################################################################
#                               Public Functions
################################################################################

"""
    function parse_icgem(filename::AbstractString)

Parse the ICGEM file `filename` and return an instance of the structure `ICGEM`
with the parsed data.

"""
function parse_icgem(filename::AbstractString)
    # Open the file and find the header.
    file = open(filename, "r")

    # Parse the header
    # ==========================================================================

    header_line_start   = 1
    header_line_end     = 0
    current_line        = 0
    begin_of_head_found = false
    end_of_head_found   = false

    # We need to first check the position of the header due to the
    # `begin_of_head` keyword that can define where the header starts.
    while !eof(file)
        current_line      +=1
        tokens = split(readline(file))

        length(tokens) < 1 && continue

        # Search for the `begin_of_head`, which is optional and makes all
        # previous lines to be ignored.
        if tokens[1] == "begin_of_head"
            # We should not have two `begin_of_head`.
            begin_of_head_found &&
            error("Invalid ICGEM format: Two `begin_of_head` keywords were found!")

            header_line_start = current_line
            begin_of_head_found = true
            continue
        elseif tokens[1] == "end_of_head"
            end_of_head_found = true
            header_line_end   = current_line
            break
        end
    end

    # `end_of_head` keyword is mandatory.
    !end_of_head_found &&
    error("Invalid ICGEM format: The mandatory keyword `end_of_head` was not found!")

    # Rewind file to read again.
    seek(file, 0)

    current_line = 0

    product_type       = ""
    modelname          = ""
    gravity_constant   = 0.0
    radius             = 0.0
    max_degree         = 0
    errors             = :formal
    tide_system        = :zero_tide
    norm               = :fully_normalized
    keywords_found     = Symbol[]
    mandatory_keywords = (:product_type, :modelname, :gravity_constant, :radius,
                          :max_degree, :errors)

    while current_line < header_line_end
        current_line += 1
        tokens = split(readline(file))

        # If we do not have at least 2 tokens, then it is not a keyword.
        length(tokens) < 2 && continue

        if tokens[1] == "product_type"
            @_keyword_found(:product_type, keywords_found, current_line)

            product_type = String(tokens[2])

            product_type != "gravity_field" &&
            error("Only ICGEM files of type \"gravity_field\" are currently supported.")

            push!(keywords_found, :product_type)

        elseif tokens[1] == "modelname"
            @_keyword_found(:modelname, keywords_found, current_line)

            modelname = String(tokens[2])
            push!(keywords_found, :modelname)

        elseif tokens[1] == "earth_gravity_constant" ||
               tokens[1] == "gravity_constant"
            @_keyword_found(:gravity_constant, keywords_found, current_line)

            try
                gravity_constant = @_parse_float(tokens[2])
            catch e
                @error(e.msg)
                error("Invalid ICGEM format (line $current_line): Invalid float number in keyword \"$(tokens[1])\".")
            end

            push!(keywords_found, :gravity_constant)

        elseif tokens[1] == "radius"
            @_keyword_found(:radius, keywords_found, current_line)

            try
                radius = @_parse_float(tokens[2])
            catch e
                @error(e.msg)
                error("Invalid ICGEM format (line $current_line): Invalid float number in keyword \"radius\".")
            end

            push!(keywords_found, :radius)

        elseif tokens[1] == "max_degree"
            @_keyword_found(:max_degree, keywords_found, current_line)

            try
                max_degree = parse(Int, tokens[2])
            catch e
                @error(e.msg)
                error("Invalid ICGEM format (line $current_line): Invalid integer number in keyword \"max_degree\".")
            end

            push!(keywords_found, :max_degree)

        elseif tokens[1] == "errors"
            @_keyword_found(:errors, keywords_found, current_line)

            errors = Symbol(tokens[2])
            push!(keywords_found, :errors)

            if !(errors in [:no, :calibrated, :formal, :calibrated_and_formal])
                @warn("Invalid ICGEM format: The error type \"$(tokens[2])\" is not valid! Assuming \"formal\".")
                errors = :formal
            end

        elseif tokens[1] == "tide_system"
            @_keyword_found(:tide_system, keywords_found, current_line)

            tide_system = Symbol(tokens[2])
            push!(keywords_found, :tide_system)

            if !(tide_system in [:zero_tide, :tide_free, :unknown])
                @warn("Invalid ICGEM format: The tide system \"$(tokens[2])\" is not valid! Assuming \"unknown\".")
                tide_system = :unknown
            end

        elseif tokens[1] == "norm"
            @_keyword_found(:norm, keywords_found, current_line)

            norm = Symbol(tokens[2])
            push!(keywords_found, :norm)

            if !(norm in [:fully_normalized, :unnormalized])
                @warn("Invalid ICGEM format: The norm \"$(tokens[2])\" is not valid! Assuming \"fully_normalized\".")
                norm = :fully_normalized
            end
        end
    end

    # Check if all mandatory keywords were found.
    found_mandatory_keywords::NTuple{6,Bool} =
        map(x->x in keywords_found, mandatory_keywords)
    all_mand_keywords_found  = true

    for i in eachindex(found_mandatory_keywords)
        if found_mandatory_keywords[i] == false
            @error("Invalid ICGEM format: The mandatory keyword \"$(mandatory_keywords[i])\" were not found!")
            all_mand_keywords_found = false
        end
    end

    !all_mand_keywords_found &&
    error("Invalid ICGEM format: Some mandatory keywords were not found!")

    # Parse the coefficients
    # ==========================================================================

    # Read the entire file to improve performance.
    raw = readdlm(filename; skipstart = header_line_end)

    # We currently only support the `gfc` coefficients.
    raw_gfc   = @views raw[ raw[:,1] .== "gfc", :]
    raw_nogfc = @views raw[ raw[:,1] .!= "gfc", :]

    !isempty(raw_nogfc) &&
    @warn("Unsupported ICGEM format: Only the keywords \"gfc\" are currently supported.")

    # Parse: gfc
    # ==========

    # Check if the number of columns is correct.
    if errors == :no
        size(raw_gfc,2) != 5 &&
        error("Invalid ICGME format: The coefficients table must have 5 columns because the keyword errors is \"$errors\".")
    elseif (errors == :calibrated) || (errors == :formal)
        size(raw_gfc,2) != 7 &&
        error("Invalid ICGME format: The coefficients table must have 7 columns because the keyword errors is \"$errors\".")
    else
        size(raw_gfc,2) != 9 &&
        error("Invalid ICGME format: The coefficients table must have 7 columns because the keyword errors is \"$errors\".")
    end

    # Check if there is any number to be converted. In this case, we assume that
    # the number is written in FORTRAN format.
    raw_gfc_coefs = @view(raw_gfc[:,4:end])

    for i in eachindex(raw_gfc_coefs)
        if !(typeof(raw_gfc_coefs[i]) <: AbstractFloat)
            try
                data             = replace(raw_gfc_coefs[i], r"[D,d]" => "e")
                raw_gfc_coefs[i] = parse(Float64, data)
            catch e
                error("Invalid ICGME format: Could not convert the coefficient \"$(raw_gfc_coefs[i])\" to a Float64.")
            end
        end
    end

    # The file will be parsed using Sparse matrices because it is much easier.
    # Perhaps, there is a method with better performance, but this function is
    # called only once per execution.
    @views I = map(x->Int(x)+1, raw_gfc[:,2])
    @views J = map(x->Int(x)+1, raw_gfc[:,3])

    @views Clm = sparse(I, J, convert(Vector{Float64}, raw_gfc[:,4]))
    @views Slm = sparse(I, J, convert(Vector{Float64}, raw_gfc[:,5]))

    if errors == :no
        return ICGEM(product_type     = product_type,
                     modelname        = modelname,
                     gravity_constant = gravity_constant,
                     radius           = radius,
                     max_degree       = max_degree,
                     errors           = errors,
                     tide_system      = tide_system,
                     norm             = norm,
                     Clm              = Clm,
                     Slm              = Slm
                    )

    elseif (errors == :calibrated) || (errors == :formal)
        @views sigmaC = sparse(I, J, convert(Vector{Float64}, raw_gfc[:,6]))
        @views sigmaS = sparse(I, J, convert(Vector{Float64}, raw_gfc[:,7]))

        ICGEM(product_type     = product_type,
              modelname        = modelname,
              gravity_constant = gravity_constant,
              radius           = radius,
              max_degree       = max_degree,
              errors           = errors,
              tide_system      = tide_system,
              norm             = norm,
              Clm              = Clm,
              sigmaC           = sigmaC,
              Slm              = Slm,
              sigmaS           = sigmaS
             )

    else
        @views sigmaC_cal    = sparse(I, J, convert(Vector{Float64}, raw_gfc[:,6]))
        @views sigmaS_cal    = sparse(I, J, convert(Vector{Float64}, raw_gfc[:,7]))
        @views sigmaC_formal = sparse(I, J, convert(Vector{Float64}, raw_gfc[:,8]))
        @views sigmaS_formal = sparse(I, J, convert(Vector{Float64}, raw_gfc[:,9]))

        ICGEM(product_type     = product_type,
              modelname        = modelname,
              gravity_constant = gravity_constant,
              radius           = radius,
              max_degree       = max_degree,
              errors           = errors,
              tide_system      = tide_system,
              norm             = norm,
              Clm              = Clm,
              sigmaC           = sigmaC_cal,
              sigmaC_formal    = sigmaC_formal,
              Slm              = Slm,
              sigmaS           = sigmaS_cal,
              sigmaS_formal    = sigmaS_formal
             )
    end
end

