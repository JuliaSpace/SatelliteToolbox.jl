#== # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   Functions to verify if the satellite is within the station visibility
#   circle.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # ==#

export ground_station_accesses, ground_station_gaps,
       list_ground_station_accesses, list_ground_station_gaps,
       ground_station_visible

"""
    ground_station_accesses(orbp, vrs_e,     Δt, ECI, ECEF, vargs...; kwargs...)
    ground_station_accesses(orbp, [(WGS84)], Δt, ECI, ECEF, vargs...; kwargs...)

Compute the accesses of a satellite with orbit propagator `orbp` (see
`init_orbit_propagator`) to the ground stations defined in the vector `vrs_e`.
The analysis interval begins in the propagator epoch and lasts `Δt` [s].

The ground stations can be specified by an array of 3×1 vectors describing the
ground stations position in an ECEF frame `vrs_e` or by an array of tuples
containing the WGS84 position of each ground station `[(WGS84)]`:

    (latitude [rad], longitude [rad], altitude [m])

# Args

* `ECI`: Earth-Centered Inertial frame in which the state vector of the
         propagator is represented.
* `ECEF`: Earth-Centered, Earth-fixed frame to be used for the analysis. It
          must be the same frame used to compute the ground station position
          vector.
* `vargs...`: list of additional arguments to be passed to the function
              `rECItoECEF` when converting the ECI frame to the ECEF.

# Keywords

* `θ`: Minimum elevation angle for communication between the satellite and the
       ground stations [rad]. (**Default** = 10ᵒ)
* `reduction`: A function that receives a boolean vector with the visibility
               between the satellite and each ground station. It must return a
               boolean value indicating if the access must be computed or not.
               This is useful to merge access time between two or more stations.
               (**Default** = `v->|(v...)` *i.e.* compute the access if at least
               one ground station is visible)

"""
ground_station_accesses(orbp, gs_wgs84::Tuple, vargs...; kwargs...) =
    ground_station_accesses(orbp, [gs_wgs84], vargs...; kwargs...)

function ground_station_accesses(orbp, vgs_wgs84::AbstractVector{T}, vargs...;
                             kwargs...) where T<:Tuple

    vrs_e = [GeodetictoECEF(gs_wgs84...) for gs_wgs84 in vgs_wgs84]
    return ground_station_accesses(orbp, vrs_e, vargs...; kwargs...)
end

ground_station_accesses(orbp, rs_e::AbstractVector{T}, vargs...; kwargs...) where
    T<:Number = ground_station_accesses(orbp, [rs_e], vargs...; kwargs...)

function ground_station_accesses(orbp, vrs_e::AbstractVector{T}, Δt::Number,
                             ECI::Union{T_ECIs, T_ECIs_IAU_2006},
                             ECEF::Union{T_ECEFs, T_ECEFs_IAU_2006}, vargs...;
                             θ::Number = 10*pi/180,
                             reduction::Function = v->|(v...),
                             step::Number = 60) where T<:AbstractVector

    # Time vector of the analysis.
    t = 0:step:Δt

    # Get the epoch of the propagator.
    JD₀ = epoch(orbp)

    # Matrix that will contain the accesses.
    accesses = Matrix{DateTime}(undef,0,2)

    # State to help the computation.
    state = :initial

    # Lambda function to check if the ground station is visible.
    f(t)::Bool = begin
        o,r_i,v_i  = propagate!(orbp, t)
        r_e        = rECItoECEF(DCM, ECI, ECEF, o.t, vargs...)*r_i
        visibility = [ground_station_visible(r_e, rs_e, θ) for rs_e in vrs_e]

        return reduction(visibility)
    end

    access_beg = DateTime(now())
    access_end = DateTime(now())

    for k in t
        # Check if the ground station is visible.
        gs_visible = f(k)

        # Handle the initial case.
        if state == :initial
            if gs_visible
                access_beg = JDtoDate(DateTime, JD₀)
                state = :visible
            else
                state = :not_visible
            end
        # Handle transitions.
        elseif (state == :not_visible) && gs_visible
            # Refine to find the edge.
            k₀ = k - step
            k₁ = k
            kc = find_crossing(f, k₀, k₁, false, true)

            state = :visible
            access_beg = JDtoDate(DateTime, JD₀ + kc/86400)

        elseif (state == :visible) && !gs_visible
            # Refine to find the edge.
            k₀ = k - step
            k₁ = k
            kc = find_crossing(f, k₀, k₁, true, false)

            state = :not_visible
            access_end = JDtoDate(DateTime, JD₀ + kc/86400)

            accesses = vcat(accesses, [access_beg access_end])
        end
    end

    # If the analysis finished during an access, then just add the end of the
    # interval as the end of the access.
    if state == :visible
        access_end = JDtoDate(DateTime, JD₀ + Δt/86400)
        accesses   = vcat(accesses, [access_beg access_end])
    end

    return accesses
end

"""
    ground_station_gaps(args...; kwargs...)

Compute the gaps between the accesses of ground stations. The arguments and
keywords are the same as the ones used in the function
`ground_station_accesses`.

Notice that the gap analysis starts in the orbit propagator epoch and ends in
the instant defined by the argument `Δt`.

"""
function ground_station_gaps(orbp, args...; kwargs...)
    # Get the epoch of the propagator.
    JD₀ = epoch(orbp)
    DT₀ = JDtoDate(DateTime, JD₀)

    # Compute the list of ground station accesses.
    accesses = ground_station_accesses(orbp, args...; kwargs...)

    # Get the analysis period, which must be the second argument in `args...`.
    Δt = args[2]

    # Compute the last propagation instant.
    JD₁ = JD₀ + Δt/86400
    DT₁ = JDtoDate(DateTime, JD₁)

    # Compute the gaps between accesses.
    gaps = Matrix{DateTime}(undef,0,2)

    # Check if the simulation did not start under the visibility of a ground
    # station.
    if accesses[1,1] != DT₀
        gaps = vcat(gaps, [JDtoDate(DateTime, JD₀) accesses[1,1]])
    end

    # Check if the simulation did not end under the visibility of a ground
    # station.
    if accesses[end,2] != DT₁
        aux  = JDtoDate(DateTime, JD₁)
        gaps = vcat(gaps, [ accesses[:,2] vcat(accesses[2:end,1], aux) ])
    else
        gaps = vcat(gaps, [ accesses[1:end-1,2] accesses[2:end,1] ])
    end

    return gaps
end

"""
    list_ground_station_accesses(io, vargs...; kwargs...)

Print the ground station accesses to the io `io`. The arguments `vargs...` and
keywords `kwargs...` are those of the function `ground_station_accesses`.

Additionally, the following keywords can be used to modify the behavior of this
function:

* `format`: If `:pretty`, then a formatted table will be printed. If `:csv`,
            then the access data will be printed using the CSV format.
            (**Default** = `:pretty`)
* `time_scale`: Select the time scale of the access duration (`:s` for seconds,
                `:m` for minutes, and `:h` for hours). (**Default** = `:m`)

"""
list_ground_station_accesses(vargs...; kwargs...) =
    list_ground_station_accesses(stdout, vargs...; kwargs...)

function list_ground_station_accesses(io::IO, vargs...; format = :pretty,
                          time_scale::Symbol = :m, kwargs...)

    # Compute the accesses.
    accesses = ground_station_accesses(vargs...; kwargs...)

    # If no access was found, then inform the user and return.
    if isempty(accesses)
        println("No access found!")
        return nothing
    end

    # Check which scale to use when printing the access duration.
    if time_scale == :s
        ts    = 1
        label = "s"
    elseif time_scale == :m
        ts    = 60
        label = "min"
    else
        ts    = 3600
        label = "h"
    end

    # Compute the accesses duration.
    access_time = map((b,e)->(e-b).value/1000/ts,
                       @view(accesses[:,1]), @view(accesses[:,2]))
    num         = length(access_time)

    # Assemble the information to be printed.
    header = ["Access #" "Beginning (UTC)" "End (UTC)" "Duration [$label]"]
    mp     = hcat(string.(collect(1:1:num)), accesses, access_time)

    # Print the information depending on the format.

    if format == :csv
        writedlm(io, vcat(header,mp), ",")
    else
        # Statistics.
        min_access, id_min = findmin(access_time)
        max_access, id_max = findmax(access_time)
        mean_access        = mean(access_time)
        total_access       = sum(access_time)

        # Mark the minimum and maximum accesses.
        mp[id_min, 1] = "MIN " * mp[id_min,1]
        mp[id_max, 1] = "MAX " * mp[id_max,1]

        # Add the statistics to the end of the table.
        mp = vcat(mp, ["" "" "         Minimum access"  min_access;
                       "" "" "            Mean access"  mean_access;
                       "" "" "         Maximum access"  max_access;
                       "" "" "           Total access"  total_access;])

        # Highlighters for the minimum and maximum accesses.
        hl_min = Highlighter((data,i,j)->i == id_min, crayon"bold red")
        hl_max = Highlighter((data,i,j)->i == id_max, crayon"bold blue")

        # Highlighter for statistics labels.
        hl = Highlighter((data,i,j)->(j == 3) && (i-num ∈ (1,2,3,4)),
                         crayon"bold")

        # Print.
        pretty_table(io, mp, header; alignment = [:r, :l, :l, :r],
                     crop = :horizontal, formatter = ft_printf("%.3f",[4]),
                     hlines = [num,num+3], highlighters = (hl,hl_min,hl_max))
    end

    return nothing
end

"""
    list_ground_station_gaps(io, vargs...; kwargs...)

Print the ground station gaps to the io `io`. The arguments `vargs...` and
keywords `kwargs...` are those of the function `ground_station_gaps`.

Additionally, the following keywords can be used to modify the behavior of this
function:

* `format`: If `:pretty`, then a formatted table will be printed. If `:csv`,
            then the access data will be printed using the CSV format.
            (**Default** = `:pretty`)
* `time_scale`: Select the time scale of the access duration (`:s` for seconds,
                `:m` for minutes, and `:h` for hours). (**Default** = `:m`)

"""
list_ground_station_gaps(vargs...; kwargs...) =
    list_ground_station_gaps(stdout, vargs...; kwargs...)

function list_ground_station_gaps(io::IO, vargs...; format = :pretty,
                                  time_scale::Symbol = :m, kwargs...)

    # Compute the gaps.
    gaps = ground_station_gaps(vargs...; kwargs...)

    # If no gap was found, then inform the user and return.
    if isempty(gaps)
        println("No gap found!")
        return nothing
    end

    # Check which scale to use when printing the gap duration.
    if time_scale == :s
        ts    = 1
        label = "s"
    elseif time_scale == :m
        ts    = 60
        label = "min"
    else
        ts    = 3600
        label = "h"
    end

    # Compute the gaps duration.
    gap_time = map((b,e)->(e-b).value/1000/ts, @view(gaps[:,1]), @view(gaps[:,2]))
    num         = length(gap_time)

    # Assemble the information to be printed.
    header = ["   Gap #" "Beginning (UTC)" "End (UTC)" "Duration [$label]"]
    mp     = hcat(string.(collect(1:1:num)), gaps, gap_time)

    # Print the information depending on the format.

    if format == :csv
        writedlm(io, vcat(header,mp), ",")
    else
        # Statistics.
        min_gap, id_min = findmin(gap_time)
        max_gap, id_max = findmax(gap_time)
        mean_gap        = mean(gap_time)
        total_gap       = sum(gap_time)

        # Mark the minimum and maximum gaps.
        mp[id_min, 1] = "MIN " * mp[id_min,1]
        mp[id_max, 1] = "MAX " * mp[id_max,1]

        # Add the statistics to the end of the table.
        mp = vcat(mp, ["" "" "            Minimum gap"  min_gap;
                       "" "" "               Mean gap"  mean_gap;
                       "" "" "            Maximum gap"  max_gap;
                       "" "" "              Total gap"  total_gap;])

        # Highlighters for the minimum and maximum gaps.
        hl_min = Highlighter((data,i,j)->i == id_min, crayon"bold red")
        hl_max = Highlighter((data,i,j)->i == id_max, crayon"bold blue")

        # Highlighter for statistics labels.
        hl = Highlighter((data,i,j)->(j == 3) && (i-num ∈ (1,2,3,4)),
                         crayon"bold")

        # Print.
        pretty_table(io, mp, header; alignment = [:r, :l, :l, :r],
                     crop = :horizontal, formatter = ft_printf("%.3f",[4]),
                     hlines = [num,num+3], highlighters = (hl,hl_min,hl_max))
    end

    return nothing
end

"""
    ground_station_visible(r_e::AbstractVector, rs_e::AbstractVector, θ::Number)

Check if the satellite with position vector `r_e` (ECEF) is inside the
visibility circle of a ground station with position vector `rs_e` (ECEF) and a
minimum elevation angle of `θ` [rad].

Notice that `r_e` and `rs_e` must be represented in the same ECEF frame, and
must have the same unit.

Returns `true` if the satellite is inside the visibility circle, or `false`
otherwise.

"""
function ground_station_visible(r_e::AbstractVector, rs_e::AbstractVector,
                                θ::Number)
    # Check if the satellite is within the visibility circle of the station.
    dr_e = r_e - rs_e
    cos_beta = dot( dr_e/norm(dr_e), rs_e/norm(rs_e) )

    return cos_beta > cos(π/2-θ)
end


"""
    ground_station_visible(r_e::AbstractVector, lat_s::Number, lon_s::Number, h_s::Number, θ::Number)

Check if the satellite with position vector `r_e` (ECEF) is inside the
visibility circle of a ground station with latitude `lat_s` [rad], longitude
`lon_s` [rad], altitude `h_s` (WGS-84), and a minimum elevation angle of
`θ` [rad].

Notice that the units of `r_e` and `h_s` must be the same.

Returns `true` if the satellite is inside the visibility circle, or `false`
otherwise.

"""
function ground_station_visible(r_e::AbstractVector, lat_s::Number,
                                lon_s::Number, h_s::Number, θ::Number)
    # Convert the ground station LLA to the ECEF frame.
    rs_e = GeodetictoECEF(lat_s, lon_s, h_s)

    return ground_station_visible(r_e, rs_e, θ)
end
