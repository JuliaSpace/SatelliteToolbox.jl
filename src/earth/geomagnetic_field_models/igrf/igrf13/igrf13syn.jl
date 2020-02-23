# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
#
#   International Geomagnetic Field Model v13.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# References
#
#   [1] https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
#   [2] https://www.ngdc.noaa.gov/IAGA/vmod/igrf13.f
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export igrf13syn

include("./igrf13syn_coefs.jl")

"""
    igrf13syn(isv::Int, date::Number, itype::Int, alt::Number, colat::Number, elong::Number; show_warns = true)

This is a Julia implementation of the official IGRF source code, which was
written in Fortran [2]. The input and output variables are exactly the same as
the ones described in the function `igrf13syn` in [2].

# Args

* `isv`: `0` if main-field values are required, `1` if secular variation values
         are required.
* `date`: Year A.D.
* `itype`: `1` if geodetic (spheroid), `2` if geocentric (sphere).
* `alt`: Height above sea level [km] if `itype = 1`, or distance from the center of
         Earth [km] if `itype = 2` (must be > 3485 km).
* `colat`: Colatitude (0 - 180) [˚].
* `elong`: East-Longitude (0 - 360) [˚].

# Keywords

* `show_warns`: Show warnings about the data (**Default** = `true`).

# Returns

* The north component [nT] if `isv = 0`, or [nT/year] if `isv = 1`.
* The east component [nT] if `isv = 0`, or [nT/year] if `isv = 1`.
* The vertical component [nT] if `isv = 0`, or [nT/year] if `isv = 1`.
* The total intensity if `isv = 0`, or rubbish if `isv = 1`.

# Remarks

* The `date` must be greater or equal to 1900 and less than or equal 2030.
  Notice that a warning message is printed for dates grated than 2025.

"""
function igrf13syn(isv::Int, date::Number, itype::Int, alt::Number,
                   colat::Number, elong::Number; show_warns = true)

    # Check the data, since this model is valid for years between 1900 and 2025.
    ( (date < 1900) || (date > 2030) ) &&
    error("This IGRF version will not work for years outside the interval [1900, 2030).")

    # Warn the user that for dates after the year 2020 the accuracy maybe reduced.
    show_warns && (date > 2025) &&
    @warn("The magnetic field computed with this IGRF version may be of reduced accuracy for years greater than 2025.")

    # Declaration of variables
    gh          = gh_igrf13
    fn::Int     = 0
    gn::Int     = 0
    kmx::Int    = 0
    ll::Int     = 0
    nc::Int     = 0
    nmx::Int    = 0
    x::Float64  = 0.0
    y::Float64  = 0.0
    z::Float64  = 0.0
    t::Float64  = 0.0
    tc::Float64 = 0.0
    cl          = zeros(13)
    sl          = zeros(13)
    p           = zeros(105)
    q           = zeros(105)

    if date < 2020
        t  = 0.2*(date - 1900)
        ll = floor(Int,t)
        t  = t - ll

        # SH models before 1995.0 are only to degree 10.
        if date < 1995
            nmx = 10

            # nc = nmx*(nmx+2)
            nc  = 120

            ll  = nc*ll

            # kmx = (nmx+1)*(nmx+2)/2
            kmx = 66
        else
            nmx = 13

            # nc  = nmx*(nmx+2)
            nc  = 195

            ll  = floor(Int,0.2*(date - 1995))

            # 19 is the number of SH models that extend to degree 10.
            ll  = 120*19 + nc*ll

            # kmx = (nmx+1)*(nmx+2)/2
            kmx = 105
        end

        tc = 1 - t

        if isv == 1
            t  = +0.2
            tc = -0.2
        end
    else
        t  = date - 2020
        tc = 1.0

        if isv == 1
            t  = 1.0
            tc = 0.0
        end

        ll  = 3255
        nmx = 13

        # nc  = nmx*(nmx+2)
        nc = 195

        # kmx = (nmx+1)*(nmx+2)/2
        kmx = 105
    end

    r::Float64  = alt
    ct          = cos(colat*pi/180)
    st          = sin(colat*pi/180)
    cl[1]       = cos(elong*pi/180)
    sl[1]       = sin(elong*pi/180)
    cd::Float64 = 1.0
    sd::Float64 = 0.0
    l::Int      = 1
    m::Int      = 1
    n::Int      = 0

    if itype != 2
        # Conversion between geodetic and geocentric coordinates using WGS84
        # spheroid.
        a2    = 40680631.6
        b2    = 40408296.0
        one   = a2*st^2
        two   = b2*ct^2
        three = one + two
        rho   = sqrt(three)
        r     = sqrt(alt*(alt + 2*rho) + (a2*one + b2*two)/three)
        cd    = (alt + rho)/r
        sd    = (a2 - b2)/rho*ct*st/r
        one   = ct
        ct    = ct*cd -  st*sd
        st    = st*cd + one*sd
    end

    ratio = 6371.2/r
    rr    = ratio^2

    # Computation of Schmidt quasi-normal coefficients p and x(=q).
    # =============================================================

    p[1] = 1.0
    p[3] = st
    q[1] = 0.0
    q[3] = ct

    for k = 2:kmx
        # There is no need to check bounds here. The code guarantees that
        # everything is inside the right bounds. This increased the source code
        # performance by 13%.
        @inbounds begin
            if n < m
                m  = 0
                n  = n+1
                rr = rr*ratio
                fn = n
                gn = n-1
            end

            fm = m

            if (m == n)
                if k != 3
                    one   = sqrt(1 - 0.5/fm)
                    j     = k - n - 1
                    p[k]  = one*st*p[j]
                    q[k]  = one*(st*q[j] + ct*p[j])
                    cl[m] = cl[m-1]*cl[1] - sl[m-1]*sl[1]
                    sl[m] = sl[m-1]*cl[1] + cl[m-1]*sl[1]
                end
            else
                gmm   = m^2
                one   = sqrt(fn^2 - gmm)
                two   = sqrt(gn^2 - gmm)/one
                three = (fn + gn)/one
                i     = k - n
                j     = i - n + 1
                p[k]  = three*ct*p[i] - two*p[j]
                q[k]  = three*(ct*q[i] - st*p[i]) - two*q[j]
            end

            # Synthesis of x, y, and z in geocentric coordinates.
            lm  = ll + l
            one = (tc*gh[lm] + t*gh[lm+nc])*rr

            if m != 0
                two   = (tc*gh[lm+1] + t*gh[lm+nc+1])*rr
                three = one*cl[m] + two*sl[m]
                x     = x + three*q[k]
                z     = z - (fn + 1)*three*p[k]

                if st != 0
                    y = y + (one*sl[m] - two*cl[m])*fm*p[k]/st
                else
                    y = y + (one*sl[m] - two*cl[m])*q[k]*ct
                end

                l = l + 2
            else
                x = x + one*q[k]
                z = z - (fn + 1)*one*p[k]
                l = l + 1
            end

            m = m + 1
        end
    end

    # Conversion to coordinate system specified by itype.
    # ===================================================
    one   = x
    x     = x*cd +   z*sd
    z     = z*cd - one*sd
    f     = sqrt(x^2 + y^2 + z^2)

    return x,y,z,f
end
