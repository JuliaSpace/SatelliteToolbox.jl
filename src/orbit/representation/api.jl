# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Description
# ==============================================================================
#
#   API functions for all the orbit representations.
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

export get_epoch, get_a, get_e, get_i, get_Ω, get_ω, get_f, get_raan, get_argp,
       get_M

"""
    get_epoch(orb)

Return the epoch of the orbit representation `orb` [m].

"""
get_epoch(orb)

"""
    get_a(orb)

Return the semi-major axis of the orbit representation `orb` [m].

"""
get_a

"""
    get_e(orb)

Return the eccentricity of the orbit representation `orb`.

"""
get_e


"""
    get_i(orb)

Return the inclination of the orbit representation `orb` [rad].

"""
get_i

"""
    get_Ω(orb)

Return the right ascention of the ascending node (RAAN) of the orbit
representation `orb` [rad].

"""
get_Ω

"""
    get_ω(orb)

Return the argument of periapsis of the representation `orb` [rad].

"""
get_ω

"""
    get_f(orb)

Return the true anomaly of the representation `orb` [rad].

"""
get_f

"""
    get_raan(orb)

Return the right ascention of the ascending node (RAAN) of the orbit
representation `orb` [rad].

"""
@inline get_raan(orb) = get_Ω(orb)

"""
    get_argp(orb)

Return the argument of periapsis of the representation `orb` [rad].

"""
@inline get_argp(orb) = get_ω(orb)

"""
    get_M(orb)

Return the mean anomaly of the representation `orb` [rad].

"""
get_M

"""
    get_r(orb)

Return the position vector of the representation `orb` [m].

"""
get_r

"""
    get_r(orb)

Return the velocity vector of the representation `orb` [m].

"""
get_v

"""
    get_rv(orb)

Return the position and velocity vector of the representation `orb` [m].

"""
get_rv
