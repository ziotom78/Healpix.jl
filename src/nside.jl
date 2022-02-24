@doc raw"""
    NSIDE_MAX

Maximum allowed value for the NSIDE resolution parameter.

"""
const NSIDE_MAX = 2^(floor(Int, 0.5 * ((sizeof(Int) * 8 - 1) - log2(12))))

########################################################################

"""
    nsideok(nside::Integer) -> Bool

Check whether `nside` is a valid `NSIDE` parameter.
"""
nsideok(nside::Integer) = (nside > 0) && ((nside) & (nside - 1) == 0)

########################################################################

"""
    nside2npix(nside::Integer) -> Integer

Return the number of pixels for a Healpix map with the specified
`NSIDE` value. If `NSIDE` is not an integer power of two, the function
throws a `DomainError` exception.
"""
function nside2npix(nside::Integer)
    (nside > 0) || throw(DomainError(nside, "`NSIDE` is not a positive number"))

    nsidelog2 = round(Int, log2(nside))
    (2^nsidelog2 == nside) ||
        throw(DomainError(nside, "`NSIDE` is not an integer power of two"))

    12(nside^2)
end

########################################################################

"""
    npix2nside(npix::Integer) -> Integer

Given the number of pixels in a Healpix map, return the `NSIDE`
resolution parameter. If the number is invalid, throw a `DomainError`
exception.
"""
function npix2nside(npix::Integer)
    (npix % 12 == 0) || throw(DomainError(npix, "Invalid number of pixels"))

    square_root = sqrt(npix / 12)
    (square_root^2 == npix / 12) || throw(DomainError(npix, "Invalid number of pixels"))

    convert(Int, round(square_root))
end

################################################################################

"""
    nside2pixarea(nside::Integer) -> Real

Return the solid angle of a pixel in a map with the specified `NSIDE` parameter.
The result is expressed in steradians.
"""
nside2pixarea(nside::Integer) = 4π / nside2npix(nside)

################################################################################

"""
    nside2resol(nside::Integer) -> Real

Return the approximate resolution of a map with the specified `NSIDE`. The
resolution is expressed in radians, and it is the square root of the pixel
size.
"""
nside2resol(nside::Integer) = sqrt(nside2pixarea(nside))

################################################################################
