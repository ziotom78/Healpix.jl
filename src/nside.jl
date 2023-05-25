@doc raw"""
    ORDER_MAX

Maximum order for the resolution of a map supported on this machine.
The value of `NSIDE_MAX` is equal to `2^ORDER_MAX`.
"""
const ORDER_MAX = floor(Int, 0.5 * log2(typemax(Int) / 12))

@doc raw"""
    NSIDE_MAX

Maximum allowed value for the NSIDE resolution parameter.

"""
const NSIDE_MAX = 2^ORDER_MAX

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

"""
    nside2order(nside::Integer)

Return the order (positive integer) associated with a given NSIDE. If
the given nside is not valid, throw a `DomainError` exception.

If you have created a [`Healpix.Resolution`](@ref) object, you can
access the order through the field `order`.

See also [`order2nside`](@ref).
"""
function nside2order(nside::Integer)
    nsideok(nside) || throw(DomainError("Invalid value for NSIDE = $nside"))

    round(Int, log2(nside))
end


"""
    order2nside(order::Integer)

Return the value of NSIDE for a given order. If the given order is
not valid, throw a `DomainError` exception.

If you have created a [`Healpix.Resolution`](@ref) object, you can
access the value of NSIDE through the field `nside`.

See also [`nside2order`](@ref).
"""
function order2nside(order::Integer)
    ((order >= 0) && (order <= ORDER_MAX)) || throw(DomainError("Invalid order = $order"))

    2^order
end
