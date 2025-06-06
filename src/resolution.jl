# Definition of the composite type "Resolution"

using Printf

"""
    struct Resolution

`Resolution` objects are needed to perform a number of pixel-related
functions, e.g., convert a direction into a pixel number and vice
versa.

The fields of a `Resolution` object are the following:

- `nside`: the NSIDE parameter
- `nsideTimesTwo`: 2 * NSIDE
- `nsideTimesFour`: 4 * NSIDE
- `numOfPixels`: number of pixels in the map
- `order`: order of the map
- `pixelsPerFace`: number of pixels in each Healpix face
- `ncap`
- `fact2`
- `fact1`
"""
struct Resolution
    nside::Int
    nsideTimesTwo::Int
    nsideTimesFour::Int
    numOfPixels::Int

    order::Int
    pixelsPerFace::Int
    ncap::Int
    fact2::Float64
    fact1::Float64
end

# Customize printing
Base.show(io::IO, r::Resolution) = @printf(io, "Healpix resolution(NSIDE = %d)", r.nside)

################################################################################

"""
    Resolution(nside) -> Resolution

Create a `Resolution` object, given a value for `NSIDE`.
"""
function Resolution(nside)
    (1 ≤ nside ≤ NSIDE_MAX) || throw(DomainError(nside))
    # The expression (nside & (nside - 1)) == 0 is a quick check for
    # detecting if nside is a power of two or not
    (nside & (nside - 1) == 0) || throw(DomainError())

    order = ilog2(nside)
    pixelsPerFace = nside^2
    numOfPixels = 12pixelsPerFace
    ncap = 2 * nside * (nside - 1)
    fact2 = 4 / numOfPixels
    fact1 = 2 * nside * fact2

    result = Resolution(
        nside,
        2nside,
        4nside,
        numOfPixels,
        order,
        pixelsPerFace,
        ncap,
        fact2,
        fact1,
    )

end

"""
    numOfRings(resol::Resolution)
    numOfRings(nside::Integer)

Return the number of horizontal rings in a Healpix map.
"""
numOfRings(resol::Resolution) = resol.nsideTimesFour - 1
numOfRings(nside::Integer) = 4*nside - 1
