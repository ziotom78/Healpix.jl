# Definition of the composite type "Resolution"

"""
`Resolution` objects are needed to perform a number of
pixel-related functions, e.g., convert a direction into a pixel number
and vice versa.
"""
struct Resolution
    nside
    nsideTimesTwo
    nsideTimesFour
    numOfPixels

    order
    pixelsPerFace
    ncap
    fact2
    fact1
end

################################################################################

"""
    Resolution(nside) -> Resolution

Create a `Resolution` object, given a value for `NSIDE`.
"""
function Resolution(nside)
    (1 ≤ nside ≤ NSIDE_MAX) || throw(DomainError())
    # The expression (nside & (nside - 1)) == 0 is a quick check for
    # detecting if nside is a power of two or not
    (nside & (nside - 1) == 0) || throw(DomainError())
    
    order          = ilog2(nside)
    pixelsPerFace  = nside^2
    numOfPixels    = 12pixelsPerFace
    ncap           = 2 * nside * (nside - 1)
    fact2          = 4 / numOfPixels
    fact1          = 2 * nside * fact2

    result = Resolution(nside,
                        2nside,
                        4nside,
                        numOfPixels,
                        order,
                        pixelsPerFace,
                        ncap,
                        fact2,
                        fact1)

end
