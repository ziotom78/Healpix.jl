################################################################################

ang2pix(map::Map{T, RingOrder}, theta, phi) where {T} = ang2pixRing(map.resolution, theta, phi)
ang2pix(map::Map{T, NestedOrder}, theta, phi) where {T} = ang2pixNest(map.resolution, theta, phi)
ang2pix(map::PolarizedMap{T, RingOrder}, theta, phi) where {T} =
    ang2pixRing(map.i.resolution, theta, phi)
ang2pix(map::PolarizedMap{T, NestedOrder}, theta, phi) where {T} =
    ang2pixNest(map.i.resolution, theta, phi)

@doc raw"""
    ang2pix{T, O <: Order}(map::Map{T, O}, theta, phi)
    ang2pix{T, O <: Order}(map::PolarizedMap{T, O}, theta, phi)

Convert the direction specified by the colatitude `theta` (∈ [0, π])
and the longitude `phi` (∈ [0, 2π]) into the index of the pixel in the
Healpix map `map`.
"""
ang2pix

################################################################################

pix2ang(map::Map{T, RingOrder}, ipix) where {T} = pix2angRing(map.resolution, ipix)
pix2ang(map::Map{T, NestedOrder}, ipix) where {T} = pix2angNest(map.resolution, ipix)
pix2ang(map::PolarizedMap{T, RingOrder}, ipix) where {T} =
    pix2angRing(map.i.resolution, ipix)
pix2ang(map::PolarizedMap{T, NestedOrder}, ipix) where {T} =
    pix2angNest(map.i.resolution, ipix)

@doc raw"""
    pix2ang{T, O <: Order}(map::Map{T, O}, ipix) -> (Float64, Float64)
    pix2ang{T, O <: Order}(map::PolarizedMap{T, O}, ipix) -> (Float64, Float64)

Return the pair (`theta`, `phi`), where `theta` is the colatitude and
`phi` the longitude of the direction of the pixel center with index
`ipix`.
"""
pix2ang
