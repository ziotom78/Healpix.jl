################################################################################

ang2pix(map::Map{T,RingOrder,AA}, theta, phi) where {T,AA} =
    ang2pixRing(map.resolution, theta, phi)
ang2pix(map::Map{T,NestedOrder,AA}, theta, phi) where {T,AA} =
    ang2pixNest(map.resolution, theta, phi)
ang2pix(map::PolarizedMap{T,RingOrder,AA}, theta, phi) where {T,AA} =
    ang2pixRing(map.i.resolution, theta, phi)
ang2pix(map::PolarizedMap{T,NestedOrder,AA}, theta, phi) where {T,AA} =
    ang2pixNest(map.i.resolution, theta, phi)

@doc raw"""
    ang2pix{T, O, AA}(map::Map{T, O}, theta, phi)
    ang2pix{T, O, AA}(map::PolarizedMap{T, O}, theta, phi)

Convert the direction specified by the colatitude `theta` (∈ [0, π])
and the longitude `phi` (∈ [0, 2π]) into the index of the pixel in the
Healpix map `map`.
"""
ang2pix

################################################################################

pix2ang(map::Map{T,RingOrder,AA}, ipix) where {T,AA} = pix2angRing(map.resolution, ipix)
pix2ang(map::Map{T,NestedOrder,AA}, ipix) where {T,AA} = pix2angNest(map.resolution, ipix)
pix2ang(map::PolarizedMap{T,RingOrder,AA}, ipix) where {T,AA} =
    pix2angRing(map.i.resolution, ipix)
pix2ang(map::PolarizedMap{T,NestedOrder,AA}, ipix) where {T,AA} =
    pix2angNest(map.i.resolution, ipix)

@doc raw"""
    pix2ang{T, O <: Order}(map::Map{T, O}, ipix) -> (Float64, Float64)
    pix2ang{T, O <: Order}(map::PolarizedMap{T, O}, ipix) -> (Float64, Float64)

Return the pair (`theta`, `phi`), where `theta` is the colatitude and
`phi` the longitude of the direction of the pixel center with index
`ipix`.
"""
pix2ang

################################################################################
# Interpolation

function interpolate(m::Map{T,RingOrder,AA}, θ, ϕ, pixbuf, weightbuf) where {T,AA}
    getinterpolRing(m.resolution, θ, ϕ, pixbuf, weightbuf)

    result = zero(weightbuf[1])
    for i = 1:4
        result += m[pixbuf[i]] * weightbuf[i]
    end

    result
end

function interpolate(m::Map{T,RingOrder,AA}, θ, ϕ) where {T,AA}
    pixbuf = Array{Int}(undef, 4)
    weightbuf = Array{Float64}(undef, 4)

    interpolate(m, θ, ϕ, pixbuf, weightbuf)
end

"""
    interpolate(m::Map{T, RingOrder, AA}, θ, ϕ) -> Value
    interpolate(m::Map{T, RingOrder, AA}, θ, ϕ, pixbuf, weightbuf) -> Value

Return an interpolated value of the map along the specified direction.

When provided, the parameters `pixbuf` and `weightbuf` must be
4-element arrays of integer and floating-point values,
respectively. They can be reused across multiple calls of
`interpolate!`, to save heap allocations:

```
pixbuf = Array{Int}(undef, 4)
weightbuf = Array{Float64}(undef, 4)

m = Map{Float64, RingOrder}(1)
for (θ, ϕ) in [(0., 0.), (π/2, π/2)]
    println(interpolate!(m, θ, ϕ, pixbuf, weightbuf))
end
```
"""
interpolate
