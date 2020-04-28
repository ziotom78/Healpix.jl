################################################################################
# Implementation of the Order and Map types

"""Abstract type representing the ordering of pixels in a Healpix map.
See also `RingOrder` and `NestedOrder`.
"""
abstract type Order end

"""The `RingOrder` type should be used when creating `Map` types in
order to specify that the pixels in the map are sorted in ring
ordering. (See also `NestedOrder`.)
"""
abstract type RingOrder <: Order end

"""The `NestedOrder` type should be used when creating `Map` types in
order to specify that the pixels in the map are sorted in ring
ordering. (See also `RingOrder`.)
"""
abstract type NestedOrder <: Order end

"""
    GenericMap{T} <: AbstractArray{T, 1}

An abstract type representing an Healpix map without a specified
ordering. This can be used to implement multiple dispatch when you
don't care about the ordering of a map."""
abstract type GenericMap{T} <: AbstractArray{T, 1} end

"""
    struct Map{T, O <: Order} <: GenericMap{T}

A Healpix map. The type `T` is used for the value of the pixels in
a map, and it can be anything (even a string!). The type `O` is used
to specify the ordering of the pixels, and it can either be
`RingOrder` or `NestedOrder`.

A `Map` type contains the following fields:

- `pixels`: array of pixels
- `resolution`: instance of a `Resolution` object

"""
mutable struct Map{T, O <: Order} <: GenericMap{T}
    pixels::Array{T,1}
    resolution::Resolution

    """
        Map{T, O <: Order}(nside) -> Map{T, O}

    Create an empty map with the specified NSIDE.
    """
    Map{T, O}(nside::Number) where {T, O <: Order} = new(zeros(T, nside2npix(nside)),
                                                 Resolution(nside))

    """
    Create a map with the specified array of pixels.
    """
    function Map{T, O}(arr::Array{T,1}) where {T, O <: Order}
        nside = npix2nside(length(arr))
        new(arr, Resolution(nside))
    end
end

import Base: +, -, *, /

+(a::Map{T,O}, b::Map{T,O}) where {T <: Number, O} = Map{T, O}(a.pixels .+ b.pixels)
-(a::Map{T,O}, b::Map{T,O}) where {T <: Number, O} = Map{T, O}(a.pixels .- b.pixels)
*(a::Map{T,O}, b::Map{T,O}) where {T <: Number, O} = Map{T, O}(a.pixels .* b.pixels)
/(a::Map{T,O}, b::Map{T,O}) where {T <: Number, O} = Map{T, O}(a.pixels ./ b.pixels)

+(a::Map{T,O}, b::Number) where {T <: Number, O} = Map{T, O}(a.pixels .+ b)
-(a::Map{T,O}, b::Number) where {T <: Number, O} = a + (-b)
*(a::Map{T,O}, b::Number) where {T <: Number, O} = Map{T, O}(a.pixels .* b)
/(a::Map{T,O}, b::Number) where {T <: Number, O} = Map{T, O}(a.pixels ./ b)

+(a::Number, b::Map{T,O}) where {T <: Number, O} = b + a
-(a::Number, b::Map{T,O}) where {T <: Number, O} = b + (-a)
*(a::Number, b::Map{T,O}) where {T <: Number, O} = b * a
/(a::Number, b::Map{T,O}) where {T <: Number, O} = Map{T, O}(a ./ b.pixels)

################################################################################
# Iterator interface

Base.size(m::Map{T, O}) where {T, O} = (m.resolution.numOfPixels,)

Base.IndexStyle(::Type{<:Map{T, O}}) where {T, O} = IndexLinear()

function getindex(m::Map{T, O}, i::Integer) where {T, O}
    1 ≤ i ≤ m.resolution.numOfPixels || throw(BoundsError(m, i))
    m.pixels[i]
end

function setindex!(m::Map{T, O}, val, i::Integer) where {T, O}
    1 ≤ i ≤ m.resolution.numOfPixels || throw(BoundsError(m, i))
    m.pixels[i] = val
end
