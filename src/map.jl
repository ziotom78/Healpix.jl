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
    AbstractHealpixMap{T} <: AbstractArray{T, 1}

An abstract type representing an Healpix map without a specified
ordering. This can be used to implement multiple dispatch when you
don't care about the ordering of a map."""
abstract type AbstractHealpixMap{T} <: AbstractArray{T,1} end

"""
    struct Map{T, O <: Order, AA <: AbstractArray{T, 1}} <: AbstractHealpixMap{T}

A Healpix map. The type `T` is used for the value of the pixels in a
map, and it can be anything (even a string!). The type `O` is used to
specify the ordering of the pixels, and it can either be `RingOrder`
or `NestedOrder`. The type `AA` is used to store the array of pixels;
typical types are `Vector`, `CUArray`, `SharedArray`, etc.

A `Map` type contains the following fields:

- `pixels`: array of pixels
- `resolution`: instance of a `Resolution` object

You can construct a map using one of the following forms:

- `Map{T, O, AA}(arr)` and `Map{T, O, AA}(nside::Number)` will use
  `AA` as basetype

- `Map{T, O}(arr)` and `Map{T, O}(nside::Number)` will use `Array{T,
  1}` as basetype

# Examples

The following example creates a map with `NSIDE=32` in `RING` order,
containing integer values starting from 1:

    mymap = Healpix.Map{Int64, Healpix.RingOrder}(1:Healpix.nside2npix(32))

The call to `collect` is required to convert the range in an array.

This example creates a map in `NESTED` order, with `NSIDE=64`, filled
with zeroes:

    mymap = Healpix.Map{Float64, Healpix.NestedOrder}(64)

Finally, the following examples show how to use `SharedArray`:

    using SharedArrays
    
    # Create a map with all pixels set to zero
    mymap = Healpix.Map{Float64, Healpix.NestedOrder, SharedArray{Float64, 1}}(64)

    # Create a map and initialize pixel values with a SharedArray
    pixels = SharedArray{Int64, 1}(1:12 |> collect)
    mymap = Healpix.Map{Int64, Healpix.RingOrder, SharedArray{Int64, 1}}(m)
"""
mutable struct Map{T,O<:Order,AA<:AbstractArray{T,1}} <: AbstractHealpixMap{T}
    pixels::AA
    resolution::Resolution

    """
        Map{T, O <: Order}(nside) -> Map{T, O}

    Create an empty map with the specified NSIDE.
    """
    function Map{T,O,AA}(nside::Number) where {T,O<:Order,AA<:AbstractArray{T,1}}
        new(zeros(T, nside2npix(nside)), Resolution(nside))
    end

    """
    Initialize a map from a generic array
    """
    function Map{T,O,AA}(arr::AA) where {T,O<:Order,AA<:AbstractArray{T,1}}
        nside = npix2nside(length(arr))
        new(arr, Resolution(nside))
    end
end

"""Convenience function that uses `Array{T, 1}` as the type used to
hold the vector of pixels.
"""
function Map{T,O}(nside::Number) where {T,O<:Order}
    Map{T,O,Array{T,1}}(nside)
end

"""Convenience function that uses `Array{T, 1}` as the type used to
hold the vector of pixels.
"""
function Map{T,O}(arr) where {T,O<:Order}
    Map{T,O,Array{T,1}}(collect(arr))
end

import Base: +, -, *, /

+(a::Map{T,O,AA}, b::Map{T,O,AA}) where {T<:Number,O,AA} = Map{T,O,AA}(a.pixels .+ b.pixels)
-(a::Map{T,O,AA}, b::Map{T,O,AA}) where {T<:Number,O,AA} = Map{T,O,AA}(a.pixels .- b.pixels)
*(a::Map{T,O,AA}, b::Map{T,O,AA}) where {T<:Number,O,AA} = Map{T,O,AA}(a.pixels .* b.pixels)
/(a::Map{T,O,AA}, b::Map{T,O,AA}) where {T<:Number,O,AA} = Map{T,O,AA}(a.pixels ./ b.pixels)

+(a::Map{T,O,AA}, b::Number) where {T<:Number,O,AA} = Map{T,O,AA}(a.pixels .+ b)
-(a::Map{T,O,AA}, b::Number) where {T<:Number,O,AA} = a + (-b)
*(a::Map{T,O,AA}, b::Number) where {T<:Number,O,AA} = Map{T,O,AA}(a.pixels .* b)
/(a::Map{T,O,AA}, b::Number) where {T<:Number,O,AA} = Map{T,O,AA}(a.pixels ./ b)

+(a::Number, b::Map{T,O,AA}) where {T<:Number,O,AA} = b + a
-(a::Number, b::Map{T,O,AA}) where {T<:Number,O,AA} = b + (-a)
*(a::Number, b::Map{T,O,AA}) where {T<:Number,O,AA} = b * a
/(a::Number, b::Map{T,O,AA}) where {T<:Number,O,AA} = Map{T,O,AA}(a ./ b.pixels)

################################################################################
# Iterator interface

Base.size(m::Map{T,O,AA}) where {T,O,AA} = (m.resolution.numOfPixels,)

Base.IndexStyle(::Type{<:Map{T,O,AA}}) where {T,O,AA} = IndexLinear()

function getindex(m::Map{T,O,AA}, i::Integer) where {T,O,AA}
    1 ≤ i ≤ m.resolution.numOfPixels || throw(BoundsError(m, i))
    m.pixels[i]
end

function setindex!(m::Map{T,O,AA}, val, i::Integer) where {T,O,AA}
    1 ≤ i ≤ m.resolution.numOfPixels || throw(BoundsError(m, i))
    m.pixels[i] = val
end
