################################################################################
# Implementation of the Order and *HealpixMap types

"""Abstract type representing the ordering of pixels in a Healpix map.
See also `RingOrder` and `NestedOrder`.
"""
abstract type Order end

"""The `RingOrder` type should be used when creating `HealpixMap` types in
order to specify that the pixels in the map are sorted in ring
ordering. (See also `NestedOrder`.)
"""
abstract type RingOrder <: Order end

"""The `NestedOrder` type should be used when creating `HealpixMap` types in
order to specify that the pixels in the map are sorted in ring
ordering. (See also `RingOrder`.)
"""
abstract type NestedOrder <: Order end

"""
    abstract type AbstractHealpixMap{T} <: AbstractVector{T}

An abstract type representing an Healpix map without a specified
ordering. This can be used to implement multiple dispatch when you
don't care about the ordering of a map.
"""
abstract type AbstractHealpixMap{T} <: AbstractVector{T} end

"""
    HealpixMap{T, O <: Order, AA <: AbstractVector{T}} <: AbstractHealpixMap{T}

A Healpix map. The type `T` is used for the value of the pixels in a
map, and it can be anything (even a string!). The type `O` is used to
specify the ordering of the pixels, and it can either be `RingOrder`
or `NestedOrder`. The type `AA` is used to store the array of pixels;
typical types are `Vector`, `CUArray`, `SharedArray`, etc.

A `HealpixMap` type contains the following fields:

- `pixels`: array of pixels
- `resolution`: instance of a `Resolution` object

You can construct a map using one of the following forms:

- `HealpixMap{T, O, AA}(arr)` and `HealpixMap{T, O, AA}(nside::Number)` will use
  `AA` as basetype

- `HealpixMap{T, O}(arr)` and `HealpixMap{T, O}(nside::Number)` will use `Vector{T}`
  as basetype

# Examples

The following example creates a map with `NSIDE=32` in `RING` order,
containing integer values starting from 1:

    mymap = Healpix.HealpixMap{Int64, Healpix.RingOrder}(1:Healpix.nside2npix(32))

The call to `collect` is required to convert the range in an array.

This example creates a map in `NESTED` order, with `NSIDE=64`, filled
with zeroes:

    mymap = Healpix.HealpixMap{Float64, Healpix.NestedOrder}(64)

Finally, the following examples show how to use `SharedArray`:

    using SharedArrays
    
    # Create a map with all pixels set to zero
    mymap = Healpix.HealpixMap{Float64, Healpix.NestedOrder, SharedArray{Float64, 1}}(64)

    # Create a map and initialize pixel values with a SharedArray
    pixels = SharedArray{Int64, 1}(1:12 |> collect)
    mymap = Healpix.HealpixMap{Int64, Healpix.RingOrder, SharedArray{Int64, 1}}(m)
"""
mutable struct HealpixMap{T,O<:Order,AA<:AbstractVector{T}} <: AbstractHealpixMap{T}
    pixels::AA
    resolution::Resolution
    
    """
        HealpixMap{Union{T, Nothing}, O <: Order}(nside) -> HealpixMap{Union{T, Nothing}, O}

    Create an empty map with the specified NSIDE. All the pixels in the map
    are set to `nothing`.
    """
    function HealpixMap{Union{T, Nothing},O,AA}(nside::Number) where {T, O <: Order, AA <: AbstractVector{Union{T, Nothing}}}
        new(Union{T, Nothing}[nothing for i in 1:nside2npix(nside)], Resolution(nside))
    end

    """
        HealpixMap{T, O <: Order}(nside) -> HealpixMap{T, O}

    Create an empty map with the specified NSIDE. All the pixels in the map
    are set to zero.
    """
    function HealpixMap{T,O,AA}(nside::Number) where {T,O<:Order,AA<:AbstractVector{T}}
        new(zeros(T, nside2npix(nside)), Resolution(nside))
    end

    """
    Initialize a map from a generic array
    """
    function HealpixMap{T,O,AA}(arr::AA) where {T,O<:Order,AA<:AbstractVector{T}}
        nside = npix2nside(length(arr))
        new(arr, Resolution(nside))
    end
end

"""
    HealpixMap{T,O}(nside::Number) where {T,O<:Order}

Convenience function that uses `Vector{T}` as the type used to
hold the vector of pixels.
"""
function HealpixMap{T,O}(nside::Number) where {T,O<:Order}
    HealpixMap{T,O,Vector{T}}(nside)
end

"""
    HealpixMap{T,O}(arr) where {T,O<:Order}

Convenience function that uses `Vector{T}` as the type used to
hold the vector of pixels.
"""
function HealpixMap{T,O}(arr) where {T,O<:Order}
    HealpixMap{T,O,Vector{T}}(collect(arr))
end

import Base: +, -, *, /

+(a::HealpixMap{T,O,AA}, b::HealpixMap{T,O,AA}) where {T<:Number,O,AA} = HealpixMap{T,O,AA}(a.pixels .+ b.pixels)
-(a::HealpixMap{T,O,AA}, b::HealpixMap{T,O,AA}) where {T<:Number,O,AA} = HealpixMap{T,O,AA}(a.pixels .- b.pixels)
*(a::HealpixMap{T,O,AA}, b::HealpixMap{T,O,AA}) where {T<:Number,O,AA} = HealpixMap{T,O,AA}(a.pixels .* b.pixels)
/(a::HealpixMap{T,O,AA}, b::HealpixMap{T,O,AA}) where {T<:Number,O,AA} = HealpixMap{T,O,AA}(a.pixels ./ b.pixels)

+(a::HealpixMap{T,O,AA}, b::Number) where {T<:Number,O,AA} = HealpixMap{T,O,AA}(a.pixels .+ b)
-(a::HealpixMap{T,O,AA}, b::Number) where {T<:Number,O,AA} = a + (-b)
*(a::HealpixMap{T,O,AA}, b::Number) where {T<:Number,O,AA} = HealpixMap{T,O,AA}(a.pixels .* b)
/(a::HealpixMap{T,O,AA}, b::Number) where {T<:Number,O,AA} = HealpixMap{T,O,AA}(a.pixels ./ b)

+(a::Number, b::HealpixMap{T,O,AA}) where {T<:Number,O,AA} = b + a
-(a::Number, b::HealpixMap{T,O,AA}) where {T<:Number,O,AA} = b + (-a)
*(a::Number, b::HealpixMap{T,O,AA}) where {T<:Number,O,AA} = b * a
/(a::Number, b::HealpixMap{T,O,AA}) where {T<:Number,O,AA} = HealpixMap{T,O,AA}(a ./ b.pixels)

################################################################################
# Iterator interface

Base.size(m::HealpixMap{T,O,AA}) where {T,O,AA} = (m.resolution.numOfPixels,)
Base.parent(m::HealpixMap) = m.pixels
Base.IndexStyle(::Type{<:HealpixMap{T,O,AA}}) where {T,O,AA} = IndexLinear()

function getindex(m::HealpixMap{T,O,AA}, i::Integer) where {T,O,AA}
    1 ≤ i ≤ m.resolution.numOfPixels || throw(BoundsError(m, i))
    m.pixels[i]
end

function setindex!(m::HealpixMap{T,O,AA}, val, i::Integer) where {T,O,AA}
    1 ≤ i ≤ m.resolution.numOfPixels || throw(BoundsError(m, i))
    m.pixels[i] = val
end
