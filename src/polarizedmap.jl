################################################################################
# Implementation of the PolarizedMap type

"""
    GenericPolarizedMap{T}

An abstract type representing an Healpix polarized map without a
specified ordering. This can be used to implement multiple dispatch
when you don't care about the ordering of a map.
"""
abstract type GenericPolarizedMap{T} end

@doc raw"""
    mutable struct PolarizedMap{T, O <: Healpix.Order, AA <: AbstractArray{T, 1}}

A polarized I/Q/U map. It contains three Healpix maps with the same NSIDE:

- `i`
- `q`
- `u`

You can create an instance of this type using the function
[`PolarizedMap{T,O}`](@ref), which comes in three flavours:

- `PolarizedMap(i::Map{T,O,AA}, q::Map{T,O,AA}, u::Map{T,O,AA})`
- `PolarizedMap{T,O}(i::AbstractVector{T}, q::AbstractVector{T}, u::AbstractVector{T})`
- `PolarizedMap{T,O}(nside::Number)`

"""
mutable struct PolarizedMap{T, O <: Order, AA <: AbstractArray{T, 1}} <: GenericPolarizedMap{T}
    i::Map{T,O,AA}
    q::Map{T,O,AA}
    u::Map{T,O,AA}

    function PolarizedMap{T, O, AA}(
        i::Map{T, O, AA},
        q::Map{T, O, AA},
        u::Map{T, O, AA},
    ) where {T, O <: Order, AA <: AbstractArray{T, 1}}
        ((length(i) != length(q)) || (length(i) != length(q))) && throw(
            ArgumentError("The three I/Q/U maps must have the same resolution"),
        )

        new(i, q, u)
    end

    function PolarizedMap{T, O, AA}(
        i::AA,
        q::AA,
        u::AA,
    ) where {T, O <: Order, AA <: AbstractArray{T, 1}}
        ((length(i) != length(q)) || (length(i) != length(q))) && throw(
            ArgumentError("The three I/Q/U vectors must have the same resolution"),
        )

        new(
            Map{T, O, AA}(i),
            Map{T, O, AA}(q),
            Map{T, O, AA}(u),
        )
    end

    function PolarizedMap{T, O, AA}(
        nside::Number
    ) where {T, O <: Order, AA <: AbstractArray{T, 1}}
        new(
            Map{T, O, AA}(nside),
            Map{T, O, AA}(nside),
            Map{T, O, AA}(nside),
        )
    end
end

function PolarizedMap{T, O}(nside::Number) where {T, O <: Order}
    PolarizedMap{T, O, Array{T, 1}}(nside)
end

function PolarizedMap{T, O}(
    i::Array{T, 1},
    q::Array{T, 1},
    u::Array{T, 1},
) where {T, O <: Order}
    PolarizedMap{T, O, Array{T, 1}}(i, q, u)
end
