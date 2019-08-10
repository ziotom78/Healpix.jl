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
    mutable struct PolarizedMap{T, O <: Healpix.Order}

A polarized I/Q/U map. It contains three Healpix maps with the same NSIDE:

- `i`
- `q`
- `u`

You can create an instance of this type using the function
[`PolarizedMap{T,O}`](@ref), which comes in three flavours:

- `PolarizedMap(i::Map{T,O}, q::Map{T,O}, u::Map{T,O})`
- `PolarizedMap{T,O}(i::AbstractVector{T}, q::AbstractVector{T}, u::AbstractVector{T})`
- `PolarizedMap{T,O}(nside::Number)`

"""
mutable struct PolarizedMap{T, O <: Order} <: GenericPolarizedMap{T}
    i::Map{T,O}
    q::Map{T,O}
    u::Map{T,O}

    function PolarizedMap{T, O}(
        i::Map{T, O},
        q::Map{T, O},
        u::Map{T, O},
    ) where {T, O <: Order}
        ((length(i) != length(q)) || (length(i) != length(q))) && throw(
            ArgumentError("The three I/Q/U maps must have the same resolution"),
        )

        new(i, q, u)
    end

    function PolarizedMap{T, O}(
        i::AbstractVector{T},
        q::AbstractVector{T},
        u::AbstractVector{T},
    ) where {T, O <: Order}

        PolarizedMap{T, O}(
            Map{T, O}(i),
            Map{T, O}(q),
            Map{T, O}(u),
        )
    end

    function PolarizedMap{T, O}(nside::Number) where {T, O <: Order}
        PolarizedMap{T, O}(
            Map{T, O}(nside),
            Map{T, O}(nside),
            Map{T, O}(nside),
        )
    end
end
