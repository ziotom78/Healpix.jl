################################################################################
# Implementation of the PolarizedHealpixMap type

"""
    AbstractPolarizedHealpixMap{T}

An abstract type representing an Healpix polarized map without a
specified ordering. This can be used to implement multiple dispatch
when you don't care about the ordering of a map.
"""
abstract type AbstractPolarizedHealpixMap{T} end

"""
    mutable struct PolarizedHealpixMap{T, O <: Healpix.Order, AA <: AbstractVector{T}}

A polarized I/Q/U map. It contains three Healpix maps with the same NSIDE:

- `i`
- `q`
- `u`

You can create an instance of this type using the function
[`PolarizedHealpixMap{T,O}`](@ref), which comes in three flavours:

- `PolarizedHealpixMap(i::HealpixMap{T,O,AA}, q::HealpixMap{T,O,AA}, u::HealpixMap{T,O,AA})`
- `PolarizedHealpixMap{T,O}(i::AbstractVector{T}, q::AbstractVector{T}, u::AbstractVector{T})`
- `PolarizedHealpixMap{T,O}(nside::Number)`

"""
mutable struct PolarizedHealpixMap{T,O<:Order,AA<:AbstractVector{T}} <: AbstractPolarizedHealpixMap{T}
    i::HealpixMap{T,O,AA}
    q::HealpixMap{T,O,AA}
    u::HealpixMap{T,O,AA}

    function PolarizedHealpixMap{T,O,AA}(
        i::HealpixMap{T,O,AA},
        q::HealpixMap{T,O,AA},
        u::HealpixMap{T,O,AA},
    ) where {T,O<:Order,AA<:AbstractVector{T}}
        ((length(i) != length(q)) || (length(i) != length(q))) &&
            throw(ArgumentError("The three I/Q/U maps must have the same resolution"),)

        new(i, q, u)
    end

    function PolarizedHealpixMap{T,O,AA}(
        i::AA,
        q::AA,
        u::AA,
    ) where {T,O<:Order,AA<:AbstractVector{T}}
        ((length(i) != length(q)) || (length(i) != length(q))) &&
            throw(ArgumentError("The three I/Q/U vectors must have the same resolution"),)

        new(HealpixMap{T,O,AA}(i), HealpixMap{T,O,AA}(q), HealpixMap{T,O,AA}(u))
    end

    function PolarizedHealpixMap{T,O,AA}(nside::Number) where {T,O<:Order,AA<:AbstractVector{T}}
        new(HealpixMap{T,O,AA}(nside), HealpixMap{T,O,AA}(nside), HealpixMap{T,O,AA}(nside))
    end
end

function PolarizedHealpixMap{T,O}(nside::Number) where {T,O<:Order}
    PolarizedHealpixMap{T,O,Vector{T}}(nside)
end

function PolarizedHealpixMap{T,O}(i::Vector{T}, q::Vector{T}, u::Vector{T}) where {T,O<:Order}
    PolarizedHealpixMap{T,O,Vector{T}}(i, q, u)
end

function PolarizedHealpixMap(
    i::HealpixMap{T,O,AA},
    q::HealpixMap{T,O,AA},
    u::HealpixMap{T,O,AA},
) where {T,O<:Order,AA<:AbstractVector{T}}
    return PolarizedHealpixMap{T,O,AA}(i, q, u)
end
