################################################################################
# Implementation of the PolarizedHealpixMap type

"""
    GenericPolarizedHealpixMap{T}

An abstract type representing an Healpix polarized map without a
specified ordering. This can be used to implement multiple dispatch
when you don't care about the ordering of a map.
"""
abstract type GenericPolarizedHealpixMap{T} end

@doc raw"""
    mutable struct PolarizedHealpixMap{T, O <: Healpix.Order, AA <: AbstractArray{T, 1}}

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
mutable struct PolarizedHealpixMap{T,O<:Order,AA<:AbstractArray{T,1}} <: GenericPolarizedHealpixMap{T}
    i::HealpixMap{T,O,AA}
    q::HealpixMap{T,O,AA}
    u::HealpixMap{T,O,AA}

    function PolarizedHealpixMap{T,O,AA}(
        i::HealpixMap{T,O,AA},
        q::HealpixMap{T,O,AA},
        u::HealpixMap{T,O,AA},
    ) where {T,O<:Order,AA<:AbstractArray{T,1}}
        ((length(i) != length(q)) || (length(i) != length(q))) &&
            throw(ArgumentError("The three I/Q/U maps must have the same resolution"),)

        new(i, q, u)
    end

    function PolarizedHealpixMap{T,O,AA}(
        i::AA,
        q::AA,
        u::AA,
    ) where {T,O<:Order,AA<:AbstractArray{T,1}}
        ((length(i) != length(q)) || (length(i) != length(q))) &&
            throw(ArgumentError("The three I/Q/U vectors must have the same resolution"),)

        new(HealpixMap{T,O,AA}(i), HealpixMap{T,O,AA}(q), HealpixMap{T,O,AA}(u))
    end

    function PolarizedHealpixMap{T,O,AA}(nside::Number) where {T,O<:Order,AA<:AbstractArray{T,1}}
        new(HealpixMap{T,O,AA}(nside), HealpixMap{T,O,AA}(nside), HealpixMap{T,O,AA}(nside))
    end
end

function PolarizedHealpixMap{T,O}(nside::Number) where {T,O<:Order}
    PolarizedHealpixMap{T,O,Array{T,1}}(nside)
end

function PolarizedHealpixMap{T,O}(i::Array{T,1}, q::Array{T,1}, u::Array{T,1}) where {T,O<:Order}
    PolarizedHealpixMap{T,O,Array{T,1}}(i, q, u)
end

function PolarizedHealpixMap(
    i::HealpixMap{T,O,AA},
    q::HealpixMap{T,O,AA},
    u::HealpixMap{T,O,AA},
) where {T,O<:Order,AA<:AbstractArray{T,1}}
    return PolarizedHealpixMap{T,O,AA}(i, q, u)
end
