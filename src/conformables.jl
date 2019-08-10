################################################################################

conformables(map1::Map{T, RingOrder}, map2::Map{S, RingOrder}) where {T, S} =
    map1.resolution.nside == map2.resolution.nside

conformables(map1::Map{T, NestedOrder}, map2::Map{S, NestedOrder}) where {T, S} =
    map1.resolution.nside == map2.resolution.nside

conformables(map1::Map{T, O1},
             map2::Map{S, O2}) where {T, S, O1 <: Order, O2 <: Order} = false

################################################################################

conformables(map1::PolarizedMap{T, RingOrder}, map2::PolarizedMap{S, RingOrder}) where {T, S} =
    map1.i.resolution.nside == map2.i.resolution.nside

conformables(map1::PolarizedMap{T, NestedOrder}, map2::PolarizedMap{S, NestedOrder}) where {T, S} =
    map1.i.resolution.nside == map2.i.resolution.nside

conformables(map1::PolarizedMap{T, O1},
             map2::PolarizedMap{S, O2}) where {T, S, O1 <: Order, O2 <: Order} = false

"""
    conformables{T, S, O1 <: Order, O2 <: Order}(map1::Map{T, O1},
                                                 map2::Map{S, O2}) -> Bool
    conformables{T, S, O1 <: Order, O2 <: Order}(map1::PolarizedMap{T, O1},
                                                 map2::PolarizedMap{S, O2}) -> Bool

Determine if two Healpix maps are "conformables", i.e., if their
shape and ordering are the same.
"""
conformables
