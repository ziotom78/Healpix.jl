################################################################################

conformables(map1::Map{T,RingOrder,AA1}, map2::Map{S,RingOrder,AA2}) where {T,S,AA1,AA2} =
    map1.resolution.nside == map2.resolution.nside

conformables(
    map1::Map{T,NestedOrder,AA1},
    map2::Map{S,NestedOrder,AA2},
) where {T,S,AA1,AA2} = map1.resolution.nside == map2.resolution.nside

conformables(map1::Map{T,O1,AA1}, map2::Map{S,O2,AA2}) where {T,S,O1,O2,AA1,AA2} = false

################################################################################

conformables(
    map1::PolarizedMap{T,RingOrder,AA1},
    map2::PolarizedMap{S,RingOrder,AA2},
) where {T,S,AA1,AA2} = map1.i.resolution.nside == map2.i.resolution.nside

conformables(
    map1::PolarizedMap{T,NestedOrder,AA1},
    map2::PolarizedMap{S,NestedOrder,AA2},
) where {T,S,AA1,AA2} = map1.i.resolution.nside == map2.i.resolution.nside

conformables(
    map1::PolarizedMap{T,O1,AA1},
    map2::PolarizedMap{S,O2,AA2},
) where {T,S,O1,O2,AA1,AA2} = false

"""
    conformables{T, S, O1, O2}(map1::Map{T, O1, AA1},
                               map2::Map{S, O2, AA2}) -> Bool
    conformables{T, S, O1, O2}(map1::PolarizedMap{T, O1, AA1},
                               map2::PolarizedMap{S, O2, AA2}) -> Bool

Determine if two Healpix maps are "conformables", i.e., if their shape
and ordering are the same. The array types `AA1` and `AA2` are not considered
in testing conformability.
"""
conformables
