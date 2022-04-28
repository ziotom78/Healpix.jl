################################################################################

conformables(map1::HealpixMap{T,RingOrder,AA1}, map2::HealpixMap{S,RingOrder,AA2}) where {T,S,AA1,AA2} =
    map1.resolution.nside == map2.resolution.nside

conformables(
    map1::HealpixMap{T,NestedOrder,AA1},
    map2::HealpixMap{S,NestedOrder,AA2},
) where {T,S,AA1,AA2} = map1.resolution.nside == map2.resolution.nside

conformables(map1::HealpixMap{T,O1,AA1}, map2::HealpixMap{S,O2,AA2}) where {T,S,O1,O2,AA1,AA2} = false

################################################################################

conformables(
    map1::PolarizedHealpixMap{T,RingOrder,AA1},
    map2::PolarizedHealpixMap{S,RingOrder,AA2},
) where {T,S,AA1,AA2} = map1.i.resolution.nside == map2.i.resolution.nside

conformables(
    map1::PolarizedHealpixMap{T,NestedOrder,AA1},
    map2::PolarizedHealpixMap{S,NestedOrder,AA2},
) where {T,S,AA1,AA2} = map1.i.resolution.nside == map2.i.resolution.nside

conformables(
    map1::PolarizedHealpixMap{T,O1,AA1},
    map2::PolarizedHealpixMap{S,O2,AA2},
) where {T,S,O1,O2,AA1,AA2} = false

"""
    conformables{T, S, O1, O2}(map1::HealpixMap{T, O1, AA1},
                               map2::HealpixMap{S, O2, AA2}) -> Bool
    conformables{T, S, O1, O2}(map1::PolarizedHealpixMap{T, O1, AA1},
                               map2::PolarizedHealpixMap{S, O2, AA2}) -> Bool

Determine if two Healpix maps are "conformables", i.e., if their shape
and ordering are the same. The array types `AA1` and `AA2` are not considered
in testing conformability.
"""
conformables
