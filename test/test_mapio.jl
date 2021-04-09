# Map loading

m = Healpix.readMapFromFITS("float_map.fits", 1, Float32)
@test typeof(m) == Healpix.Map{Float32,Healpix.RingOrder,Array{Float32,1}}
@test m.resolution.nside == 4
@test m.pixels == [Float32(x) for x = 0:(12*4^2-1)]

m = Healpix.readMapFromFITS("int_map.fits", 1, Int8)
@test typeof(m) == Healpix.Map{Int8,Healpix.RingOrder,Array{Int8,1}}
@test m.resolution.nside == 1
@test m.pixels == [Int8(x) for x = 0:11]

# Map saving

const mapFileName = tempname()
print("Saving $mapFileName\n")
Healpix.saveToFITS(m, "!$mapFileName", typechar = "I")
m2 = Healpix.readMapFromFITS(mapFileName, 1, Int8)
@test m.pixels == m2.pixels


# Map reordering
m_ring = Healpix.readMapFromFITS("float_map.fits", 1, Float32)
m_nest = Healpix.readMapFromFITS("float_map_nest.fits", 1, Float32)
@test typeof(m_ring) == Healpix.Map{Float32,Healpix.RingOrder,Array{Float32,1}}
@test typeof(m_nest) == Healpix.Map{Float32,Healpix.NestedOrder,Array{Float32,1}}

m_converted = Healpix.nest2ring(m_nest)
@test all(m_converted.pixels .≈ m_ring.pixels)
m_reconverted = Healpix.ring2nest(m_converted)
@test all(m_reconverted.pixels .≈ m_nest.pixels)

m_converted = Healpix.ring2nest(m_ring)
@test all(m_converted.pixels .≈ m_nest.pixels)
