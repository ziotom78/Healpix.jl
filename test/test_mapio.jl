# Map loading

m = Healpix.readMapFromFITS("float_map.fits", 1, Float32)
@test typeof(m) == Healpix.HealpixMap{Float32,Healpix.RingOrder,Array{Float32,1}}
@test m.resolution.nside == 4
@test m.pixels == [Float32(x) for x = 0:(12*4^2-1)]

m = Healpix.readMapFromFITS("int_map.fits", 1, Int8)
@test typeof(m) == Healpix.HealpixMap{Int8,Healpix.RingOrder,Array{Int8,1}}
@test m.resolution.nside == 1
@test m.pixels == [Int8(x) for x = 0:11]

# Map saving

mapFileName = tempname()
Healpix.saveToFITS(m, "!$mapFileName", typechar = "I")
m2 = Healpix.readMapFromFITS(mapFileName, 1, Int8)
@test m.pixels == m2.pixels

# Polarized map loading

m = Healpix.readPolarizedMapFromFITS("pol_float_map.fits", 1, Float32)
@test typeof(m) == Healpix.PolarizedHealpixMap{Float32,Healpix.RingOrder,Array{Float32,1}}
@test m.i.resolution.nside == 1
@test m.q.resolution.nside == 1
@test m.u.resolution.nside == 1
@test m.i == [Float32(x) for x = 0:11]
@test m.q == [Float32(x) for x = 12:23]
@test m.u == [Float32(x) for x = 24:35]

m = Healpix.readPolarizedMapFromFITS("pol_int_map.fits", 1, Int8)
@test typeof(m) == Healpix.PolarizedHealpixMap{Int8,Healpix.RingOrder,Array{Int8,1}}
@test m.i.resolution.nside == 1
@test m.q.resolution.nside == 1
@test m.u.resolution.nside == 1
@test m.i == [Int8(x) for x = 0:11]
@test m.q == [Int8(x) for x = 12:23]
@test m.u == [Int8(x) for x = 24:35]

# Polarized map saving

mapFileName = tempname()
Healpix.saveToFITS(m, "!$mapFileName", typechar = "I")
m2 = Healpix.readPolarizedMapFromFITS(mapFileName, 1, Int8)
@test m.i == m2.i
@test m.q == m2.q
@test m.u == m2.u

# Map reordering
m_ring = Healpix.readMapFromFITS("float_map.fits", 1, Float32)
m_nest = Healpix.readMapFromFITS("float_map_nest.fits", 1, Float32)
@test typeof(m_ring) == Healpix.HealpixMap{Float32,Healpix.RingOrder,Array{Float32,1}}
@test typeof(m_nest) == Healpix.HealpixMap{Float32,Healpix.NestedOrder,Array{Float32,1}}

m_converted = Healpix.nest2ring(m_nest)
@test all(m_converted.pixels .≈ m_ring.pixels)
m_reconverted = Healpix.ring2nest(m_converted)
@test all(m_reconverted.pixels .≈ m_nest.pixels)

m_converted = Healpix.ring2nest(m_ring)
@test all(m_converted.pixels .≈ m_nest.pixels)

# Columns (see issue #61, https://github.com/ziotom78/Healpix.jl/issues/61)

map = Healpix.readPolarizedMapFromFITS("many_columns.fits", 2, Float32)
@test m.i == [Float32(x) for x = 0:11]
@test m.q == [Float32(x) for x = 12:23]
@test m.u == [Float32(x) for x = 24:35]
