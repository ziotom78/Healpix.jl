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
