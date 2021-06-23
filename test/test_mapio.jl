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

const mapFileName = tempname()
print("Saving $mapFileName\n")
Healpix.saveToFITS(m, "!$mapFileName", typechar = "I")
m2 = Healpix.readMapFromFITS(mapFileName, 1, Int8)
@test m.pixels == m2.pixels


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

# Nothingness/UNSEEN on a *temperature* map

# Create a map with no data associated to pixels
m_nothing = Healpix.HealpixMap{Union{Float32, Nothing}, Healpix.RingOrder}(1)
for i in length(m_nothing)
    @test isnothing(m_nothing[i])
end

nothing_name = tempname(; cleanup=false)
Healpix.saveToFITS(m_nothing, nothing_name)

m_read = Healpix.readMapFromFITS(nothing_name, 1, Union{Float32, Nothing})
@test length(m_read) == length(m_nothing)
for i in length(m_read)
    @test isnothing(m_read[i])
end

# Not using Union{Nothing, T} means that we get UNSEEN pixels
m_unseen_read = Healpix.readMapFromFITS(nothing_name, 1, Float32)
for i in length(m_unseen_read)
    @test m_unseen_read[i] ≈ Healpix.UNSEEN
end

# Nothingness/UNSEEN on a *polarization* map

# Create a map with no data associated to pixels
m_nothing = Healpix.PolarizedHealpixMap{Union{Float32, Nothing}, Healpix.RingOrder}(1)
for i in length(m_nothing.i)
    @test isnothing(m_nothing.i[i])
    @test isnothing(m_nothing.q[i])
    @test isnothing(m_nothing.u[i])
end

nothing_name = tempname(; cleanup=false)
Healpix.saveToFITS(m_nothing, nothing_name)

m_read = Healpix.readPolarizedMapFromFITS(nothing_name, 1, Union{Float32, Nothing})
@test length(m_read.i) == length(m_nothing.i)
for i in length(m_read.i)
    @test isnothing(m_read.i[i])
    @test isnothing(m_read.q[i])
    @test isnothing(m_read.u[i])
end

# Not using Union{Nothing, T} means that we get UNSEEN pixels
m_unseen_read = Healpix.readPolarizedMapFromFITS(nothing_name, 1, Float32)
for i in length(m_unseen_read.i)
    @test m_unseen_read.i[i] ≈ Healpix.UNSEEN
    @test m_unseen_read.q[i] ≈ Healpix.UNSEEN
    @test m_unseen_read.u[i] ≈ Healpix.UNSEEN
end
