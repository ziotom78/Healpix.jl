pixels_nside1 = collect(1:12)
pixels_nside2 = collect(1:(12 * 2^2))
pixels_nside4 = collect(1:(12 * 4^2))

polmap = Healpix.PolarizedMap{Int, Healpix.RingOrder}(
    pixels_nside1,
    pixels_nside1,
    pixels_nside1,
)
@test polmap.i == pixels_nside1
@test polmap.q == pixels_nside1
@test polmap.u == pixels_nside1
                        
polmap = Healpix.PolarizedMap{Int8, Healpix.RingOrder}(128)

@test length(polmap.i) == length(polmap.q) == length(polmap.u)
@test length(polmap.i) == Healpix.nside2npix(128)

@test_throws ArgumentError("The three I/Q/U maps must have the same resolution") Healpix.PolarizedMap{Int, Healpix.RingOrder}(
    pixels_nside1,
    pixels_nside2,
    pixels_nside4,
)
