pixels_nside1 = collect(1:12)
pixels_nside2 = collect(1:(12*2^2))
pixels_nside4 = collect(1:(12*4^2))

polmap =
    Healpix.PolarizedMap{Int,Healpix.RingOrder}(pixels_nside1, pixels_nside1, pixels_nside1)
@test polmap.i == pixels_nside1
@test polmap.q == pixels_nside1
@test polmap.u == pixels_nside1


# test the untyped constructor
polmap_new = Healpix.PolarizedMap(polmap.i, polmap.q, polmap.u)
@test polmap.i == polmap_new.i
@test polmap.q == polmap_new.q
@test polmap.u == polmap_new.u


polmap = Healpix.PolarizedMap{Int8,Healpix.RingOrder}(128)

@test length(polmap.i) == length(polmap.q) == length(polmap.u)
@test length(polmap.i) == Healpix.nside2npix(128)

@test_throws ArgumentError("The three I/Q/U vectors must have the same resolution") Healpix.PolarizedMap{
    Int,
    Healpix.RingOrder,
}(
    pixels_nside1,
    pixels_nside2,
    pixels_nside4,
)
