pixidx = [1, 5, 3, 2, 2, 3, 3]
tod = [0.0, missing, 2.5, 1.5, 2.5, 3.5, missing]
(binmap, hitmap) = Healpix.tod2map(pixidx, tod, nside = 1, ordering = Healpix.RingOrder)
@test binmap.pixels ≈ [0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test hitmap.pixels == [1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0]

pixidx1 = [1, 3, 2, 2, 3]
tod1 = [0.0, 2.5, 1.5, 2.5, 3.5]
(binmap1, hitmap1) = Healpix.tod2map(pixidx1, tod1, nside = 1, ordering = Healpix.RingOrder)
@test binmap1.pixels ≈ [0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test hitmap1.pixels == [1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0]

binmap2 = Healpix.Map{Float64,Healpix.RingOrder}(1)
binmap2.pixels = [1.0, 0.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
hitmap2 = Healpix.Map{Int,Healpix.RingOrder}(1)
hitmap2.pixels = [1, 6, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0]
Healpix.combinemaps!(binmap1, hitmap1, binmap2, hitmap2)
@test binmap1.pixels ≈ [0.5, 0.5, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test hitmap1.pixels == [2, 8, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0]
