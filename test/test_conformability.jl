@test Healpix.conformables(
    Healpix.HealpixMap{Int16,Healpix.RingOrder}(4),
    Healpix.HealpixMap{Float32,Healpix.RingOrder}(4),
)
@test Healpix.conformables(
    Healpix.PolarizedHealpixMap{Int16,Healpix.RingOrder}(4),
    Healpix.PolarizedHealpixMap{Float32,Healpix.RingOrder}(4),
)

# nside mismatch
@test !Healpix.conformables(
    Healpix.HealpixMap{Int16,Healpix.RingOrder}(8),
    Healpix.HealpixMap{Int16,Healpix.RingOrder}(4),
)
@test !Healpix.conformables(
    Healpix.PolarizedHealpixMap{Int16,Healpix.RingOrder}(8),
    Healpix.PolarizedHealpixMap{Int16,Healpix.RingOrder}(4),
)
# order mismatch
@test !Healpix.conformables(
    Healpix.HealpixMap{Float32,Healpix.RingOrder}(4),
    Healpix.HealpixMap{Float32,Healpix.NestedOrder}(4),
)
@test !Healpix.conformables(
    Healpix.PolarizedHealpixMap{Float32,Healpix.RingOrder}(4),
    Healpix.PolarizedHealpixMap{Float32,Healpix.NestedOrder}(4),
)
