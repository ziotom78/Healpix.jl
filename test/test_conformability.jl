@test Healpix.conformables(
    Healpix.Map{Int16,Healpix.RingOrder}(4),
    Healpix.Map{Float32,Healpix.RingOrder}(4),
)
@test Healpix.conformables(
    Healpix.PolarizedMap{Int16,Healpix.RingOrder}(4),
    Healpix.PolarizedMap{Float32,Healpix.RingOrder}(4),
)

# nside mismatch
@test !Healpix.conformables(
    Healpix.Map{Int16,Healpix.RingOrder}(8),
    Healpix.Map{Int16,Healpix.RingOrder}(4),
)
@test !Healpix.conformables(
    Healpix.PolarizedMap{Int16,Healpix.RingOrder}(8),
    Healpix.PolarizedMap{Int16,Healpix.RingOrder}(4),
)
# order mismatch
@test !Healpix.conformables(
    Healpix.Map{Float32,Healpix.RingOrder}(4),
    Healpix.Map{Float32,Healpix.NestedOrder}(4),
)
@test !Healpix.conformables(
    Healpix.PolarizedMap{Float32,Healpix.RingOrder}(4),
    Healpix.PolarizedMap{Float32,Healpix.NestedOrder}(4),
)
