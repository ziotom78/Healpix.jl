@test Healpix.ilog2(1) == 0
@test Healpix.ilog2(6) == 2
@test Healpix.ilog2(1023) == 9
@test Healpix.ilog2(1024) == 10
@test Healpix.ilog2(8194) == 13
@test Healpix.ilog2(131124) == 17
