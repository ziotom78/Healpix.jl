@test Healpix.numOfRings(Healpix.Resolution(2)) == 7
@test Healpix.numOfRings(Healpix.Resolution(4)) == 15
@test Healpix.numOfRings(2) == 7
@test Healpix.numOfRings(4) == 15
