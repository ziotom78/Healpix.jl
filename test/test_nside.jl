@test Healpix.nsideok(16)
@test Healpix.nsideok(1024)
@test !Healpix.nsideok(17)
@test !Healpix.nsideok(-4)

@test Healpix.nside2pixarea(128) ≈ 6.391586616190171e-5
@test Healpix.nside2resol(128) ≈ 0.007994739905831941
@test_throws DomainError(-1, "`NSIDE` is not a positive number") Healpix.nside2pixarea(-1)
@test_throws DomainError(-1, "`NSIDE` is not a positive number") Healpix.nside2resol(-1)
@test_throws DomainError(513, "`NSIDE` is not an integer power of two") Healpix.nside2pixarea(513)
@test_throws DomainError(513, "`NSIDE` is not an integer power of two") Healpix.nside2resol(513)

@test Healpix.nside2npix(4) == 192
@test Healpix.npix2nside(192) == 4
@test_throws DomainError(15, "`NSIDE` is not an integer power of two") Healpix.nside2npix(15)
@test_throws DomainError(7, "Invalid number of pixels") Healpix.npix2nside(7)
@test_throws DomainError(12 * 8 * 9, "Invalid number of pixels") Healpix.npix2nside(12 * 8 * 9)
