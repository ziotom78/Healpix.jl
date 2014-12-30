using Healpix
using Base.Test

@test Healpix.ilog2(convert(Uint32, 1024)) == 10

# nside2npix and npix2nside

@test Healpix.nside2npix(4) == 192
@test Healpix.npix2nside(192) == 4
@test_throws DomainError Healpix.nside2npix(15)
@test_throws DomainError Healpix.npix2nside(7)

# ang2pixNest

resol = Healpix.Resolution(256)

@test Healpix.ang2pixNest(resol, 0.0000000000000000, 0.0000000000000000) ==  65536
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 1.2566370614359172) ==  65536
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 2.5132741228718345) == 131072
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 3.7699111843077517) == 196608
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 5.0265482457436690) == 262144
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 6.2831853071795862) ==  65536

@test Healpix.ang2pixNest(resol, 0.6283185307179586, 0.0000000000000000) ==  45055
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 1.2566370614359172) ==  31074
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 2.5132741228718345) == 116111
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 3.7699111843077517) == 182862
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 5.0265482457436690) == 243347
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 6.2831853071795862) ==  45055

@test Healpix.ang2pixNest(resol, 1.2566370614359172, 0.0000000000000000) == 315344
@test Healpix.ang2pixNest(resol, 1.2566370614359172, 1.2566370614359172) == 387305
@test Healpix.ang2pixNest(resol, 1.2566370614359172, 2.5132741228718345) ==  71955
@test Healpix.ang2pixNest(resol, 1.2566370614359172, 3.7699111843077517) == 140834
@test Healpix.ang2pixNest(resol, 1.2566370614359172, 5.0265482457436690) == 513237
@test Healpix.ang2pixNest(resol, 1.2566370614359172, 6.2831853071795862) == 315344

@test Healpix.ang2pixNest(resol, 1.8849555921538759, 0.0000000000000000) == 274481
@test Healpix.ang2pixNest(resol, 1.8849555921538759, 1.2566370614359172) == 338732
@test Healpix.ang2pixNest(resol, 1.8849555921538759, 2.5132741228718345) == 645599
@test Healpix.ang2pixNest(resol, 1.8849555921538759, 3.7699111843077517) == 714478
@test Healpix.ang2pixNest(resol, 1.8849555921538759, 5.0265482457436690) == 464664
@test Healpix.ang2pixNest(resol, 1.8849555921538759, 6.2831853071795862) == 274481

@test Healpix.ang2pixNest(resol, 2.5132741228718345, 0.0000000000000000) == 565251
@test Healpix.ang2pixNest(resol, 2.5132741228718345, 1.2566370614359172) == 543086
@test Healpix.ang2pixNest(resol, 2.5132741228718345, 2.5132741228718345) == 603571
@test Healpix.ang2pixNest(resol, 2.5132741228718345, 3.7699111843077517) == 670322
@test Healpix.ang2pixNest(resol, 2.5132741228718345, 5.0265482457436690) == 755359
@test Healpix.ang2pixNest(resol, 2.5132741228718345, 6.2831853071795862) == 565251

@test Healpix.ang2pixNest(resol, 3.1415926535897931, 0.0000000000000000) == 524289
@test Healpix.ang2pixNest(resol, 3.1415926535897931, 1.2566370614359172) == 524289
@test Healpix.ang2pixNest(resol, 3.1415926535897931, 2.5132741228718345) == 589825
@test Healpix.ang2pixNest(resol, 3.1415926535897931, 3.7699111843077517) == 655361
@test Healpix.ang2pixNest(resol, 3.1415926535897931, 5.0265482457436690) == 720897
@test Healpix.ang2pixNest(resol, 3.1415926535897931, 6.2831853071795862) == 524289

# Map loading

m = Healpix.Map("int_map.fits", 1, Int8)
@test m.resolution.nside == 1
@test m.ordering == Healpix.Ring
@test m.pixels == [int8(x) for x in 0:11]

m = Healpix.Map("float_map.fits", 1, Float32)
@test m.resolution.nside == 4
@test m.ordering == Healpix.Ring
@test m.pixels == [float32(x) for x in 0:(12*4^2 - 1)]

# Alm creation

@test Healpix.numberOfAlms(10, 5) == 51
@test Healpix.numberOfAlms(10, 7) == 60
@test Healpix.numberOfAlms(12, 7) == 76
@test Healpix.numberOfAlms(12, 12) == 91
@test_throws DomainError Healpix.numberOfAlms(-1, 1)
@test_throws DomainError Healpix.numberOfAlms(4, -1)
@test_throws DomainError Healpix.numberOfAlms(5, 7)

alm = Healpix.Alm{Float64}(10, 10)
