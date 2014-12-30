using Healpix
using Base.Test

@test Healpix.ilog2(convert(Uint32, 1)) == 0
@test Healpix.ilog2(convert(Uint32, 6)) == 2
@test Healpix.ilog2(convert(Uint32, 1023)) == 9
@test Healpix.ilog2(convert(Uint32, 1024)) == 10
@test Healpix.ilog2(convert(Uint32, 8194)) == 13
@test Healpix.ilog2(convert(Uint32, 131124)) == 17

# nside2npix and npix2nside

@test Healpix.nside2npix(4) == 192
@test Healpix.npix2nside(192) == 4
@test_throws DomainError Healpix.nside2npix(15)
@test_throws DomainError Healpix.npix2nside(7)
@test_throws DomainError Healpix.npix2nside(12 * 8 * 9)

# ang2pixNest

@test_throws DomainError Healpix.Resolution(-5)
@test_throws DomainError Healpix.Resolution(100000000)
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

# ang2pixRing

@test Healpix.ang2pixRing(resol, 0.0000000000000000, 0.0000000000000000) ==      1
@test Healpix.ang2pixRing(resol, 0.0000000000000000, 1.2566370614359172) ==      1
@test Healpix.ang2pixRing(resol, 0.0000000000000000, 2.5132741228718345) ==      2
@test Healpix.ang2pixRing(resol, 0.0000000000000000, 3.7699111843077517) ==      3
@test Healpix.ang2pixRing(resol, 0.0000000000000000, 5.0265482457436690) ==      4
@test Healpix.ang2pixRing(resol, 0.0000000000000000, 6.2831853071795862) ==      1
@test Healpix.ang2pixRing(resol, 0.6283185307179586, 0.0000000000000000) ==  74885
@test Healpix.ang2pixRing(resol, 0.6283185307179586, 1.2566370614359172) ==  75040
@test Healpix.ang2pixRing(resol, 0.6283185307179586, 2.5132741228718345) ==  75195
@test Healpix.ang2pixRing(resol, 0.6283185307179586, 3.7699111843077517) ==  75350
@test Healpix.ang2pixRing(resol, 0.6283185307179586, 5.0265482457436690) ==  75505
@test Healpix.ang2pixRing(resol, 0.6283185307179586, 6.2831853071795862) ==  74885
@test Healpix.ang2pixRing(resol, 1.2566370614359172, 0.0000000000000000) == 270849
@test Healpix.ang2pixRing(resol, 1.2566370614359172, 1.2566370614359172) == 271054
@test Healpix.ang2pixRing(resol, 1.2566370614359172, 2.5132741228718345) == 272282
@test Healpix.ang2pixRing(resol, 1.2566370614359172, 3.7699111843077517) == 272487
@test Healpix.ang2pixRing(resol, 1.2566370614359172, 5.0265482457436690) == 271668
@test Healpix.ang2pixRing(resol, 1.2566370614359172, 6.2831853071795862) == 270849
@test Healpix.ang2pixRing(resol, 1.8849555921538759, 0.0000000000000000) == 514561
@test Healpix.ang2pixRing(resol, 1.8849555921538759, 1.2566370614359172) == 514766
@test Healpix.ang2pixRing(resol, 1.8849555921538759, 2.5132741228718345) == 513946
@test Healpix.ang2pixRing(resol, 1.8849555921538759, 3.7699111843077517) == 514151
@test Healpix.ang2pixRing(resol, 1.8849555921538759, 5.0265482457436690) == 515380
@test Healpix.ang2pixRing(resol, 1.8849555921538759, 6.2831853071795862) == 514561
@test Healpix.ang2pixRing(resol, 2.5132741228718345, 0.0000000000000000) == 710773
@test Healpix.ang2pixRing(resol, 2.5132741228718345, 1.2566370614359172) == 710928
@test Healpix.ang2pixRing(resol, 2.5132741228718345, 2.5132741228718345) == 711083
@test Healpix.ang2pixRing(resol, 2.5132741228718345, 3.7699111843077517) == 711238
@test Healpix.ang2pixRing(resol, 2.5132741228718345, 5.0265482457436690) == 711393
@test Healpix.ang2pixRing(resol, 2.5132741228718345, 6.2831853071795862) == 710773
@test Healpix.ang2pixRing(resol, 3.1415926535897931, 0.0000000000000000) == 786429
@test Healpix.ang2pixRing(resol, 3.1415926535897931, 1.2566370614359172) == 786429
@test Healpix.ang2pixRing(resol, 3.1415926535897931, 2.5132741228718345) == 786430
@test Healpix.ang2pixRing(resol, 3.1415926535897931, 3.7699111843077517) == 786431
@test Healpix.ang2pixRing(resol, 3.1415926535897931, 5.0265482457436690) == 786432
@test Healpix.ang2pixRing(resol, 3.1415926535897931, 6.2831853071795862) == 786429

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

alm = Healpix.Alm{Complex64}(10, 8)
@test almIndex(alm, 4, 2) == 24
@test almIndex(alm, 5, 2) == 25
@test almIndex(alm, 5, 3) == 33
@test almIndex(alm, [4, 6, 5], [3, 4, 5]) == [32, 41, 46]

alm = Healpix.Alm("alm.fits", Complex128)
@test_approx_eq alm[1] (5.443205775735e+03 + 0.000000000000e+00im)
@test_approx_eq alm[2] (-3.143659646589e+03 + 0.000000000000e+00im)
@test_approx_eq alm[3] (-8.445976910202e-07 + 0.000000000000e+00im)
@test_approx_eq alm[4] (3.003475555079e-07 + 0.000000000000e+00im)
@test_approx_eq alm[5] (-1.094164444296e-06 + 0.000000000000e+00im)
@test_approx_eq alm[6] (3.745732939005e-07 + 0.000000000000e+00im)
@test_approx_eq alm[7] (-1.344818023454e-06 + 0.000000000000e+00im)
@test_approx_eq alm[8] (6.658742467775e-01 + -3.280201017201e+01im)
@test_approx_eq alm[9] (1.200156696497e-15 + -1.483545355539e-15im)
@test_approx_eq alm[10] (-1.959362683381e-01 + -3.069662938925e+00im)
@test_approx_eq alm[11] (-4.068874968688e-15 + -5.766800364761e-16im)
@test_approx_eq alm[12] (-1.134265758838e-01 + 1.102692148708e+00im)
@test_approx_eq alm[13] (-2.153997451326e-15 + -5.152134573109e-15im)
@test_approx_eq alm[14] (-6.895464505703e-01 + 1.597235776935e+01im)
@test_approx_eq alm[15] (4.344648310121e-15 + -1.892253965354e-15im)
@test_approx_eq alm[16] (6.488192068017e-02 + 3.967268430232e+00im)
@test_approx_eq alm[17] (1.340722726247e-15 + -9.822172672993e-15im)
@test_approx_eq alm[18] (1.485579846001e-01 + 7.942515768792e-01im)
@test_approx_eq alm[19] (6.938894269949e-01 + -1.028821633749e+01im)
@test_approx_eq alm[20] (3.412745093853e-15 + -4.420258002639e-15im)
@test_approx_eq alm[21] (3.712279444471e-02 + -3.564529788432e+00im)
@test_approx_eq alm[22] (-6.205890749143e-15 + 4.855843625358e-16im)
@test_approx_eq alm[23] (-6.893793878861e-01 + 7.453657075068e+00im)
@test_approx_eq alm[24] (4.644726834988e-11 + -4.685776289892e-15im)
@test_approx_eq alm[25] (-1.162363120645e-01 + 3.084870414802e+00im)
@test_approx_eq alm[26] (6.806207052149e-01 + -5.769903744427e+00im)
@test_approx_eq alm[27] (-1.307082489340e-15 + 1.143824238899e-14im)
@test_approx_eq alm[28] (-6.698490836781e-01 + 4.661675665246e+00im)
