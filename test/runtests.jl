using Healpix
using Base.Test

const eps = 1e-10

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

@test_throws DomainError Healpix.Resolution(0)
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

# pix2angNest

@test_approx_eq_eps Healpix.pix2angNest(resol,      1)[1] 1.5681921571847817 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,      1)[2] 0.7853981633974483 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,      2)[1] 1.5655879699137618 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,      2)[2] 0.7884661249732196 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,      3)[1] 1.5655879699137618 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,      3)[2] 0.7823302018216770 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,      4)[1] 1.5629837473198540 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,      4)[2] 0.7853981633974483 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,  74885)[1] 1.2884125182769750 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,  74885)[2] 2.2396119503130363 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,  75040)[1] 1.2447369771100412 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,  75040)[2] 2.3193789512830896 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,  75195)[1] 1.2198889832038156 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,  75195)[2] 2.2733595276465204 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,  75350)[1] 1.2309594173407747 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,  75350)[2] 2.2549517581918925 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,  75505)[1] 1.2059873663963379 eps
@test_approx_eq_eps Healpix.pix2angNest(resol,  75505)[2] 2.2089323345553233 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 270849)[1] 2.0439875486296080 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 270849)[2] 6.0377483811178836 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 271054)[1] 1.9834500066417189 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 271054)[2] 6.0408163426936552 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 272282)[1] 1.8667651301538279 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 272282)[2] 6.1696707288760484 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 272487)[1] 1.9551931012905357 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 272487)[2] 6.0040008037843995 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 271668)[1] 1.9244782426347979 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 271668)[2] 6.2340979219672459 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 514561)[1] 1.2721038462777619 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 514561)[2] 4.7614763655970300 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 514766)[1] 1.2143369935354960 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 514766)[2] 4.7645443271728016 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 513946)[1] 1.1863995522992576 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 513946)[2] 4.9915734837798764 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 514151)[1] 1.2748275234359654 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 514151)[2] 4.8259035586882275 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 515380)[1] 1.1552980808319346 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 515380)[2] 4.9578259064463923 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 710773)[1] 1.9467798202990929 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 710773)[2] 4.0558452031696355 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 710928)[1] 1.9356052871934553 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 710928)[2] 4.0742529726242633 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 711083)[1] 1.9106332362490186 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 711083)[2] 4.0282335489876937 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 711238)[1] 1.9217036703859778 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 711238)[2] 4.0098257795330658 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 711393)[1] 1.8968556764797520 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 711393)[2] 3.9638063558964967 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 710773)[1] 1.9467798202990929 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 710773)[2] 4.0558452031696355 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 786429)[1] 1.5786089062699391 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 786429)[2] 5.4977871437821380 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 786429)[1] 1.5786089062699391 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 786429)[2] 5.4977871437821380 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 786430)[1] 1.5760046836760313 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 786430)[2] 5.5008551053579096 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 786431)[1] 1.5760046836760313 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 786431)[2] 5.4947191822063663 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 786432)[1] 1.5734004964050115 eps
@test_approx_eq_eps Healpix.pix2angNest(resol, 786432)[2] 5.4977871437821380 eps

# pix2angRing

@test_approx_eq_eps Healpix.pix2angRing(resol,      1)[1] 0.0031894411211228 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,      1)[2] 0.7853981633974483 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,      2)[1] 0.0031894411211113 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,      2)[2] 2.3561944901923448 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,      3)[1] 0.0031894411211113 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,      3)[2] 3.9269908169872414 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,      4)[1] 0.0031894411211113 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,      4)[2] 5.4977871437821380 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,  74885)[1] 0.6290728186524177 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,  74885)[2] 0.0040484441412240 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,  75040)[1] 0.6290728186524177 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,  75040)[2] 1.2590661279206516 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,  75195)[1] 0.6290728186524177 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,  75195)[2] 2.5140838117000794 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,  75350)[1] 0.6290728186524177 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,  75350)[2] 3.7691014954795072 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,  75505)[1] 0.6290728186524177 eps
@test_approx_eq_eps Healpix.pix2angRing(resol,  75505)[2] 5.0241191792589346 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 270849)[1] 1.2557128565725952 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 270849)[2] 0.0000000000000000 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 271054)[1] 1.2557128565725952 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 271054)[2] 1.2578642460662257 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 272282)[1] 1.2584506449956214 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 272282)[2] 2.5126605305566803 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 272487)[1] 1.2584506449956214 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 272487)[2] 3.7705247766229060 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 271668)[1] 1.2557128565725952 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 271668)[2] 5.0253210611133605 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 514561)[1] 1.8858797970171981 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 514561)[2] 0.0000000000000000 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 514766)[1] 1.8858797970171981 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 514766)[2] 1.2578642460662257 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 513946)[1] 1.8831420085941717 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 513946)[2] 2.5126605305566803 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 514151)[1] 1.8831420085941717 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 514151)[2] 3.7705247766229060 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 515380)[1] 1.8858797970171981 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 515380)[2] 5.0253210611133605 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 710773)[1] 2.5125198349373754 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 710773)[2] 0.0040484441412240 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 710928)[1] 2.5125198349373754 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 710928)[2] 1.2590661279206516 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 711083)[1] 2.5125198349373754 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 711083)[2] 2.5140838117000794 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 711238)[1] 2.5125198349373754 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 711238)[2] 3.7691014954795072 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 711393)[1] 2.5125198349373754 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 711393)[2] 5.0241191792589346 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 710773)[1] 2.5125198349373754 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 710773)[2] 0.0040484441412240 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 786429)[1] 3.1384032124686820 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 786429)[2] 0.7853981633974483 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 786429)[1] 3.1384032124686820 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 786429)[2] 0.7853981633974483 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 786430)[1] 3.1384032124686820 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 786430)[2] 2.3561944901923448 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 786431)[1] 3.1384032124686820 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 786431)[2] 3.9269908169872414 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 786432)[1] 3.1384032124686820 eps
@test_approx_eq_eps Healpix.pix2angRing(resol, 786432)[2] 5.4977871437821380 eps

# Conformability

@test Healpix.conformables(Map{Int16, RingOrder}(4),
                           Map{Float32, RingOrder}(4))
@test ! Healpix.conformables(Map{Int16, RingOrder}(8),
                             Map{Int16, RingOrder}(4)) # nside mismatch
@test ! Healpix.conformables(Map{Float32, RingOrder}(4),
                             Map{Float32, NestedOrder}(4)) # order mismatch

# Map loading

m = Healpix.readMapFromFITS("float_map.fits", 1, Float32)
@test typeof(m) == Healpix.Map{Float32, Healpix.RingOrder}
@test m.resolution.nside == 4
@test m.pixels == [Float32(x) for x in 0:(12*4^2 - 1)]

m = Healpix.readMapFromFITS("int_map.fits", 1, Int8)
@test typeof(m) == Healpix.Map{Int8, Healpix.RingOrder}
@test m.resolution.nside == 1
@test m.pixels == [Int8(x) for x in 0:11]

# Map saving

const mapFileName = tempname()
print("Saving $mapFileName\n")
Healpix.saveToFITS(m, "!$mapFileName", "I")
m2 = Healpix.readMapFromFITS(mapFileName, 1, Int8)
@test m.pixels == m2.pixels

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

alm = Healpix.readAlmFromFITS("alm.fits", Complex128)
@test_approx_eq_eps alm[1] (5.443205775735e+03 + 0.000000000000e+00im) eps
@test_approx_eq_eps alm[2] (-3.143659646589e+03 + 0.000000000000e+00im) eps
@test_approx_eq_eps alm[3] (-8.445976910202e-07 + 0.000000000000e+00im) eps
@test_approx_eq_eps alm[4] (3.003475555079e-07 + 0.000000000000e+00im) eps
@test_approx_eq_eps alm[5] (-1.094164444296e-06 + 0.000000000000e+00im) eps
@test_approx_eq_eps alm[6] (3.745732939005e-07 + 0.000000000000e+00im) eps
@test_approx_eq_eps alm[7] (-1.344818023454e-06 + 0.000000000000e+00im) eps
@test_approx_eq_eps alm[8] (6.658742467775e-01 + -3.280201017201e+01im) eps
@test_approx_eq_eps alm[9] (1.200156696497e-15 + -1.483545355539e-15im) eps
@test_approx_eq_eps alm[10] (-1.959362683381e-01 + -3.069662938925e+00im) eps
@test_approx_eq_eps alm[11] (-4.068874968688e-15 + -5.766800364761e-16im) eps
@test_approx_eq_eps alm[12] (-1.134265758838e-01 + 1.102692148708e+00im) eps
@test_approx_eq_eps alm[13] (-2.153997451326e-15 + -5.152134573109e-15im) eps
@test_approx_eq_eps alm[14] (-6.895464505703e-01 + 1.597235776935e+01im) eps
@test_approx_eq_eps alm[15] (4.344648310121e-15 + -1.892253965354e-15im) eps
@test_approx_eq_eps alm[16] (6.488192068017e-02 + 3.967268430232e+00im) eps
@test_approx_eq_eps alm[17] (1.340722726247e-15 + -9.822172672993e-15im) eps
@test_approx_eq_eps alm[18] (1.485579846001e-01 + 7.942515768792e-01im) eps
@test_approx_eq_eps alm[19] (6.938894269949e-01 + -1.028821633749e+01im) eps
@test_approx_eq_eps alm[20] (3.412745093853e-15 + -4.420258002639e-15im) eps
@test_approx_eq_eps alm[21] (3.712279444471e-02 + -3.564529788432e+00im) eps
@test_approx_eq_eps alm[22] (-6.205890749143e-15 + 4.855843625358e-16im) eps
@test_approx_eq_eps alm[23] (-6.893793878861e-01 + 7.453657075068e+00im) eps
@test_approx_eq_eps alm[24] (4.644726834988e-11 + -4.685776289892e-15im) eps
@test_approx_eq_eps alm[25] (-1.162363120645e-01 + 3.084870414802e+00im) eps
@test_approx_eq_eps alm[26] (6.806207052149e-01 + -5.769903744427e+00im) eps
@test_approx_eq_eps alm[27] (-1.307082489340e-15 + 1.143824238899e-14im) eps
@test_approx_eq_eps alm[28] (-6.698490836781e-01 + 4.661675665246e+00im) eps
