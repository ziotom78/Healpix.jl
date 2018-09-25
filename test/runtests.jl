# -*- encoding: utf-8 -*-

import Healpix
import Cairo
using Test

const eps = 1e-10

@test Healpix.ilog2(1) == 0
@test Healpix.ilog2(6) == 2
@test Healpix.ilog2(1023) == 9
@test Healpix.ilog2(1024) == 10
@test Healpix.ilog2(8194) == 13
@test Healpix.ilog2(131124) == 17

# nsideok, nside2pixarea, nside2resol, nside2npix, npix2nside

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

# ang2vec

x, y, z = Healpix.ang2vec(2.556481306546430, 1.573953117024192)
@test x ≈ -0.001743467757429 atol = eps
@test y ≈ 0.552289458399410 atol = eps
@test z ≈ -0.833650594950345 atol = eps

x, y, z = Healpix.ang2vec(0.140282766372135, 0.233893209450599)
@test x ≈ 0.136015928900898 atol = eps
@test y ≈ 0.032406308842171 atol = eps
@test z ≈ 0.990176498525617 atol = eps

x, y, z = Healpix.ang2vec(0.441131335644146, 3.464317397342697)
@test x ≈ -0.404920724073465 atol = eps
@test y ≈ -0.135412016604473 atol = eps
@test z ≈ 0.904269203818714 atol = eps

x, y, z = Healpix.ang2vec(2.339362129743034, 2.415312683457592)
@test x ≈ -0.537491821098482 atol = eps
@test y ≈ 0.477421955393414 atol = eps
@test z ≈ -0.695104897666939 atol = eps

x, y, z = Healpix.ang2vec(2.012884650691451, 0.429579124174392)
@test x ≈ 0.821736466784090 atol = eps
@test y ≈ 0.376447107630774 atol = eps
@test z ≈ -0.427827949430170 atol = eps

x, y, z = Healpix.ang2vec(0.698567538240779, 2.951569899700352)
@test x ≈ -0.631545193788522 atol = eps
@test y ≈ 0.121473570759142 atol = eps
@test z ≈ 0.765764219462911 atol = eps

x, y, z = Healpix.ang2vec(0.602746245881256, 2.363106505915441)
@test x ≈ -0.403624903017147 atol = eps
@test y ≈ 0.398083394409416 atol = eps
@test z ≈ 0.823781857508331 atol = eps

x, y, z = Healpix.ang2vec(3.106081268871312, 4.592838732778776)
@test x ≈ -0.004234399242699 atol = eps
@test y ≈ -0.035250507901764 atol = eps
@test z ≈ -0.999369537036086 atol = eps

x, y, z = Healpix.ang2vec(1.906789620733219, 4.322731382956660)
@test x ≈ -0.358630521788042 atol = eps
@test y ≈ -0.873314025552965 atol = eps
@test z ≈ -0.329707084568277 atol = eps

x, y, z = Healpix.ang2vec(0.637907993514304, 4.877925523463614)
@test x ≈ 0.098130081527976 atol = eps
@test y ≈ -0.587375491097860 atol = eps
@test z ≈ 0.803343338527719 atol = eps

# vec2ang

theta, phi = Healpix.vec2ang(2.479973695958578, 2.540405094768749, 1.360043653263107)
@test theta ≈ 1.204952691410965 atol = eps
@test phi ≈ 0.797434801040588 atol = eps

theta, phi = Healpix.vec2ang(1.095479074939707, 0.898675220740957, 0.147081790022168)
@test theta ≈ 1.467363643308421 atol = eps
@test phi ≈ 0.687026327338059 atol = eps

theta, phi = Healpix.vec2ang(0.576316747680224, 0.592397935683216, 0.829922061194829)
@test theta ≈ 0.783322736466122 atol = eps
@test phi ≈ 0.799157016147229 atol = eps

theta, phi = Healpix.vec2ang(2.502045018710797, 1.228308248258689, 0.292075568023260)
@test theta ≈ 1.466388897914957 atol = eps
@test phi ≈ 0.456358646456269 atol = eps

theta, phi = Healpix.vec2ang(0.784144561032690, 2.203400212521606, 2.645182960828146)
@test theta ≈ 0.723995999159878 atol = eps
@test phi ≈ 1.228893440982002 atol = eps

theta, phi = Healpix.vec2ang(2.313385295284291, 2.310829312494675, 0.909839585819006)
@test theta ≈ 1.299407206361556 atol = eps
@test phi ≈ 0.784845424663175 atol = eps

theta, phi = Healpix.vec2ang(2.028549527136607, 2.806016077249088, 0.732470589106458)
@test theta ≈ 1.362324539502480 atol = eps
@test phi ≈ 0.944847150958497 atol = eps

theta, phi = Healpix.vec2ang(2.806692561465244, 1.750787884540839, 1.617079425838591)
@test theta ≈ 1.116115824440890 atol = eps
@test phi ≈ 0.557729020009978 atol = eps

theta, phi = Healpix.vec2ang(1.973555660995911, 1.086043406239103, 2.959437359715485)
@test theta ≈ 0.650613947401187 atol = eps
@test phi ≈ 0.503071845787347 atol = eps

theta, phi = Healpix.vec2ang(2.922856075911509, 2.304551510739625, 0.874687489108078)
@test theta ≈ 1.339986032596264 atol = eps
@test phi ≈ 0.667663882320130 atol = eps

# ang2pixNest

@test_throws DomainError(0) Healpix.Resolution(0)
@test_throws DomainError(100000000) Healpix.Resolution(100000000)
lowresol = Healpix.Resolution(2)
resol = Healpix.Resolution(256)

@test Healpix.ang2pixNest(resol, 0.0000000000000000, 0.0000000000000000) ==  65536
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 1.2566370614359172) ==  65536
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 2.5132741228718345) == 131072
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 3.7699111843077517) == 196608
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 5.0265482457436690) == 262144
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 6.2831853071795862) == 262144

@test Healpix.ang2pixNest(resol, 0.6283185307179586, 0.0000000000000000) ==  45055
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 1.2566370614359172) ==  31074
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 2.5132741228718345) == 116111
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 3.7699111843077517) == 182862
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 5.0265482457436690) == 243347
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 6.2831853071795862) == 221182

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
@test Healpix.ang2pixNest(resol, 2.5132741228718345, 6.2831853071795862) == 741378

@test Healpix.ang2pixNest(resol, 3.1415926535897931, 0.0000000000000000) == 524289
@test Healpix.ang2pixNest(resol, 3.1415926535897931, 1.2566370614359172) == 524289
@test Healpix.ang2pixNest(resol, 3.1415926535897931, 2.5132741228718345) == 589825
@test Healpix.ang2pixNest(resol, 3.1415926535897931, 3.7699111843077517) == 655361
@test Healpix.ang2pixNest(resol, 3.1415926535897931, 5.0265482457436690) == 720897
@test Healpix.ang2pixNest(resol, 3.1415926535897931, 6.2831853071795862) == 720897

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

@test Healpix.pix2angNest(lowresol,  1)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angNest(lowresol,  1)[2] ≈ 0.785398163397448 atol = eps
@test Healpix.pix2angNest(lowresol,  2)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angNest(lowresol,  2)[2] ≈ 1.178097245096172 atol = eps
@test Healpix.pix2angNest(lowresol,  3)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angNest(lowresol,  3)[2] ≈ 0.392699081698724 atol = eps
@test Healpix.pix2angNest(lowresol,  4)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angNest(lowresol,  4)[2] ≈ 0.785398163397448 atol = eps
@test Healpix.pix2angNest(lowresol,  5)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angNest(lowresol,  5)[2] ≈ 2.356194490192345 atol = eps
@test Healpix.pix2angNest(lowresol,  6)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angNest(lowresol,  6)[2] ≈ 2.748893571891069 atol = eps
@test Healpix.pix2angNest(lowresol,  7)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angNest(lowresol,  7)[2] ≈ 1.963495408493621 atol = eps
@test Healpix.pix2angNest(lowresol,  8)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angNest(lowresol,  8)[2] ≈ 2.356194490192345 atol = eps
@test Healpix.pix2angNest(lowresol,  9)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angNest(lowresol,  9)[2] ≈ 3.926990816987241 atol = eps
@test Healpix.pix2angNest(lowresol, 10)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angNest(lowresol, 10)[2] ≈ 4.319689898685965 atol = eps
@test Healpix.pix2angNest(lowresol, 11)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angNest(lowresol, 11)[2] ≈ 3.534291735288517 atol = eps
@test Healpix.pix2angNest(lowresol, 12)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angNest(lowresol, 12)[2] ≈ 3.926990816987241 atol = eps
@test Healpix.pix2angNest(lowresol, 13)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angNest(lowresol, 13)[2] ≈ 5.497787143782138 atol = eps
@test Healpix.pix2angNest(lowresol, 14)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angNest(lowresol, 14)[2] ≈ 5.890486225480862 atol = eps
@test Healpix.pix2angNest(lowresol, 15)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angNest(lowresol, 15)[2] ≈ 5.105088062083414 atol = eps
@test Healpix.pix2angNest(lowresol, 16)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angNest(lowresol, 16)[2] ≈ 5.497787143782138 atol = eps
@test Healpix.pix2angNest(lowresol, 17)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angNest(lowresol, 17)[2] ≈ 0.000000000000000 atol = eps
@test Healpix.pix2angNest(lowresol, 18)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angNest(lowresol, 18)[2] ≈ 0.392699081698724 atol = eps
@test Healpix.pix2angNest(lowresol, 19)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angNest(lowresol, 19)[2] ≈ 5.890486225480862 atol = eps
@test Healpix.pix2angNest(lowresol, 20)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angNest(lowresol, 20)[2] ≈ 0.000000000000000 atol = eps
@test Healpix.pix2angNest(lowresol, 21)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angNest(lowresol, 21)[2] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angNest(lowresol, 22)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angNest(lowresol, 22)[2] ≈ 1.963495408493621 atol = eps
@test Healpix.pix2angNest(lowresol, 23)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angNest(lowresol, 23)[2] ≈ 1.178097245096172 atol = eps
@test Healpix.pix2angNest(lowresol, 24)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angNest(lowresol, 24)[2] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angNest(lowresol, 25)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angNest(lowresol, 25)[2] ≈ 3.141592653589793 atol = eps
@test Healpix.pix2angNest(lowresol, 26)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angNest(lowresol, 26)[2] ≈ 3.534291735288517 atol = eps
@test Healpix.pix2angNest(lowresol, 27)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angNest(lowresol, 27)[2] ≈ 2.748893571891069 atol = eps
@test Healpix.pix2angNest(lowresol, 28)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angNest(lowresol, 28)[2] ≈ 3.141592653589793 atol = eps
@test Healpix.pix2angNest(lowresol, 29)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angNest(lowresol, 29)[2] ≈ 4.712388980384690 atol = eps
@test Healpix.pix2angNest(lowresol, 30)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angNest(lowresol, 30)[2] ≈ 5.105088062083414 atol = eps
@test Healpix.pix2angNest(lowresol, 31)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angNest(lowresol, 31)[2] ≈ 4.319689898685965 atol = eps
@test Healpix.pix2angNest(lowresol, 32)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angNest(lowresol, 32)[2] ≈ 4.712388980384690 atol = eps
@test Healpix.pix2angNest(lowresol, 33)[1] ≈ 2.730454791267445 atol = eps
@test Healpix.pix2angNest(lowresol, 33)[2] ≈ 0.785398163397448 atol = eps
@test Healpix.pix2angNest(lowresol, 34)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angNest(lowresol, 34)[2] ≈ 1.178097245096172 atol = eps
@test Healpix.pix2angNest(lowresol, 35)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angNest(lowresol, 35)[2] ≈ 0.392699081698724 atol = eps
@test Healpix.pix2angNest(lowresol, 36)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angNest(lowresol, 36)[2] ≈ 0.785398163397448 atol = eps
@test Healpix.pix2angNest(lowresol, 37)[1] ≈ 2.730454791267445 atol = eps
@test Healpix.pix2angNest(lowresol, 37)[2] ≈ 2.356194490192345 atol = eps
@test Healpix.pix2angNest(lowresol, 38)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angNest(lowresol, 38)[2] ≈ 2.748893571891069 atol = eps
@test Healpix.pix2angNest(lowresol, 39)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angNest(lowresol, 39)[2] ≈ 1.963495408493621 atol = eps
@test Healpix.pix2angNest(lowresol, 40)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angNest(lowresol, 40)[2] ≈ 2.356194490192345 atol = eps
@test Healpix.pix2angNest(lowresol, 41)[1] ≈ 2.730454791267445 atol = eps
@test Healpix.pix2angNest(lowresol, 41)[2] ≈ 3.926990816987241 atol = eps
@test Healpix.pix2angNest(lowresol, 42)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angNest(lowresol, 42)[2] ≈ 4.319689898685965 atol = eps
@test Healpix.pix2angNest(lowresol, 43)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angNest(lowresol, 43)[2] ≈ 3.534291735288517 atol = eps
@test Healpix.pix2angNest(lowresol, 44)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angNest(lowresol, 44)[2] ≈ 3.926990816987241 atol = eps
@test Healpix.pix2angNest(lowresol, 45)[1] ≈ 2.730454791267445 atol = eps
@test Healpix.pix2angNest(lowresol, 45)[2] ≈ 5.497787143782138 atol = eps
@test Healpix.pix2angNest(lowresol, 46)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angNest(lowresol, 46)[2] ≈ 5.890486225480862 atol = eps
@test Healpix.pix2angNest(lowresol, 47)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angNest(lowresol, 47)[2] ≈ 5.105088062083414 atol = eps

@test Healpix.pix2angNest(resol,      1)[1] ≈ 1.5681921571847817 atol = eps
@test Healpix.pix2angNest(resol,      1)[2] ≈ 0.7853981633974483 atol = eps
@test Healpix.pix2angNest(resol,      2)[1] ≈ 1.5655879699137618 atol = eps
@test Healpix.pix2angNest(resol,      2)[2] ≈ 0.7884661249732196 atol = eps
@test Healpix.pix2angNest(resol,      3)[1] ≈ 1.5655879699137618 atol = eps
@test Healpix.pix2angNest(resol,      3)[2] ≈ 0.7823302018216770 atol = eps
@test Healpix.pix2angNest(resol,      4)[1] ≈ 1.5629837473198540 atol = eps
@test Healpix.pix2angNest(resol,      4)[2] ≈ 0.7853981633974483 atol = eps
@test Healpix.pix2angNest(resol,  74885)[1] ≈ 1.2884125182769750 atol = eps
@test Healpix.pix2angNest(resol,  74885)[2] ≈ 2.2396119503130363 atol = eps
@test Healpix.pix2angNest(resol,  75040)[1] ≈ 1.2447369771100412 atol = eps
@test Healpix.pix2angNest(resol,  75040)[2] ≈ 2.3193789512830896 atol = eps
@test Healpix.pix2angNest(resol,  75195)[1] ≈ 1.2198889832038156 atol = eps
@test Healpix.pix2angNest(resol,  75195)[2] ≈ 2.2733595276465204 atol = eps
@test Healpix.pix2angNest(resol,  75350)[1] ≈ 1.2309594173407747 atol = eps
@test Healpix.pix2angNest(resol,  75350)[2] ≈ 2.2549517581918925 atol = eps
@test Healpix.pix2angNest(resol,  75505)[1] ≈ 1.2059873663963379 atol = eps
@test Healpix.pix2angNest(resol,  75505)[2] ≈ 2.2089323345553233 atol = eps
@test Healpix.pix2angNest(resol, 270849)[1] ≈ 2.0439875486296080 atol = eps
@test Healpix.pix2angNest(resol, 270849)[2] ≈ 6.0377483811178836 atol = eps
@test Healpix.pix2angNest(resol, 271054)[1] ≈ 1.9834500066417189 atol = eps
@test Healpix.pix2angNest(resol, 271054)[2] ≈ 6.0408163426936552 atol = eps
@test Healpix.pix2angNest(resol, 272282)[1] ≈ 1.8667651301538279 atol = eps
@test Healpix.pix2angNest(resol, 272282)[2] ≈ 6.1696707288760484 atol = eps
@test Healpix.pix2angNest(resol, 272487)[1] ≈ 1.9551931012905357 atol = eps
@test Healpix.pix2angNest(resol, 272487)[2] ≈ 6.0040008037843995 atol = eps
@test Healpix.pix2angNest(resol, 271668)[1] ≈ 1.9244782426347979 atol = eps
@test Healpix.pix2angNest(resol, 271668)[2] ≈ 6.2340979219672459 atol = eps
@test Healpix.pix2angNest(resol, 514561)[1] ≈ 1.2721038462777619 atol = eps
@test Healpix.pix2angNest(resol, 514561)[2] ≈ 4.7614763655970300 atol = eps
@test Healpix.pix2angNest(resol, 514766)[1] ≈ 1.2143369935354960 atol = eps
@test Healpix.pix2angNest(resol, 514766)[2] ≈ 4.7645443271728016 atol = eps
@test Healpix.pix2angNest(resol, 513946)[1] ≈ 1.1863995522992576 atol = eps
@test Healpix.pix2angNest(resol, 513946)[2] ≈ 4.9915734837798764 atol = eps
@test Healpix.pix2angNest(resol, 514151)[1] ≈ 1.2748275234359654 atol = eps
@test Healpix.pix2angNest(resol, 514151)[2] ≈ 4.8259035586882275 atol = eps
@test Healpix.pix2angNest(resol, 515380)[1] ≈ 1.1552980808319346 atol = eps
@test Healpix.pix2angNest(resol, 515380)[2] ≈ 4.9578259064463923 atol = eps
@test Healpix.pix2angNest(resol, 710773)[1] ≈ 1.9467798202990929 atol = eps
@test Healpix.pix2angNest(resol, 710773)[2] ≈ 4.0558452031696355 atol = eps
@test Healpix.pix2angNest(resol, 710928)[1] ≈ 1.9356052871934553 atol = eps
@test Healpix.pix2angNest(resol, 710928)[2] ≈ 4.0742529726242633 atol = eps
@test Healpix.pix2angNest(resol, 711083)[1] ≈ 1.9106332362490186 atol = eps
@test Healpix.pix2angNest(resol, 711083)[2] ≈ 4.0282335489876937 atol = eps
@test Healpix.pix2angNest(resol, 711238)[1] ≈ 1.9217036703859778 atol = eps
@test Healpix.pix2angNest(resol, 711238)[2] ≈ 4.0098257795330658 atol = eps
@test Healpix.pix2angNest(resol, 711393)[1] ≈ 1.8968556764797520 atol = eps
@test Healpix.pix2angNest(resol, 711393)[2] ≈ 3.9638063558964967 atol = eps
@test Healpix.pix2angNest(resol, 710773)[1] ≈ 1.9467798202990929 atol = eps
@test Healpix.pix2angNest(resol, 710773)[2] ≈ 4.0558452031696355 atol = eps
@test Healpix.pix2angNest(resol, 786429)[1] ≈ 1.5786089062699391 atol = eps
@test Healpix.pix2angNest(resol, 786429)[2] ≈ 5.4977871437821380 atol = eps
@test Healpix.pix2angNest(resol, 786429)[1] ≈ 1.5786089062699391 atol = eps
@test Healpix.pix2angNest(resol, 786429)[2] ≈ 5.4977871437821380 atol = eps
@test Healpix.pix2angNest(resol, 786430)[1] ≈ 1.5760046836760313 atol = eps
@test Healpix.pix2angNest(resol, 786430)[2] ≈ 5.5008551053579096 atol = eps
@test Healpix.pix2angNest(resol, 786431)[1] ≈ 1.5760046836760313 atol = eps
@test Healpix.pix2angNest(resol, 786431)[2] ≈ 5.4947191822063663 atol = eps
@test Healpix.pix2angNest(resol, 786432)[1] ≈ 1.5734004964050115 atol = eps
@test Healpix.pix2angNest(resol, 786432)[2] ≈ 5.4977871437821380 atol = eps

# pix2angRing

@test Healpix.pix2angRing(lowresol,  1)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angRing(lowresol,  1)[2] ≈ 0.785398163397448 atol = eps
@test Healpix.pix2angRing(lowresol,  2)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angRing(lowresol,  2)[2] ≈ 2.356194490192345 atol = eps
@test Healpix.pix2angRing(lowresol,  3)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angRing(lowresol,  3)[2] ≈ 3.926990816987241 atol = eps
@test Healpix.pix2angRing(lowresol,  4)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angRing(lowresol,  4)[2] ≈ 5.497787143782138 atol = eps
@test Healpix.pix2angRing(lowresol,  5)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol,  5)[2] ≈ 0.392699081698724 atol = eps
@test Healpix.pix2angRing(lowresol,  6)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol,  6)[2] ≈ 1.178097245096172 atol = eps
@test Healpix.pix2angRing(lowresol,  7)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol,  7)[2] ≈ 1.963495408493621 atol = eps
@test Healpix.pix2angRing(lowresol,  8)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol,  8)[2] ≈ 2.748893571891069 atol = eps
@test Healpix.pix2angRing(lowresol,  9)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol,  9)[2] ≈ 3.534291735288517 atol = eps
@test Healpix.pix2angRing(lowresol, 10)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol, 10)[2] ≈ 4.319689898685965 atol = eps
@test Healpix.pix2angRing(lowresol, 11)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol, 11)[2] ≈ 5.105088062083414 atol = eps
@test Healpix.pix2angRing(lowresol, 12)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol, 12)[2] ≈ 5.890486225480862 atol = eps
@test Healpix.pix2angRing(lowresol, 13)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angRing(lowresol, 13)[2] ≈ 0.000000000000000 atol = eps
@test Healpix.pix2angRing(lowresol, 14)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angRing(lowresol, 14)[2] ≈ 0.785398163397448 atol = eps
@test Healpix.pix2angRing(lowresol, 15)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angRing(lowresol, 15)[2] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angRing(lowresol, 16)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angRing(lowresol, 16)[2] ≈ 2.356194490192345 atol = eps
@test Healpix.pix2angRing(lowresol, 17)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angRing(lowresol, 17)[2] ≈ 3.141592653589793 atol = eps
@test Healpix.pix2angRing(lowresol, 18)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angRing(lowresol, 18)[2] ≈ 3.926990816987241 atol = eps
@test Healpix.pix2angRing(lowresol, 19)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angRing(lowresol, 19)[2] ≈ 4.712388980384690 atol = eps
@test Healpix.pix2angRing(lowresol, 20)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angRing(lowresol, 20)[2] ≈ 5.497787143782138 atol = eps
@test Healpix.pix2angRing(lowresol, 21)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angRing(lowresol, 21)[2] ≈ 0.392699081698724 atol = eps
@test Healpix.pix2angRing(lowresol, 22)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angRing(lowresol, 22)[2] ≈ 1.178097245096172 atol = eps
@test Healpix.pix2angRing(lowresol, 23)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angRing(lowresol, 23)[2] ≈ 1.963495408493621 atol = eps
@test Healpix.pix2angRing(lowresol, 24)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angRing(lowresol, 24)[2] ≈ 2.748893571891069 atol = eps
@test Healpix.pix2angRing(lowresol, 25)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angRing(lowresol, 25)[2] ≈ 3.534291735288517 atol = eps
@test Healpix.pix2angRing(lowresol, 26)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angRing(lowresol, 26)[2] ≈ 4.319689898685965 atol = eps
@test Healpix.pix2angRing(lowresol, 27)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angRing(lowresol, 27)[2] ≈ 5.105088062083414 atol = eps
@test Healpix.pix2angRing(lowresol, 28)[1] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angRing(lowresol, 28)[2] ≈ 5.890486225480862 atol = eps
@test Healpix.pix2angRing(lowresol, 29)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angRing(lowresol, 29)[2] ≈ 0.000000000000000 atol = eps
@test Healpix.pix2angRing(lowresol, 30)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angRing(lowresol, 30)[2] ≈ 0.785398163397448 atol = eps
@test Healpix.pix2angRing(lowresol, 31)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angRing(lowresol, 31)[2] ≈ 1.570796326794897 atol = eps
@test Healpix.pix2angRing(lowresol, 32)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angRing(lowresol, 32)[2] ≈ 2.356194490192345 atol = eps
@test Healpix.pix2angRing(lowresol, 33)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angRing(lowresol, 33)[2] ≈ 3.141592653589793 atol = eps
@test Healpix.pix2angRing(lowresol, 34)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angRing(lowresol, 34)[2] ≈ 3.926990816987241 atol = eps
@test Healpix.pix2angRing(lowresol, 35)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angRing(lowresol, 35)[2] ≈ 4.712388980384690 atol = eps
@test Healpix.pix2angRing(lowresol, 36)[1] ≈ 1.910633236249019 atol = eps
@test Healpix.pix2angRing(lowresol, 36)[2] ≈ 5.497787143782138 atol = eps
@test Healpix.pix2angRing(lowresol, 37)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angRing(lowresol, 37)[2] ≈ 0.392699081698724 atol = eps
@test Healpix.pix2angRing(lowresol, 38)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angRing(lowresol, 38)[2] ≈ 1.178097245096172 atol = eps
@test Healpix.pix2angRing(lowresol, 39)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angRing(lowresol, 39)[2] ≈ 1.963495408493621 atol = eps
@test Healpix.pix2angRing(lowresol, 40)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angRing(lowresol, 40)[2] ≈ 2.748893571891069 atol = eps
@test Healpix.pix2angRing(lowresol, 41)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angRing(lowresol, 41)[2] ≈ 3.534291735288517 atol = eps
@test Healpix.pix2angRing(lowresol, 42)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angRing(lowresol, 42)[2] ≈ 4.319689898685965 atol = eps
@test Healpix.pix2angRing(lowresol, 43)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angRing(lowresol, 43)[2] ≈ 5.105088062083414 atol = eps
@test Healpix.pix2angRing(lowresol, 44)[1] ≈ 2.300523983021863 atol = eps
@test Healpix.pix2angRing(lowresol, 44)[2] ≈ 5.890486225480862 atol = eps
@test Healpix.pix2angRing(lowresol, 45)[1] ≈ 2.730454791267445 atol = eps
@test Healpix.pix2angRing(lowresol, 45)[2] ≈ 0.785398163397448 atol = eps
@test Healpix.pix2angRing(lowresol, 46)[1] ≈ 2.730454791267445 atol = eps
@test Healpix.pix2angRing(lowresol, 46)[2] ≈ 2.356194490192345 atol = eps
@test Healpix.pix2angRing(lowresol, 47)[1] ≈ 2.730454791267445 atol = eps
@test Healpix.pix2angRing(lowresol, 47)[2] ≈ 3.926990816987241 atol = eps

@test Healpix.pix2angRing(resol,      1)[1] ≈ 0.0031894411211228 atol = eps
@test Healpix.pix2angRing(resol,      1)[2] ≈ 0.7853981633974483 atol = eps
@test Healpix.pix2angRing(resol,      2)[1] ≈ 0.0031894411211113 atol = eps
@test Healpix.pix2angRing(resol,      2)[2] ≈ 2.3561944901923448 atol = eps
@test Healpix.pix2angRing(resol,      3)[1] ≈ 0.0031894411211113 atol = eps
@test Healpix.pix2angRing(resol,      3)[2] ≈ 3.9269908169872414 atol = eps
@test Healpix.pix2angRing(resol,      4)[1] ≈ 0.0031894411211113 atol = eps
@test Healpix.pix2angRing(resol,      4)[2] ≈ 5.4977871437821380 atol = eps
@test Healpix.pix2angRing(resol,  74885)[1] ≈ 0.6290728186524177 atol = eps
@test Healpix.pix2angRing(resol,  74885)[2] ≈ 0.0040484441412240 atol = eps
@test Healpix.pix2angRing(resol,  75040)[1] ≈ 0.6290728186524177 atol = eps
@test Healpix.pix2angRing(resol,  75040)[2] ≈ 1.2590661279206516 atol = eps
@test Healpix.pix2angRing(resol,  75195)[1] ≈ 0.6290728186524177 atol = eps
@test Healpix.pix2angRing(resol,  75195)[2] ≈ 2.5140838117000794 atol = eps
@test Healpix.pix2angRing(resol,  75350)[1] ≈ 0.6290728186524177 atol = eps
@test Healpix.pix2angRing(resol,  75350)[2] ≈ 3.7691014954795072 atol = eps
@test Healpix.pix2angRing(resol,  75505)[1] ≈ 0.6290728186524177 atol = eps
@test Healpix.pix2angRing(resol,  75505)[2] ≈ 5.0241191792589346 atol = eps
@test Healpix.pix2angRing(resol, 270849)[1] ≈ 1.2557128565725952 atol = eps
@test Healpix.pix2angRing(resol, 270849)[2] ≈ 0.0000000000000000 atol = eps
@test Healpix.pix2angRing(resol, 271054)[1] ≈ 1.2557128565725952 atol = eps
@test Healpix.pix2angRing(resol, 271054)[2] ≈ 1.2578642460662257 atol = eps
@test Healpix.pix2angRing(resol, 272282)[1] ≈ 1.2584506449956214 atol = eps
@test Healpix.pix2angRing(resol, 272282)[2] ≈ 2.5126605305566803 atol = eps
@test Healpix.pix2angRing(resol, 272487)[1] ≈ 1.2584506449956214 atol = eps
@test Healpix.pix2angRing(resol, 272487)[2] ≈ 3.7705247766229060 atol = eps
@test Healpix.pix2angRing(resol, 271668)[1] ≈ 1.2557128565725952 atol = eps
@test Healpix.pix2angRing(resol, 271668)[2] ≈ 5.0253210611133605 atol = eps
@test Healpix.pix2angRing(resol, 514561)[1] ≈ 1.8858797970171981 atol = eps
@test Healpix.pix2angRing(resol, 514561)[2] ≈ 0.0000000000000000 atol = eps
@test Healpix.pix2angRing(resol, 514766)[1] ≈ 1.8858797970171981 atol = eps
@test Healpix.pix2angRing(resol, 514766)[2] ≈ 1.2578642460662257 atol = eps
@test Healpix.pix2angRing(resol, 513946)[1] ≈ 1.8831420085941717 atol = eps
@test Healpix.pix2angRing(resol, 513946)[2] ≈ 2.5126605305566803 atol = eps
@test Healpix.pix2angRing(resol, 514151)[1] ≈ 1.8831420085941717 atol = eps
@test Healpix.pix2angRing(resol, 514151)[2] ≈ 3.7705247766229060 atol = eps
@test Healpix.pix2angRing(resol, 515380)[1] ≈ 1.8858797970171981 atol = eps
@test Healpix.pix2angRing(resol, 515380)[2] ≈ 5.0253210611133605 atol = eps
@test Healpix.pix2angRing(resol, 710773)[1] ≈ 2.5125198349373754 atol = eps
@test Healpix.pix2angRing(resol, 710773)[2] ≈ 0.0040484441412240 atol = eps
@test Healpix.pix2angRing(resol, 710928)[1] ≈ 2.5125198349373754 atol = eps
@test Healpix.pix2angRing(resol, 710928)[2] ≈ 1.2590661279206516 atol = eps
@test Healpix.pix2angRing(resol, 711083)[1] ≈ 2.5125198349373754 atol = eps
@test Healpix.pix2angRing(resol, 711083)[2] ≈ 2.5140838117000794 atol = eps
@test Healpix.pix2angRing(resol, 711238)[1] ≈ 2.5125198349373754 atol = eps
@test Healpix.pix2angRing(resol, 711238)[2] ≈ 3.7691014954795072 atol = eps
@test Healpix.pix2angRing(resol, 711393)[1] ≈ 2.5125198349373754 atol = eps
@test Healpix.pix2angRing(resol, 711393)[2] ≈ 5.0241191792589346 atol = eps
@test Healpix.pix2angRing(resol, 710773)[1] ≈ 2.5125198349373754 atol = eps
@test Healpix.pix2angRing(resol, 710773)[2] ≈ 0.0040484441412240 atol = eps
@test Healpix.pix2angRing(resol, 786429)[1] ≈ 3.1384032124686820 atol = eps
@test Healpix.pix2angRing(resol, 786429)[2] ≈ 0.7853981633974483 atol = eps
@test Healpix.pix2angRing(resol, 786429)[1] ≈ 3.1384032124686820 atol = eps
@test Healpix.pix2angRing(resol, 786429)[2] ≈ 0.7853981633974483 atol = eps
@test Healpix.pix2angRing(resol, 786430)[1] ≈ 3.1384032124686820 atol = eps
@test Healpix.pix2angRing(resol, 786430)[2] ≈ 2.3561944901923448 atol = eps
@test Healpix.pix2angRing(resol, 786431)[1] ≈ 3.1384032124686820 atol = eps
@test Healpix.pix2angRing(resol, 786431)[2] ≈ 3.9269908169872414 atol = eps
@test Healpix.pix2angRing(resol, 786432)[1] ≈ 3.1384032124686820 atol = eps
@test Healpix.pix2angRing(resol, 786432)[2] ≈ 5.4977871437821380 atol = eps

include("interp.jl")
include("xyf.jl")

# Conformability

@test Healpix.conformables(Healpix.Map{Int16,Healpix.RingOrder}(4), 
                           Healpix.Map{Float32,Healpix.RingOrder}(4))
# nside mismatch
@test !Healpix.conformables(Healpix.Map{Int16,Healpix.RingOrder}(8),
                            Healpix.Map{Int16,Healpix.RingOrder}(4))
# order mismatch
@test !Healpix.conformables(Healpix.Map{Float32,Healpix.RingOrder}(4),
                            Healpix.Map{Float32,Healpix.NestedOrder}(4))

# Map making

pixidx1 = [1, 3, 2, 2, 3]
tod1 = [0.0, 2.5, 1.5, 2.5, 3.5]
(binmap1, hitmap1) = Healpix.tod2map(pixidx1, tod1, nside=1, ordering=Healpix.RingOrder)
@test binmap1.pixels ≈ [0.0, 2.0, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test hitmap1.pixels == [1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0]

binmap2 = Healpix.Map{Float64,Healpix.RingOrder}(1)
binmap2.pixels = [1.0, 0.0, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
hitmap2 = Healpix.Map{Int,Healpix.RingOrder}(1)
hitmap2.pixels = [1, 6, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0]
Healpix.combinemaps!(binmap1, hitmap1, binmap2, hitmap2)
@test binmap1.pixels ≈ [0.5, 0.5, 3.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
@test hitmap1.pixels == [2, 8, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0]

# Map loading

m = Healpix.readMapFromFITS("float_map.fits", 1, Float32)
@test typeof(m) == Healpix.Map{Float32,Healpix.RingOrder}
@test m.resolution.nside == 4
@test m.pixels == [Float32(x) for x in 0:(12 * 4^2 - 1)]

m = Healpix.readMapFromFITS("int_map.fits", 1, Int8)
@test typeof(m) == Healpix.Map{Int8,Healpix.RingOrder}
@test m.resolution.nside == 1
@test m.pixels == [Int8(x) for x in 0:11]

# Map saving

const mapFileName = tempname()
print("Saving $mapFileName\n")
Healpix.saveToFITS(m, "!$mapFileName", typechar="I")
m2 = Healpix.readMapFromFITS(mapFileName, 1, Int8)
@test m.pixels == m2.pixels

# Map projections

m = Healpix.Map{Float64,Healpix.RingOrder}(1)
m.pixels = 1.0:12.0

fig = Healpix.mollweide(m)
figname = tempname()
Cairo.write_to_png(fig, figname)
println("Mollweide projection saved in file $figname")

fig = Healpix.equirectangular(m)
figname = tempname()
Cairo.write_to_png(fig, figname)
println("Equirectangular projection saved in file $figname")

fig = Healpix.orthographic(m, 0.0, 0.0)
figname = tempname()
Cairo.write_to_png(fig, figname)
println("Orthographic projection saved in file $figname")

# Alm creation

@test Healpix.numberOfAlms(10, 5) == 51
@test Healpix.numberOfAlms(10, 7) == 60
@test Healpix.numberOfAlms(12, 7) == 76
@test Healpix.numberOfAlms(12, 12) == 91
@test_throws DomainError(-1, "`lmax` is not positive or zero") Healpix.numberOfAlms(-1, 1)
@test_throws DomainError(-1, "`mmax` is not positive or zero") Healpix.numberOfAlms(4, -1)
@test_throws DomainError((5, 7), "`lmax` and `mmax` are inconsistent") Healpix.numberOfAlms(5, 7)

alm = Healpix.Alm{ComplexF32}(10, 8)
@test Healpix.almIndex(alm, 4, 2) == 24
@test Healpix.almIndex(alm, 5, 2) == 25
@test Healpix.almIndex(alm, 5, 3) == 33
@test Healpix.almIndex(alm, [4, 6, 5], [3, 4, 5]) == [32, 41, 46]

alm = Healpix.readAlmFromFITS("alm.fits", ComplexF64)
@test alm[1]  ≈ (5.443205775735e+03 + 0.000000000000e+00im) atol = eps
@test alm[2]  ≈ (-3.143659646589e+03 + 0.000000000000e+00im) atol = eps
@test alm[3]  ≈ (-8.445976910202e-07 + 0.000000000000e+00im) atol = eps
@test alm[4]  ≈ (3.003475555079e-07 + 0.000000000000e+00im) atol = eps
@test alm[5]  ≈ (-1.094164444296e-06 + 0.000000000000e+00im) atol = eps
@test alm[6]  ≈ (3.745732939005e-07 + 0.000000000000e+00im) atol = eps
@test alm[7]  ≈ (-1.344818023454e-06 + 0.000000000000e+00im) atol = eps
@test alm[8]  ≈ (6.658742467775e-01 + -3.280201017201e+01im) atol = eps
@test alm[9]  ≈ (1.200156696497e-15 + -1.483545355539e-15im) atol = eps
@test alm[10] ≈ (-1.959362683381e-01 + -3.069662938925e+00im) atol = eps
@test alm[11] ≈ (-4.068874968688e-15 + -5.766800364761e-16im) atol = eps
@test alm[12] ≈ (-1.134265758838e-01 + 1.102692148708e+00im) atol = eps
@test alm[13] ≈ (-2.153997451326e-15 + -5.152134573109e-15im) atol = eps
@test alm[14] ≈ (-6.895464505703e-01 + 1.597235776935e+01im) atol = eps
@test alm[15] ≈ (4.344648310121e-15 + -1.892253965354e-15im) atol = eps
@test alm[16] ≈ (6.488192068017e-02 + 3.967268430232e+00im) atol = eps
@test alm[17] ≈ (1.340722726247e-15 + -9.822172672993e-15im) atol = eps
@test alm[18] ≈ (1.485579846001e-01 + 7.942515768792e-01im) atol = eps
@test alm[19] ≈ (6.938894269949e-01 + -1.028821633749e+01im) atol = eps
@test alm[20] ≈ (3.412745093853e-15 + -4.420258002639e-15im) atol = eps
@test alm[21] ≈ (3.712279444471e-02 + -3.564529788432e+00im) atol = eps
@test alm[22] ≈ (-6.205890749143e-15 + 4.855843625358e-16im) atol = eps
@test alm[23] ≈ (-6.893793878861e-01 + 7.453657075068e+00im) atol = eps
@test alm[24] ≈ (4.644726834988e-11 + -4.685776289892e-15im) atol = eps
@test alm[25] ≈ (-1.162363120645e-01 + 3.084870414802e+00im) atol = eps
@test alm[26] ≈ (6.806207052149e-01 + -5.769903744427e+00im) atol = eps
@test alm[27] ≈ (-1.307082489340e-15 + 1.143824238899e-14im) atol = eps
@test alm[28] ≈ (-6.698490836781e-01 + 4.661675665246e+00im) atol = eps
