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
@test_throws DomainError(2 * Healpix.NSIDE_MAX) Healpix.Resolution(2 * Healpix.NSIDE_MAX)

lowresol = Healpix.Resolution(2)
resol = Healpix.Resolution(256)
highresol = Healpix.Resolution(2^29)

@test Healpix.ang2pixNest(resol, 0.0000000000000000, 0.0000000000000000) == 65536
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 1.2566370614359172) == 65536
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 2.5132741228718345) == 131072
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 3.7699111843077517) == 196608
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 5.0265482457436690) == 262144
@test Healpix.ang2pixNest(resol, 0.0000000000000000, 6.2831853071795862) == 262144

@test Healpix.ang2pixNest(resol, 0.6283185307179586, 0.0000000000000000) == 45055
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 1.2566370614359172) == 31074
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 2.5132741228718345) == 116111
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 3.7699111843077517) == 182862
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 5.0265482457436690) == 243347
@test Healpix.ang2pixNest(resol, 0.6283185307179586, 6.2831853071795862) == 221182

@test Healpix.ang2pixNest(resol, 1.2566370614359172, 0.0000000000000000) == 315344
@test Healpix.ang2pixNest(resol, 1.2566370614359172, 1.2566370614359172) == 387305
@test Healpix.ang2pixNest(resol, 1.2566370614359172, 2.5132741228718345) == 71955
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


@test Healpix.ang2pixNest(highresol, 1.570796325553133199, 0.785398163397448279) == 1
@test Healpix.ang2pixNest(highresol, 1.570796324311369840, 0.785398164860366288) == 2
@test Healpix.ang2pixNest(highresol, 1.570796324311369840, 0.785398161934530159) == 3
@test Healpix.ang2pixNest(highresol, 1.570796323069606260, 0.785398163397448279) == 4
@test Healpix.ang2pixNest(highresol, 1.570796323069606260, 0.785398166323284408) == 5
@test Healpix.ang2pixNest(highresol, 1.570796321827842901, 0.785398167786202528) == 6
@test Healpix.ang2pixNest(highresol, 1.570796321827842901, 0.785398164860366288) == 7
@test Healpix.ang2pixNest(highresol, 1.570796320586079542, 0.785398166323284408) == 8
@test Healpix.ang2pixNest(highresol, 1.570796323069606260, 0.785398160471612039) == 9
@test Healpix.ang2pixNest(highresol, 1.570796321827842901, 0.785398161934530159) == 10
@test Healpix.ang2pixNest(highresol, 1.570796325553133199, 5.497787143782137953) == 864691128455135233
@test Healpix.ang2pixNest(highresol, 1.570796324311369840, 5.497787145245055740) == 864691128455135234
@test Healpix.ang2pixNest(highresol, 1.570796324311369840, 5.497787142319219278) == 864691128455135235
@test Healpix.ang2pixNest(highresol, 1.570796323069606260, 5.497787143782137953) == 864691128455135236
@test Healpix.ang2pixNest(highresol, 1.570796323069606260, 5.497787146707974415) == 864691128455135237
@test Healpix.ang2pixNest(highresol, 1.570796321827842901, 5.497787148170891314) == 864691128455135238
@test Healpix.ang2pixNest(highresol, 1.570796321827842901, 5.497787145245055740) == 864691128455135239
@test Healpix.ang2pixNest(highresol, 1.570796320586079542, 5.497787146707974415) == 864691128455135240
@test Healpix.ang2pixNest(highresol, 1.570796323069606260, 5.497787140856301491) == 864691128455135241
@test Healpix.ang2pixNest(highresol, 1.570796321827842901, 5.497787142319219278) == 864691128455135242
@test Healpix.ang2pixNest(highresol, 2.300523981355862446, 3.141592653589793116) == 1729382256910270465
@test Healpix.ang2pixNest(highresol, 2.300523979689862220, 3.141592655052710903) == 1729382256910270466
@test Healpix.ang2pixNest(highresol, 2.300523979689862220, 3.141592652126874885) == 1729382256910270467
@test Healpix.ang2pixNest(highresol, 2.300523978023861549, 3.141592653589793116) == 1729382256910270468
@test Healpix.ang2pixNest(highresol, 2.300523978023861549, 3.141592656515629134) == 1729382256910270469
@test Healpix.ang2pixNest(highresol, 2.300523976357860878, 3.141592657978547365) == 1729382256910270470
@test Healpix.ang2pixNest(highresol, 2.300523976357860878, 3.141592655052710903) == 1729382256910270471
@test Healpix.ang2pixNest(highresol, 2.300523974691860651, 3.141592656515629134) == 1729382256910270472
@test Healpix.ang2pixNest(highresol, 2.300523978023861549, 3.141592650663956654) == 1729382256910270473
@test Healpix.ang2pixNest(highresol, 2.300523976357860878, 3.141592652126874885) == 1729382256910270474
@test Healpix.ang2pixNest(highresol, 3.141592652068949665, 2.356194490192344837) == 2594073385365405697
@test Healpix.ang2pixNest(highresol, 3.141592650548106658, 2.748893571891068976) == 2594073385365405698
@test Healpix.ang2pixNest(highresol, 3.141592650548106658, 1.963495408493620697) == 2594073385365405699
@test Healpix.ang2pixNest(highresol, 3.141592649027263207, 2.356194490192344837) == 2594073385365405700
@test Healpix.ang2pixNest(highresol, 3.141592649027263207, 2.879793265790643542) == 2594073385365405701
@test Healpix.ang2pixNest(highresol, 3.141592647506419755, 2.945243112740430824) == 2594073385365405702
@test Healpix.ang2pixNest(highresol, 3.141592647506419755, 2.552544031041707129) == 2594073385365405703
@test Healpix.ang2pixNest(highresol, 3.141592645985576304, 2.670353755551324149) == 2594073385365405704
@test Healpix.ang2pixNest(highresol, 3.141592649027263207, 1.832595714594045910) == 2594073385365405705
@test Healpix.ang2pixNest(highresol, 3.141592647506419755, 2.159844949342982545) == 2594073385365405706
@test Healpix.ang2pixNest(highresol, 1.570796331761950437, 5.497787148170891314) == 3458764513820540918
@test Healpix.ang2pixNest(highresol, 1.570796331761950437, 5.497787145245055740) == 3458764513820540919
@test Healpix.ang2pixNest(highresol, 1.570796330520186856, 5.497787146707974415) == 3458764513820540920
@test Healpix.ang2pixNest(highresol, 1.570796333003713796, 5.497787140856301491) == 3458764513820540921
@test Healpix.ang2pixNest(highresol, 1.570796331761950437, 5.497787142319219278) == 3458764513820540922
@test Healpix.ang2pixNest(highresol, 1.570796331761950437, 5.497787139393382816) == 3458764513820540923
@test Healpix.ang2pixNest(highresol, 1.570796330520186856, 5.497787140856301491) == 3458764513820540924
@test Healpix.ang2pixNest(highresol, 1.570796330520186856, 5.497787143782137953) == 3458764513820540925
@test Healpix.ang2pixNest(highresol, 1.570796329278423498, 5.497787145245055740) == 3458764513820540926
@test Healpix.ang2pixNest(highresol, 1.570796329278423498, 5.497787142319219278) == 3458764513820540927

# ang2pixRing

@test Healpix.ang2pixRing(resol, 0.0000000000000000, 0.0000000000000000) == 1
@test Healpix.ang2pixRing(resol, 0.0000000000000000, 1.2566370614359172) == 1
@test Healpix.ang2pixRing(resol, 0.0000000000000000, 2.5132741228718345) == 2
@test Healpix.ang2pixRing(resol, 0.0000000000000000, 3.7699111843077517) == 3
@test Healpix.ang2pixRing(resol, 0.0000000000000000, 5.0265482457436690) == 4
@test Healpix.ang2pixRing(resol, 0.0000000000000000, 6.2831853071795862) == 1
@test Healpix.ang2pixRing(resol, 0.6283185307179586, 0.0000000000000000) == 74885
@test Healpix.ang2pixRing(resol, 0.6283185307179586, 1.2566370614359172) == 75040
@test Healpix.ang2pixRing(resol, 0.6283185307179586, 2.5132741228718345) == 75195
@test Healpix.ang2pixRing(resol, 0.6283185307179586, 3.7699111843077517) == 75350
@test Healpix.ang2pixRing(resol, 0.6283185307179586, 5.0265482457436690) == 75505
@test Healpix.ang2pixRing(resol, 0.6283185307179586, 6.2831853071795862) == 74885
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

@test Healpix.ang2pixRing(highresol, 0.000000001520843396, 0.785398163397448279) == 1
@test Healpix.ang2pixRing(highresol, 0.000000001520843396, 2.356194490192344837) == 2
@test Healpix.ang2pixRing(highresol, 0.000000001520843396, 3.926990816987241395) == 3
@test Healpix.ang2pixRing(highresol, 0.000000001520843396, 5.497787143782137953) == 4
@test Healpix.ang2pixRing(highresol, 0.000000003041686792, 0.392699081698724139) == 5
@test Healpix.ang2pixRing(highresol, 0.000000003041686792, 1.178097245096172418) == 6
@test Healpix.ang2pixRing(highresol, 0.000000003041686792, 1.963495408493620697) == 7
@test Healpix.ang2pixRing(highresol, 0.000000003041686792, 2.748893571891068976) == 8
@test Healpix.ang2pixRing(highresol, 0.000000003041686792, 3.534291735288517255) == 9
@test Healpix.ang2pixRing(highresol, 0.000000003041686792, 4.319689898685965090) == 10
@test Healpix.ang2pixRing(highresol, 1.047197551196597853, 3.141592655052710903) == 864691128455135233
@test Healpix.ang2pixRing(highresol, 1.047197551196597853, 3.141592657978547365) == 864691128455135234
@test Healpix.ang2pixRing(highresol, 1.047197551196597853, 3.141592660904383383) == 864691128455135235
@test Healpix.ang2pixRing(highresol, 1.047197551196597853, 3.141592663830219401) == 864691128455135236
@test Healpix.ang2pixRing(highresol, 1.047197551196597853, 3.141592666756055863) == 864691128455135237
@test Healpix.ang2pixRing(highresol, 1.047197551196597853, 3.141592669681891881) == 864691128455135238
@test Healpix.ang2pixRing(highresol, 1.047197551196597853, 3.141592672607727899) == 864691128455135239
@test Healpix.ang2pixRing(highresol, 1.047197551196597853, 3.141592675533564361) == 864691128455135240
@test Healpix.ang2pixRing(highresol, 1.047197551196597853, 3.141592678459400378) == 864691128455135241
@test Healpix.ang2pixRing(highresol, 1.047197551196597853, 3.141592681385236396) == 864691128455135242
@test Healpix.ang2pixRing(highresol, 1.570796326794896558, 3.141592655052710903) == 1729382256910270465
@test Healpix.ang2pixRing(highresol, 1.570796326794896558, 3.141592657978547365) == 1729382256910270466
@test Healpix.ang2pixRing(highresol, 1.570796326794896558, 3.141592660904383383) == 1729382256910270467
@test Healpix.ang2pixRing(highresol, 1.570796326794896558, 3.141592663830219401) == 1729382256910270468
@test Healpix.ang2pixRing(highresol, 1.570796326794896558, 3.141592666756055863) == 1729382256910270469
@test Healpix.ang2pixRing(highresol, 1.570796326794896558, 3.141592669681891881) == 1729382256910270470
@test Healpix.ang2pixRing(highresol, 1.570796326794896558, 3.141592672607727899) == 1729382256910270471
@test Healpix.ang2pixRing(highresol, 1.570796326794896558, 3.141592675533564361) == 1729382256910270472
@test Healpix.ang2pixRing(highresol, 1.570796326794896558, 3.141592678459400378) == 1729382256910270473
@test Healpix.ang2pixRing(highresol, 1.570796326794896558, 3.141592681385236396) == 1729382256910270474
@test Healpix.ang2pixRing(highresol, 2.094395102393195707, 3.141592655052710903) == 2594073385365405697
@test Healpix.ang2pixRing(highresol, 2.094395102393195707, 3.141592657978547365) == 2594073385365405698
@test Healpix.ang2pixRing(highresol, 2.094395102393195707, 3.141592660904383383) == 2594073385365405699
@test Healpix.ang2pixRing(highresol, 2.094395102393195707, 3.141592663830219401) == 2594073385365405700
@test Healpix.ang2pixRing(highresol, 2.094395102393195707, 3.141592666756055863) == 2594073385365405701
@test Healpix.ang2pixRing(highresol, 2.094395102393195707, 3.141592669681891881) == 2594073385365405702
@test Healpix.ang2pixRing(highresol, 2.094395102393195707, 3.141592672607727899) == 2594073385365405703
@test Healpix.ang2pixRing(highresol, 2.094395102393195707, 3.141592675533564361) == 2594073385365405704
@test Healpix.ang2pixRing(highresol, 2.094395102393195707, 3.141592678459400378) == 2594073385365405705
@test Healpix.ang2pixRing(highresol, 2.094395102393195707, 3.141592681385236396) == 2594073385365405706
@test Healpix.ang2pixRing(highresol, 3.141592650548106658, 1.178097245096172418) == 3458764513820540918
@test Healpix.ang2pixRing(highresol, 3.141592650548106658, 1.963495408493620697) == 3458764513820540919
@test Healpix.ang2pixRing(highresol, 3.141592650548106658, 2.748893571891068976) == 3458764513820540920
@test Healpix.ang2pixRing(highresol, 3.141592650548106658, 3.534291735288517255) == 3458764513820540921
@test Healpix.ang2pixRing(highresol, 3.141592650548106658, 4.319689898685965090) == 3458764513820540922
@test Healpix.ang2pixRing(highresol, 3.141592650548106658, 5.105088062083414258) == 3458764513820540923
@test Healpix.ang2pixRing(highresol, 3.141592650548106658, 5.890486225480861648) == 3458764513820540924
@test Healpix.ang2pixRing(highresol, 3.141592652068949665, 0.785398163397448279) == 3458764513820540925
@test Healpix.ang2pixRing(highresol, 3.141592652068949665, 2.356194490192344837) == 3458764513820540926
@test Healpix.ang2pixRing(highresol, 3.141592652068949665, 3.926990816987241395) == 3458764513820540927

# pix2angNest

@test Healpix.pix2angNest(lowresol, 1)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angNest(lowresol, 1)[2] ≈ 0.785398163397448 atol = eps
@test Healpix.pix2angNest(lowresol, 2)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angNest(lowresol, 2)[2] ≈ 1.178097245096172 atol = eps
@test Healpix.pix2angNest(lowresol, 3)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angNest(lowresol, 3)[2] ≈ 0.392699081698724 atol = eps
@test Healpix.pix2angNest(lowresol, 4)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angNest(lowresol, 4)[2] ≈ 0.785398163397448 atol = eps
@test Healpix.pix2angNest(lowresol, 5)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angNest(lowresol, 5)[2] ≈ 2.356194490192345 atol = eps
@test Healpix.pix2angNest(lowresol, 6)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angNest(lowresol, 6)[2] ≈ 2.748893571891069 atol = eps
@test Healpix.pix2angNest(lowresol, 7)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angNest(lowresol, 7)[2] ≈ 1.963495408493621 atol = eps
@test Healpix.pix2angNest(lowresol, 8)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angNest(lowresol, 8)[2] ≈ 2.356194490192345 atol = eps
@test Healpix.pix2angNest(lowresol, 9)[1] ≈ 1.230959417340775 atol = eps
@test Healpix.pix2angNest(lowresol, 9)[2] ≈ 3.926990816987241 atol = eps
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

@test Healpix.pix2angNest(resol, 1)[1] ≈ 1.5681921571847817 atol = eps
@test Healpix.pix2angNest(resol, 1)[2] ≈ 0.7853981633974483 atol = eps
@test Healpix.pix2angNest(resol, 2)[1] ≈ 1.5655879699137618 atol = eps
@test Healpix.pix2angNest(resol, 2)[2] ≈ 0.7884661249732196 atol = eps
@test Healpix.pix2angNest(resol, 3)[1] ≈ 1.5655879699137618 atol = eps
@test Healpix.pix2angNest(resol, 3)[2] ≈ 0.7823302018216770 atol = eps
@test Healpix.pix2angNest(resol, 4)[1] ≈ 1.5629837473198540 atol = eps
@test Healpix.pix2angNest(resol, 4)[2] ≈ 0.7853981633974483 atol = eps
@test Healpix.pix2angNest(resol, 74885)[1] ≈ 1.2884125182769750 atol = eps
@test Healpix.pix2angNest(resol, 74885)[2] ≈ 2.2396119503130363 atol = eps
@test Healpix.pix2angNest(resol, 75040)[1] ≈ 1.2447369771100412 atol = eps
@test Healpix.pix2angNest(resol, 75040)[2] ≈ 2.3193789512830896 atol = eps
@test Healpix.pix2angNest(resol, 75195)[1] ≈ 1.2198889832038156 atol = eps
@test Healpix.pix2angNest(resol, 75195)[2] ≈ 2.2733595276465204 atol = eps
@test Healpix.pix2angNest(resol, 75350)[1] ≈ 1.2309594173407747 atol = eps
@test Healpix.pix2angNest(resol, 75350)[2] ≈ 2.2549517581918925 atol = eps
@test Healpix.pix2angNest(resol, 75505)[1] ≈ 1.2059873663963379 atol = eps
@test Healpix.pix2angNest(resol, 75505)[2] ≈ 2.2089323345553233 atol = eps
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

@test Healpix.pix2angNest(highresol, 1)[1] ≈ 1.570796325553133199 atol = eps
@test Healpix.pix2angNest(highresol, 1)[2] ≈ 0.785398163397448279 atol = eps
@test Healpix.pix2angNest(highresol, 2)[1] ≈ 1.570796324311369840 atol = eps
@test Healpix.pix2angNest(highresol, 2)[2] ≈ 0.785398164860366288 atol = eps
@test Healpix.pix2angNest(highresol, 3)[1] ≈ 1.570796324311369840 atol = eps
@test Healpix.pix2angNest(highresol, 3)[2] ≈ 0.785398161934530159 atol = eps
@test Healpix.pix2angNest(highresol, 4)[1] ≈ 1.570796323069606260 atol = eps
@test Healpix.pix2angNest(highresol, 4)[2] ≈ 0.785398163397448279 atol = eps
@test Healpix.pix2angNest(highresol, 5)[1] ≈ 1.570796323069606260 atol = eps
@test Healpix.pix2angNest(highresol, 5)[2] ≈ 0.785398166323284408 atol = eps
@test Healpix.pix2angNest(highresol, 6)[1] ≈ 1.570796321827842901 atol = eps
@test Healpix.pix2angNest(highresol, 6)[2] ≈ 0.785398167786202528 atol = eps
@test Healpix.pix2angNest(highresol, 7)[1] ≈ 1.570796321827842901 atol = eps
@test Healpix.pix2angNest(highresol, 7)[2] ≈ 0.785398164860366288 atol = eps
@test Healpix.pix2angNest(highresol, 8)[1] ≈ 1.570796320586079542 atol = eps
@test Healpix.pix2angNest(highresol, 8)[2] ≈ 0.785398166323284408 atol = eps
@test Healpix.pix2angNest(highresol, 9)[1] ≈ 1.570796323069606260 atol = eps
@test Healpix.pix2angNest(highresol, 9)[2] ≈ 0.785398160471612039 atol = eps
@test Healpix.pix2angNest(highresol, 10)[1] ≈ 1.570796321827842901 atol = eps
@test Healpix.pix2angNest(highresol, 10)[2] ≈ 0.785398161934530159 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135233)[1] ≈ 1.570796325553133199 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135233)[2] ≈ 5.497787143782137953 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135234)[1] ≈ 1.570796324311369840 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135234)[2] ≈ 5.497787145245055740 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135235)[1] ≈ 1.570796324311369840 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135235)[2] ≈ 5.497787142319219278 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135236)[1] ≈ 1.570796323069606260 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135236)[2] ≈ 5.497787143782137953 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135237)[1] ≈ 1.570796323069606260 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135237)[2] ≈ 5.497787146707974415 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135238)[1] ≈ 1.570796321827842901 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135238)[2] ≈ 5.497787148170891314 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135239)[1] ≈ 1.570796321827842901 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135239)[2] ≈ 5.497787145245055740 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135240)[1] ≈ 1.570796320586079542 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135240)[2] ≈ 5.497787146707974415 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135241)[1] ≈ 1.570796323069606260 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135241)[2] ≈ 5.497787140856301491 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135242)[1] ≈ 1.570796321827842901 atol = eps
@test Healpix.pix2angNest(highresol, 864691128455135242)[2] ≈ 5.497787142319219278 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270465)[1] ≈ 2.300523981355862446 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270465)[2] ≈ 3.141592653589793116 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270466)[1] ≈ 2.300523979689862220 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270466)[2] ≈ 3.141592655052710903 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270467)[1] ≈ 2.300523979689862220 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270467)[2] ≈ 3.141592652126874885 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270468)[1] ≈ 2.300523978023861549 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270468)[2] ≈ 3.141592653589793116 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270469)[1] ≈ 2.300523978023861549 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270469)[2] ≈ 3.141592656515629134 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270470)[1] ≈ 2.300523976357860878 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270470)[2] ≈ 3.141592657978547365 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270471)[1] ≈ 2.300523976357860878 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270471)[2] ≈ 3.141592655052710903 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270472)[1] ≈ 2.300523974691860651 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270472)[2] ≈ 3.141592656515629134 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270473)[1] ≈ 2.300523978023861549 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270473)[2] ≈ 3.141592650663956654 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270474)[1] ≈ 2.300523976357860878 atol = eps
@test Healpix.pix2angNest(highresol, 1729382256910270474)[2] ≈ 3.141592652126874885 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405697)[1] ≈ 3.141592652068949665 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405697)[2] ≈ 2.356194490192344837 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405698)[1] ≈ 3.141592650548106658 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405698)[2] ≈ 2.748893571891068976 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405699)[1] ≈ 3.141592650548106658 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405699)[2] ≈ 1.963495408493620697 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405700)[1] ≈ 3.141592649027263207 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405700)[2] ≈ 2.356194490192344837 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405701)[1] ≈ 3.141592649027263207 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405701)[2] ≈ 2.879793265790643542 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405702)[1] ≈ 3.141592647506419755 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405702)[2] ≈ 2.945243112740430824 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405703)[1] ≈ 3.141592647506419755 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405703)[2] ≈ 2.552544031041707129 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405704)[1] ≈ 3.141592645985576304 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405704)[2] ≈ 2.670353755551324149 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405705)[1] ≈ 3.141592649027263207 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405705)[2] ≈ 1.832595714594045910 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405706)[1] ≈ 3.141592647506419755 atol = eps
@test Healpix.pix2angNest(highresol, 2594073385365405706)[2] ≈ 2.159844949342982545 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540918)[1] ≈ 1.570796331761950437 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540918)[2] ≈ 5.497787148170891314 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540919)[1] ≈ 1.570796331761950437 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540919)[2] ≈ 5.497787145245055740 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540920)[1] ≈ 1.570796330520186856 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540920)[2] ≈ 5.497787146707974415 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540921)[1] ≈ 1.570796333003713796 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540921)[2] ≈ 5.497787140856301491 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540922)[1] ≈ 1.570796331761950437 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540922)[2] ≈ 5.497787142319219278 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540923)[1] ≈ 1.570796331761950437 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540923)[2] ≈ 5.497787139393382816 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540924)[1] ≈ 1.570796330520186856 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540924)[2] ≈ 5.497787140856301491 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540925)[1] ≈ 1.570796330520186856 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540925)[2] ≈ 5.497787143782137953 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540926)[1] ≈ 1.570796329278423498 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540926)[2] ≈ 5.497787145245055740 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540927)[1] ≈ 1.570796329278423498 atol = eps
@test Healpix.pix2angNest(highresol, 3458764513820540927)[2] ≈ 5.497787142319219278 atol = eps

# pix2angRing

@test Healpix.pix2angRing(lowresol, 1)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angRing(lowresol, 1)[2] ≈ 0.785398163397448 atol = eps
@test Healpix.pix2angRing(lowresol, 2)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angRing(lowresol, 2)[2] ≈ 2.356194490192345 atol = eps
@test Healpix.pix2angRing(lowresol, 3)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angRing(lowresol, 3)[2] ≈ 3.926990816987241 atol = eps
@test Healpix.pix2angRing(lowresol, 4)[1] ≈ 0.411137862322348 atol = eps
@test Healpix.pix2angRing(lowresol, 4)[2] ≈ 5.497787143782138 atol = eps
@test Healpix.pix2angRing(lowresol, 5)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol, 5)[2] ≈ 0.392699081698724 atol = eps
@test Healpix.pix2angRing(lowresol, 6)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol, 6)[2] ≈ 1.178097245096172 atol = eps
@test Healpix.pix2angRing(lowresol, 7)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol, 7)[2] ≈ 1.963495408493621 atol = eps
@test Healpix.pix2angRing(lowresol, 8)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol, 8)[2] ≈ 2.748893571891069 atol = eps
@test Healpix.pix2angRing(lowresol, 9)[1] ≈ 0.841068670567930 atol = eps
@test Healpix.pix2angRing(lowresol, 9)[2] ≈ 3.534291735288517 atol = eps
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

@test Healpix.pix2angRing(resol, 1)[1] ≈ 0.0031894411211228 atol = eps
@test Healpix.pix2angRing(resol, 1)[2] ≈ 0.7853981633974483 atol = eps
@test Healpix.pix2angRing(resol, 2)[1] ≈ 0.0031894411211113 atol = eps
@test Healpix.pix2angRing(resol, 2)[2] ≈ 2.3561944901923448 atol = eps
@test Healpix.pix2angRing(resol, 3)[1] ≈ 0.0031894411211113 atol = eps
@test Healpix.pix2angRing(resol, 3)[2] ≈ 3.9269908169872414 atol = eps
@test Healpix.pix2angRing(resol, 4)[1] ≈ 0.0031894411211113 atol = eps
@test Healpix.pix2angRing(resol, 4)[2] ≈ 5.4977871437821380 atol = eps
@test Healpix.pix2angRing(resol, 74885)[1] ≈ 0.6290728186524177 atol = eps
@test Healpix.pix2angRing(resol, 74885)[2] ≈ 0.0040484441412240 atol = eps
@test Healpix.pix2angRing(resol, 75040)[1] ≈ 0.6290728186524177 atol = eps
@test Healpix.pix2angRing(resol, 75040)[2] ≈ 1.2590661279206516 atol = eps
@test Healpix.pix2angRing(resol, 75195)[1] ≈ 0.6290728186524177 atol = eps
@test Healpix.pix2angRing(resol, 75195)[2] ≈ 2.5140838117000794 atol = eps
@test Healpix.pix2angRing(resol, 75350)[1] ≈ 0.6290728186524177 atol = eps
@test Healpix.pix2angRing(resol, 75350)[2] ≈ 3.7691014954795072 atol = eps
@test Healpix.pix2angRing(resol, 75505)[1] ≈ 0.6290728186524177 atol = eps
@test Healpix.pix2angRing(resol, 75505)[2] ≈ 5.0241191792589346 atol = eps
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

@test Healpix.pix2angRing(highresol, 1)[1] ≈ 0.000000001520843396 atol = eps
@test Healpix.pix2angRing(highresol, 1)[2] ≈ 0.785398163397448279 atol = eps
@test Healpix.pix2angRing(highresol, 2)[1] ≈ 0.000000001520843396 atol = eps
@test Healpix.pix2angRing(highresol, 2)[2] ≈ 2.356194490192344837 atol = eps
@test Healpix.pix2angRing(highresol, 3)[1] ≈ 0.000000001520843396 atol = eps
@test Healpix.pix2angRing(highresol, 3)[2] ≈ 3.926990816987241395 atol = eps
@test Healpix.pix2angRing(highresol, 4)[1] ≈ 0.000000001520843396 atol = eps
@test Healpix.pix2angRing(highresol, 4)[2] ≈ 5.497787143782137953 atol = eps
@test Healpix.pix2angRing(highresol, 5)[1] ≈ 0.000000003041686792 atol = eps
@test Healpix.pix2angRing(highresol, 5)[2] ≈ 0.392699081698724139 atol = eps
@test Healpix.pix2angRing(highresol, 6)[1] ≈ 0.000000003041686792 atol = eps
@test Healpix.pix2angRing(highresol, 6)[2] ≈ 1.178097245096172418 atol = eps
@test Healpix.pix2angRing(highresol, 7)[1] ≈ 0.000000003041686792 atol = eps
@test Healpix.pix2angRing(highresol, 7)[2] ≈ 1.963495408493620697 atol = eps
@test Healpix.pix2angRing(highresol, 8)[1] ≈ 0.000000003041686792 atol = eps
@test Healpix.pix2angRing(highresol, 8)[2] ≈ 2.748893571891068976 atol = eps
@test Healpix.pix2angRing(highresol, 9)[1] ≈ 0.000000003041686792 atol = eps
@test Healpix.pix2angRing(highresol, 9)[2] ≈ 3.534291735288517255 atol = eps
@test Healpix.pix2angRing(highresol, 10)[1] ≈ 0.000000003041686792 atol = eps
@test Healpix.pix2angRing(highresol, 10)[2] ≈ 4.319689898685965090 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135233)[1] ≈ 1.047197551196597853 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135233)[2] ≈ 3.141592655052710903 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135234)[1] ≈ 1.047197551196597853 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135234)[2] ≈ 3.141592657978547365 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135235)[1] ≈ 1.047197551196597853 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135235)[2] ≈ 3.141592660904383383 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135236)[1] ≈ 1.047197551196597853 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135236)[2] ≈ 3.141592663830219401 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135237)[1] ≈ 1.047197551196597853 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135237)[2] ≈ 3.141592666756055863 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135238)[1] ≈ 1.047197551196597853 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135238)[2] ≈ 3.141592669681891881 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135239)[1] ≈ 1.047197551196597853 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135239)[2] ≈ 3.141592672607727899 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135240)[1] ≈ 1.047197551196597853 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135240)[2] ≈ 3.141592675533564361 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135241)[1] ≈ 1.047197551196597853 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135241)[2] ≈ 3.141592678459400378 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135242)[1] ≈ 1.047197551196597853 atol = eps
@test Healpix.pix2angRing(highresol, 864691128455135242)[2] ≈ 3.141592681385236396 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270465)[1] ≈ 1.570796326794896558 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270465)[2] ≈ 3.141592655052710903 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270466)[1] ≈ 1.570796326794896558 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270466)[2] ≈ 3.141592657978547365 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270467)[1] ≈ 1.570796326794896558 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270467)[2] ≈ 3.141592660904383383 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270468)[1] ≈ 1.570796326794896558 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270468)[2] ≈ 3.141592663830219401 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270469)[1] ≈ 1.570796326794896558 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270469)[2] ≈ 3.141592666756055863 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270470)[1] ≈ 1.570796326794896558 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270470)[2] ≈ 3.141592669681891881 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270471)[1] ≈ 1.570796326794896558 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270471)[2] ≈ 3.141592672607727899 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270472)[1] ≈ 1.570796326794896558 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270472)[2] ≈ 3.141592675533564361 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270473)[1] ≈ 1.570796326794896558 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270473)[2] ≈ 3.141592678459400378 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270474)[1] ≈ 1.570796326794896558 atol = eps
@test Healpix.pix2angRing(highresol, 1729382256910270474)[2] ≈ 3.141592681385236396 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405697)[1] ≈ 2.094395102393195707 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405697)[2] ≈ 3.141592655052710903 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405698)[1] ≈ 2.094395102393195707 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405698)[2] ≈ 3.141592657978547365 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405699)[1] ≈ 2.094395102393195707 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405699)[2] ≈ 3.141592660904383383 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405700)[1] ≈ 2.094395102393195707 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405700)[2] ≈ 3.141592663830219401 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405701)[1] ≈ 2.094395102393195707 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405701)[2] ≈ 3.141592666756055863 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405702)[1] ≈ 2.094395102393195707 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405702)[2] ≈ 3.141592669681891881 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405703)[1] ≈ 2.094395102393195707 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405703)[2] ≈ 3.141592672607727899 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405704)[1] ≈ 2.094395102393195707 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405704)[2] ≈ 3.141592675533564361 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405705)[1] ≈ 2.094395102393195707 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405705)[2] ≈ 3.141592678459400378 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405706)[1] ≈ 2.094395102393195707 atol = eps
@test Healpix.pix2angRing(highresol, 2594073385365405706)[2] ≈ 3.141592681385236396 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540918)[1] ≈ 3.141592650548106658 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540918)[2] ≈ 1.178097245096172418 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540919)[1] ≈ 3.141592650548106658 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540919)[2] ≈ 1.963495408493620697 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540920)[1] ≈ 3.141592650548106658 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540920)[2] ≈ 2.748893571891068976 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540921)[1] ≈ 3.141592650548106658 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540921)[2] ≈ 3.534291735288517255 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540922)[1] ≈ 3.141592650548106658 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540922)[2] ≈ 4.319689898685965090 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540923)[1] ≈ 3.141592650548106658 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540923)[2] ≈ 5.105088062083414258 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540924)[1] ≈ 3.141592650548106658 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540924)[2] ≈ 5.890486225480861648 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540925)[1] ≈ 3.141592652068949665 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540925)[2] ≈ 0.785398163397448279 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540926)[1] ≈ 3.141592652068949665 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540926)[2] ≈ 2.356194490192344837 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540927)[1] ≈ 3.141592652068949665 atol = eps
@test Healpix.pix2angRing(highresol, 3458764513820540927)[2] ≈ 3.926990816987241395 atol = eps

@testset "pix2zphiRing, NSIDE=1" begin
    res = Healpix.Resolution(1)

    (z, phi) = Healpix.pix2zphiRing(res, 1)
    @test z ≈ 0.6666666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 2)
    @test z ≈ 0.6666666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 3)
    @test z ≈ 0.6666666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 4)
    @test z ≈ 0.6666666666667
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiRing(res, 5)
    @test z ≈ 0
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiRing(res, 6)
    @test z ≈ 0
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiRing(res, 7)
    @test z ≈ 0
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiRing(res, 8)
    @test z ≈ 0
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiRing(res, 9)
    @test z ≈ -0.6666666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 10)
    @test z ≈ -0.6666666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 11)
    @test z ≈ -0.6666666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 12)
    @test z ≈ -0.6666666666667
    @test phi ≈ 5.497787143782
end

@testset "pix2zphiNest, NSIDE=1" begin
    res = Healpix.Resolution(1)

    (z, phi) = Healpix.pix2zphiNest(res, 1)
    @test z ≈ 0.6666666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 2)
    @test z ≈ 0.6666666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 3)
    @test z ≈ 0.6666666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 4)
    @test z ≈ 0.6666666666667
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiNest(res, 5)
    @test z ≈ 0
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiNest(res, 6)
    @test z ≈ 0
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiNest(res, 7)
    @test z ≈ 0
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiNest(res, 8)
    @test z ≈ 0
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiNest(res, 9)
    @test z ≈ -0.6666666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 10)
    @test z ≈ -0.6666666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 11)
    @test z ≈ -0.6666666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 12)
    @test z ≈ -0.6666666666667
    @test phi ≈ 5.497787143782
end

@testset "Healpix.pix2zphiRing, NSIDE=2" begin
    res = Healpix.Resolution(2)

    (z, phi) = Healpix.pix2zphiRing(res, 1)
    @test z ≈ 0.9166666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 2)
    @test z ≈ 0.9166666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 3)
    @test z ≈ 0.9166666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 4)
    @test z ≈ 0.9166666666667
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiRing(res, 5)
    @test z ≈ 0.6666666666667
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 6)
    @test z ≈ 0.6666666666667
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiRing(res, 7)
    @test z ≈ 0.6666666666667
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 8)
    @test z ≈ 0.6666666666667
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiRing(res, 9)
    @test z ≈ 0.6666666666667
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiRing(res, 10)
    @test z ≈ 0.6666666666667
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiRing(res, 11)
    @test z ≈ 0.6666666666667
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiRing(res, 12)
    @test z ≈ 0.6666666666667
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 13)
    @test z ≈ 0.3333333333333
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiRing(res, 14)
    @test z ≈ 0.3333333333333
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 15)
    @test z ≈ 0.3333333333333
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiRing(res, 16)
    @test z ≈ 0.3333333333333
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 17)
    @test z ≈ 0.3333333333333
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiRing(res, 18)
    @test z ≈ 0.3333333333333
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 19)
    @test z ≈ 0.3333333333333
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiRing(res, 20)
    @test z ≈ 0.3333333333333
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiRing(res, 21)
    @test z ≈ 0
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 22)
    @test z ≈ 0
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiRing(res, 23)
    @test z ≈ 0
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 24)
    @test z ≈ 0
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiRing(res, 25)
    @test z ≈ 0
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiRing(res, 26)
    @test z ≈ 0
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiRing(res, 27)
    @test z ≈ 0
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiRing(res, 28)
    @test z ≈ 0
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 29)
    @test z ≈ -0.3333333333333
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiRing(res, 30)
    @test z ≈ -0.3333333333333
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 31)
    @test z ≈ -0.3333333333333
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiRing(res, 32)
    @test z ≈ -0.3333333333333
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 33)
    @test z ≈ -0.3333333333333
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiRing(res, 34)
    @test z ≈ -0.3333333333333
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 35)
    @test z ≈ -0.3333333333333
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiRing(res, 36)
    @test z ≈ -0.3333333333333
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiRing(res, 37)
    @test z ≈ -0.6666666666667
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 38)
    @test z ≈ -0.6666666666667
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiRing(res, 39)
    @test z ≈ -0.6666666666667
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 40)
    @test z ≈ -0.6666666666667
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiRing(res, 41)
    @test z ≈ -0.6666666666667
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiRing(res, 42)
    @test z ≈ -0.6666666666667
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiRing(res, 43)
    @test z ≈ -0.6666666666667
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiRing(res, 44)
    @test z ≈ -0.6666666666667
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 45)
    @test z ≈ -0.9166666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 46)
    @test z ≈ -0.9166666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 47)
    @test z ≈ -0.9166666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 48)
    @test z ≈ -0.9166666666667
    @test phi ≈ 5.497787143782
end

@testset "Healpix.pix2zphiNest, NSIDE=2" begin
    res = Healpix.Resolution(2)

    (z, phi) = Healpix.pix2zphiNest(res, 1)
    @test z ≈ 0.3333333333333
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 2)
    @test z ≈ 0.6666666666667
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiNest(res, 3)
    @test z ≈ 0.6666666666667
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 4)
    @test z ≈ 0.9166666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 5)
    @test z ≈ 0.3333333333333
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 6)
    @test z ≈ 0.6666666666667
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiNest(res, 7)
    @test z ≈ 0.6666666666667
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 8)
    @test z ≈ 0.9166666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 9)
    @test z ≈ 0.3333333333333
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 10)
    @test z ≈ 0.6666666666667
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiNest(res, 11)
    @test z ≈ 0.6666666666667
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiNest(res, 12)
    @test z ≈ 0.9166666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 13)
    @test z ≈ 0.3333333333333
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiNest(res, 14)
    @test z ≈ 0.6666666666667
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 15)
    @test z ≈ 0.6666666666667
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiNest(res, 16)
    @test z ≈ 0.9166666666667
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiNest(res, 17)
    @test z ≈ -0.3333333333333
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiNest(res, 18)
    @test z ≈ 0
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 19)
    @test z ≈ 0
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 20)
    @test z ≈ 0.3333333333333
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiNest(res, 21)
    @test z ≈ -0.3333333333333
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiNest(res, 22)
    @test z ≈ 0
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 23)
    @test z ≈ 0
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiNest(res, 24)
    @test z ≈ 0.3333333333333
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiNest(res, 25)
    @test z ≈ -0.3333333333333
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiNest(res, 26)
    @test z ≈ 0
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiNest(res, 27)
    @test z ≈ 0
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiNest(res, 28)
    @test z ≈ 0.3333333333333
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiNest(res, 29)
    @test z ≈ -0.3333333333333
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiNest(res, 30)
    @test z ≈ 0
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiNest(res, 31)
    @test z ≈ 0
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiNest(res, 32)
    @test z ≈ 0.3333333333333
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiNest(res, 33)
    @test z ≈ -0.9166666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 34)
    @test z ≈ -0.6666666666667
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiNest(res, 35)
    @test z ≈ -0.6666666666667
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 36)
    @test z ≈ -0.3333333333333
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 37)
    @test z ≈ -0.9166666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 38)
    @test z ≈ -0.6666666666667
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiNest(res, 39)
    @test z ≈ -0.6666666666667
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 40)
    @test z ≈ -0.3333333333333
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 41)
    @test z ≈ -0.9166666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 42)
    @test z ≈ -0.6666666666667
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiNest(res, 43)
    @test z ≈ -0.6666666666667
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiNest(res, 44)
    @test z ≈ -0.3333333333333
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 45)
    @test z ≈ -0.9166666666667
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiNest(res, 46)
    @test z ≈ -0.6666666666667
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 47)
    @test z ≈ -0.6666666666667
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiNest(res, 48)
    @test z ≈ -0.3333333333333
    @test phi ≈ 5.497787143782
end

@testset "Healpix.pix2zphiRing, NSIDE=4" begin
    res = Healpix.Resolution(4)

    (z, phi) = Healpix.pix2zphiRing(res, 1)
    @test z ≈ 0.9791666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 2)
    @test z ≈ 0.9791666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 3)
    @test z ≈ 0.9791666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 4)
    @test z ≈ 0.9791666666667
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiRing(res, 5)
    @test z ≈ 0.9166666666667
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 6)
    @test z ≈ 0.9166666666667
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiRing(res, 7)
    @test z ≈ 0.9166666666667
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 8)
    @test z ≈ 0.9166666666667
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiRing(res, 9)
    @test z ≈ 0.9166666666667
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiRing(res, 10)
    @test z ≈ 0.9166666666667
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiRing(res, 11)
    @test z ≈ 0.9166666666667
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiRing(res, 12)
    @test z ≈ 0.9166666666667
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 13)
    @test z ≈ 0.8125
    @test phi ≈ 0.2617993877991
    (z, phi) = Healpix.pix2zphiRing(res, 14)
    @test z ≈ 0.8125
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 15)
    @test z ≈ 0.8125
    @test phi ≈ 1.308996938996
    (z, phi) = Healpix.pix2zphiRing(res, 16)
    @test z ≈ 0.8125
    @test phi ≈ 1.832595714594
    (z, phi) = Healpix.pix2zphiRing(res, 17)
    @test z ≈ 0.8125
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 18)
    @test z ≈ 0.8125
    @test phi ≈ 2.879793265791
    (z, phi) = Healpix.pix2zphiRing(res, 19)
    @test z ≈ 0.8125
    @test phi ≈ 3.403392041389
    (z, phi) = Healpix.pix2zphiRing(res, 20)
    @test z ≈ 0.8125
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 21)
    @test z ≈ 0.8125
    @test phi ≈ 4.450589592586
    (z, phi) = Healpix.pix2zphiRing(res, 22)
    @test z ≈ 0.8125
    @test phi ≈ 4.974188368184
    (z, phi) = Healpix.pix2zphiRing(res, 23)
    @test z ≈ 0.8125
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiRing(res, 24)
    @test z ≈ 0.8125
    @test phi ≈ 6.02138591938
    (z, phi) = Healpix.pix2zphiRing(res, 25)
    @test z ≈ 0.6666666666667
    @test phi ≈ 0.1963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 26)
    @test z ≈ 0.6666666666667
    @test phi ≈ 0.5890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 27)
    @test z ≈ 0.6666666666667
    @test phi ≈ 0.9817477042468
    (z, phi) = Healpix.pix2zphiRing(res, 28)
    @test z ≈ 0.6666666666667
    @test phi ≈ 1.374446785946
    (z, phi) = Healpix.pix2zphiRing(res, 29)
    @test z ≈ 0.6666666666667
    @test phi ≈ 1.767145867644
    (z, phi) = Healpix.pix2zphiRing(res, 30)
    @test z ≈ 0.6666666666667
    @test phi ≈ 2.159844949343
    (z, phi) = Healpix.pix2zphiRing(res, 31)
    @test z ≈ 0.6666666666667
    @test phi ≈ 2.552544031042
    (z, phi) = Healpix.pix2zphiRing(res, 32)
    @test z ≈ 0.6666666666667
    @test phi ≈ 2.94524311274
    (z, phi) = Healpix.pix2zphiRing(res, 33)
    @test z ≈ 0.6666666666667
    @test phi ≈ 3.337942194439
    (z, phi) = Healpix.pix2zphiRing(res, 34)
    @test z ≈ 0.6666666666667
    @test phi ≈ 3.730641276138
    (z, phi) = Healpix.pix2zphiRing(res, 35)
    @test z ≈ 0.6666666666667
    @test phi ≈ 4.123340357837
    (z, phi) = Healpix.pix2zphiRing(res, 36)
    @test z ≈ 0.6666666666667
    @test phi ≈ 4.516039439535
    (z, phi) = Healpix.pix2zphiRing(res, 37)
    @test z ≈ 0.6666666666667
    @test phi ≈ 4.908738521234
    (z, phi) = Healpix.pix2zphiRing(res, 38)
    @test z ≈ 0.6666666666667
    @test phi ≈ 5.301437602933
    (z, phi) = Healpix.pix2zphiRing(res, 39)
    @test z ≈ 0.6666666666667
    @test phi ≈ 5.694136684632
    (z, phi) = Healpix.pix2zphiRing(res, 40)
    @test z ≈ 0.6666666666667
    @test phi ≈ 6.08683576633
    (z, phi) = Healpix.pix2zphiRing(res, 41)
    @test z ≈ 0.5
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiRing(res, 42)
    @test z ≈ 0.5
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 43)
    @test z ≈ 0.5
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 44)
    @test z ≈ 0.5
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiRing(res, 45)
    @test z ≈ 0.5
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiRing(res, 46)
    @test z ≈ 0.5
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 47)
    @test z ≈ 0.5
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 48)
    @test z ≈ 0.5
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiRing(res, 49)
    @test z ≈ 0.5
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiRing(res, 50)
    @test z ≈ 0.5
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiRing(res, 51)
    @test z ≈ 0.5
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 52)
    @test z ≈ 0.5
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiRing(res, 53)
    @test z ≈ 0.5
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiRing(res, 54)
    @test z ≈ 0.5
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiRing(res, 55)
    @test z ≈ 0.5
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiRing(res, 56)
    @test z ≈ 0.5
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 57)
    @test z ≈ 0.3333333333333
    @test phi ≈ 0.1963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 58)
    @test z ≈ 0.3333333333333
    @test phi ≈ 0.5890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 59)
    @test z ≈ 0.3333333333333
    @test phi ≈ 0.9817477042468
    (z, phi) = Healpix.pix2zphiRing(res, 60)
    @test z ≈ 0.3333333333333
    @test phi ≈ 1.374446785946
    (z, phi) = Healpix.pix2zphiRing(res, 61)
    @test z ≈ 0.3333333333333
    @test phi ≈ 1.767145867644
    (z, phi) = Healpix.pix2zphiRing(res, 62)
    @test z ≈ 0.3333333333333
    @test phi ≈ 2.159844949343
    (z, phi) = Healpix.pix2zphiRing(res, 63)
    @test z ≈ 0.3333333333333
    @test phi ≈ 2.552544031042
    (z, phi) = Healpix.pix2zphiRing(res, 64)
    @test z ≈ 0.3333333333333
    @test phi ≈ 2.94524311274
    (z, phi) = Healpix.pix2zphiRing(res, 65)
    @test z ≈ 0.3333333333333
    @test phi ≈ 3.337942194439
    (z, phi) = Healpix.pix2zphiRing(res, 66)
    @test z ≈ 0.3333333333333
    @test phi ≈ 3.730641276138
    (z, phi) = Healpix.pix2zphiRing(res, 67)
    @test z ≈ 0.3333333333333
    @test phi ≈ 4.123340357837
    (z, phi) = Healpix.pix2zphiRing(res, 68)
    @test z ≈ 0.3333333333333
    @test phi ≈ 4.516039439535
    (z, phi) = Healpix.pix2zphiRing(res, 69)
    @test z ≈ 0.3333333333333
    @test phi ≈ 4.908738521234
    (z, phi) = Healpix.pix2zphiRing(res, 70)
    @test z ≈ 0.3333333333333
    @test phi ≈ 5.301437602933
    (z, phi) = Healpix.pix2zphiRing(res, 71)
    @test z ≈ 0.3333333333333
    @test phi ≈ 5.694136684632
    (z, phi) = Healpix.pix2zphiRing(res, 72)
    @test z ≈ 0.3333333333333
    @test phi ≈ 6.08683576633
    (z, phi) = Healpix.pix2zphiRing(res, 73)
    @test z ≈ 0.1666666666667
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiRing(res, 74)
    @test z ≈ 0.1666666666667
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 75)
    @test z ≈ 0.1666666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 76)
    @test z ≈ 0.1666666666667
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiRing(res, 77)
    @test z ≈ 0.1666666666667
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiRing(res, 78)
    @test z ≈ 0.1666666666667
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 79)
    @test z ≈ 0.1666666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 80)
    @test z ≈ 0.1666666666667
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiRing(res, 81)
    @test z ≈ 0.1666666666667
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiRing(res, 82)
    @test z ≈ 0.1666666666667
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiRing(res, 83)
    @test z ≈ 0.1666666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 84)
    @test z ≈ 0.1666666666667
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiRing(res, 85)
    @test z ≈ 0.1666666666667
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiRing(res, 86)
    @test z ≈ 0.1666666666667
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiRing(res, 87)
    @test z ≈ 0.1666666666667
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiRing(res, 88)
    @test z ≈ 0.1666666666667
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 89)
    @test z ≈ 0
    @test phi ≈ 0.1963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 90)
    @test z ≈ 0
    @test phi ≈ 0.5890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 91)
    @test z ≈ 0
    @test phi ≈ 0.9817477042468
    (z, phi) = Healpix.pix2zphiRing(res, 92)
    @test z ≈ 0
    @test phi ≈ 1.374446785946
    (z, phi) = Healpix.pix2zphiRing(res, 93)
    @test z ≈ 0
    @test phi ≈ 1.767145867644
    (z, phi) = Healpix.pix2zphiRing(res, 94)
    @test z ≈ 0
    @test phi ≈ 2.159844949343
    (z, phi) = Healpix.pix2zphiRing(res, 95)
    @test z ≈ 0
    @test phi ≈ 2.552544031042
    (z, phi) = Healpix.pix2zphiRing(res, 96)
    @test z ≈ 0
    @test phi ≈ 2.94524311274
    (z, phi) = Healpix.pix2zphiRing(res, 97)
    @test z ≈ 0
    @test phi ≈ 3.337942194439
    (z, phi) = Healpix.pix2zphiRing(res, 98)
    @test z ≈ 0
    @test phi ≈ 3.730641276138
    (z, phi) = Healpix.pix2zphiRing(res, 99)
    @test z ≈ 0
    @test phi ≈ 4.123340357837
    (z, phi) = Healpix.pix2zphiRing(res, 100)
    @test z ≈ 0
    @test phi ≈ 4.516039439535
    (z, phi) = Healpix.pix2zphiRing(res, 101)
    @test z ≈ 0
    @test phi ≈ 4.908738521234
    (z, phi) = Healpix.pix2zphiRing(res, 102)
    @test z ≈ 0
    @test phi ≈ 5.301437602933
    (z, phi) = Healpix.pix2zphiRing(res, 103)
    @test z ≈ 0
    @test phi ≈ 5.694136684632
    (z, phi) = Healpix.pix2zphiRing(res, 104)
    @test z ≈ 0
    @test phi ≈ 6.08683576633
    (z, phi) = Healpix.pix2zphiRing(res, 105)
    @test z ≈ -0.1666666666667
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiRing(res, 106)
    @test z ≈ -0.1666666666667
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 107)
    @test z ≈ -0.1666666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 108)
    @test z ≈ -0.1666666666667
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiRing(res, 109)
    @test z ≈ -0.1666666666667
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiRing(res, 110)
    @test z ≈ -0.1666666666667
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 111)
    @test z ≈ -0.1666666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 112)
    @test z ≈ -0.1666666666667
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiRing(res, 113)
    @test z ≈ -0.1666666666667
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiRing(res, 114)
    @test z ≈ -0.1666666666667
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiRing(res, 115)
    @test z ≈ -0.1666666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 116)
    @test z ≈ -0.1666666666667
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiRing(res, 117)
    @test z ≈ -0.1666666666667
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiRing(res, 118)
    @test z ≈ -0.1666666666667
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiRing(res, 119)
    @test z ≈ -0.1666666666667
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiRing(res, 120)
    @test z ≈ -0.1666666666667
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 121)
    @test z ≈ -0.3333333333333
    @test phi ≈ 0.1963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 122)
    @test z ≈ -0.3333333333333
    @test phi ≈ 0.5890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 123)
    @test z ≈ -0.3333333333333
    @test phi ≈ 0.9817477042468
    (z, phi) = Healpix.pix2zphiRing(res, 124)
    @test z ≈ -0.3333333333333
    @test phi ≈ 1.374446785946
    (z, phi) = Healpix.pix2zphiRing(res, 125)
    @test z ≈ -0.3333333333333
    @test phi ≈ 1.767145867644
    (z, phi) = Healpix.pix2zphiRing(res, 126)
    @test z ≈ -0.3333333333333
    @test phi ≈ 2.159844949343
    (z, phi) = Healpix.pix2zphiRing(res, 127)
    @test z ≈ -0.3333333333333
    @test phi ≈ 2.552544031042
    (z, phi) = Healpix.pix2zphiRing(res, 128)
    @test z ≈ -0.3333333333333
    @test phi ≈ 2.94524311274
    (z, phi) = Healpix.pix2zphiRing(res, 129)
    @test z ≈ -0.3333333333333
    @test phi ≈ 3.337942194439
    (z, phi) = Healpix.pix2zphiRing(res, 130)
    @test z ≈ -0.3333333333333
    @test phi ≈ 3.730641276138
    (z, phi) = Healpix.pix2zphiRing(res, 131)
    @test z ≈ -0.3333333333333
    @test phi ≈ 4.123340357837
    (z, phi) = Healpix.pix2zphiRing(res, 132)
    @test z ≈ -0.3333333333333
    @test phi ≈ 4.516039439535
    (z, phi) = Healpix.pix2zphiRing(res, 133)
    @test z ≈ -0.3333333333333
    @test phi ≈ 4.908738521234
    (z, phi) = Healpix.pix2zphiRing(res, 134)
    @test z ≈ -0.3333333333333
    @test phi ≈ 5.301437602933
    (z, phi) = Healpix.pix2zphiRing(res, 135)
    @test z ≈ -0.3333333333333
    @test phi ≈ 5.694136684632
    (z, phi) = Healpix.pix2zphiRing(res, 136)
    @test z ≈ -0.3333333333333
    @test phi ≈ 6.08683576633
    (z, phi) = Healpix.pix2zphiRing(res, 137)
    @test z ≈ -0.5
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiRing(res, 138)
    @test z ≈ -0.5
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 139)
    @test z ≈ -0.5
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 140)
    @test z ≈ -0.5
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiRing(res, 141)
    @test z ≈ -0.5
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiRing(res, 142)
    @test z ≈ -0.5
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 143)
    @test z ≈ -0.5
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 144)
    @test z ≈ -0.5
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiRing(res, 145)
    @test z ≈ -0.5
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiRing(res, 146)
    @test z ≈ -0.5
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiRing(res, 147)
    @test z ≈ -0.5
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 148)
    @test z ≈ -0.5
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiRing(res, 149)
    @test z ≈ -0.5
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiRing(res, 150)
    @test z ≈ -0.5
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiRing(res, 151)
    @test z ≈ -0.5
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiRing(res, 152)
    @test z ≈ -0.5
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 153)
    @test z ≈ -0.6666666666667
    @test phi ≈ 0.1963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 154)
    @test z ≈ -0.6666666666667
    @test phi ≈ 0.5890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 155)
    @test z ≈ -0.6666666666667
    @test phi ≈ 0.9817477042468
    (z, phi) = Healpix.pix2zphiRing(res, 156)
    @test z ≈ -0.6666666666667
    @test phi ≈ 1.374446785946
    (z, phi) = Healpix.pix2zphiRing(res, 157)
    @test z ≈ -0.6666666666667
    @test phi ≈ 1.767145867644
    (z, phi) = Healpix.pix2zphiRing(res, 158)
    @test z ≈ -0.6666666666667
    @test phi ≈ 2.159844949343
    (z, phi) = Healpix.pix2zphiRing(res, 159)
    @test z ≈ -0.6666666666667
    @test phi ≈ 2.552544031042
    (z, phi) = Healpix.pix2zphiRing(res, 160)
    @test z ≈ -0.6666666666667
    @test phi ≈ 2.94524311274
    (z, phi) = Healpix.pix2zphiRing(res, 161)
    @test z ≈ -0.6666666666667
    @test phi ≈ 3.337942194439
    (z, phi) = Healpix.pix2zphiRing(res, 162)
    @test z ≈ -0.6666666666667
    @test phi ≈ 3.730641276138
    (z, phi) = Healpix.pix2zphiRing(res, 163)
    @test z ≈ -0.6666666666667
    @test phi ≈ 4.123340357837
    (z, phi) = Healpix.pix2zphiRing(res, 164)
    @test z ≈ -0.6666666666667
    @test phi ≈ 4.516039439535
    (z, phi) = Healpix.pix2zphiRing(res, 165)
    @test z ≈ -0.6666666666667
    @test phi ≈ 4.908738521234
    (z, phi) = Healpix.pix2zphiRing(res, 166)
    @test z ≈ -0.6666666666667
    @test phi ≈ 5.301437602933
    (z, phi) = Healpix.pix2zphiRing(res, 167)
    @test z ≈ -0.6666666666667
    @test phi ≈ 5.694136684632
    (z, phi) = Healpix.pix2zphiRing(res, 168)
    @test z ≈ -0.6666666666667
    @test phi ≈ 6.08683576633
    (z, phi) = Healpix.pix2zphiRing(res, 169)
    @test z ≈ -0.8125
    @test phi ≈ 0.2617993877991
    (z, phi) = Healpix.pix2zphiRing(res, 170)
    @test z ≈ -0.8125
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 171)
    @test z ≈ -0.8125
    @test phi ≈ 1.308996938996
    (z, phi) = Healpix.pix2zphiRing(res, 172)
    @test z ≈ -0.8125
    @test phi ≈ 1.832595714594
    (z, phi) = Healpix.pix2zphiRing(res, 173)
    @test z ≈ -0.8125
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 174)
    @test z ≈ -0.8125
    @test phi ≈ 2.879793265791
    (z, phi) = Healpix.pix2zphiRing(res, 175)
    @test z ≈ -0.8125
    @test phi ≈ 3.403392041389
    (z, phi) = Healpix.pix2zphiRing(res, 176)
    @test z ≈ -0.8125
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 177)
    @test z ≈ -0.8125
    @test phi ≈ 4.450589592586
    (z, phi) = Healpix.pix2zphiRing(res, 178)
    @test z ≈ -0.8125
    @test phi ≈ 4.974188368184
    (z, phi) = Healpix.pix2zphiRing(res, 179)
    @test z ≈ -0.8125
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiRing(res, 180)
    @test z ≈ -0.8125
    @test phi ≈ 6.02138591938
    (z, phi) = Healpix.pix2zphiRing(res, 181)
    @test z ≈ -0.9166666666667
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 182)
    @test z ≈ -0.9166666666667
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiRing(res, 183)
    @test z ≈ -0.9166666666667
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiRing(res, 184)
    @test z ≈ -0.9166666666667
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiRing(res, 185)
    @test z ≈ -0.9166666666667
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiRing(res, 186)
    @test z ≈ -0.9166666666667
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiRing(res, 187)
    @test z ≈ -0.9166666666667
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiRing(res, 188)
    @test z ≈ -0.9166666666667
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiRing(res, 189)
    @test z ≈ -0.9791666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiRing(res, 190)
    @test z ≈ -0.9791666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiRing(res, 191)
    @test z ≈ -0.9791666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiRing(res, 192)
    @test z ≈ -0.9791666666667
    @test phi ≈ 5.497787143782
end

@testset "Healpix.pix2zphiNest, NSIDE=4" begin
    res = Healpix.Resolution(4)

    (z, phi) = Healpix.pix2zphiNest(res, 1)
    @test z ≈ 0.1666666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 2)
    @test z ≈ 0.3333333333333
    @test phi ≈ 0.9817477042468
    (z, phi) = Healpix.pix2zphiNest(res, 3)
    @test z ≈ 0.3333333333333
    @test phi ≈ 0.5890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 4)
    @test z ≈ 0.5
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 5)
    @test z ≈ 0.5
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiNest(res, 6)
    @test z ≈ 0.6666666666667
    @test phi ≈ 1.374446785946
    (z, phi) = Healpix.pix2zphiNest(res, 7)
    @test z ≈ 0.6666666666667
    @test phi ≈ 0.9817477042468
    (z, phi) = Healpix.pix2zphiNest(res, 8)
    @test z ≈ 0.8125
    @test phi ≈ 1.308996938996
    (z, phi) = Healpix.pix2zphiNest(res, 9)
    @test z ≈ 0.5
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 10)
    @test z ≈ 0.6666666666667
    @test phi ≈ 0.5890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 11)
    @test z ≈ 0.6666666666667
    @test phi ≈ 0.1963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 12)
    @test z ≈ 0.8125
    @test phi ≈ 0.2617993877991
    (z, phi) = Healpix.pix2zphiNest(res, 13)
    @test z ≈ 0.8125
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 14)
    @test z ≈ 0.9166666666667
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiNest(res, 15)
    @test z ≈ 0.9166666666667
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 16)
    @test z ≈ 0.9791666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 17)
    @test z ≈ 0.1666666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 18)
    @test z ≈ 0.3333333333333
    @test phi ≈ 2.552544031042
    (z, phi) = Healpix.pix2zphiNest(res, 19)
    @test z ≈ 0.3333333333333
    @test phi ≈ 2.159844949343
    (z, phi) = Healpix.pix2zphiNest(res, 20)
    @test z ≈ 0.5
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 21)
    @test z ≈ 0.5
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiNest(res, 22)
    @test z ≈ 0.6666666666667
    @test phi ≈ 2.94524311274
    (z, phi) = Healpix.pix2zphiNest(res, 23)
    @test z ≈ 0.6666666666667
    @test phi ≈ 2.552544031042
    (z, phi) = Healpix.pix2zphiNest(res, 24)
    @test z ≈ 0.8125
    @test phi ≈ 2.879793265791
    (z, phi) = Healpix.pix2zphiNest(res, 25)
    @test z ≈ 0.5
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 26)
    @test z ≈ 0.6666666666667
    @test phi ≈ 2.159844949343
    (z, phi) = Healpix.pix2zphiNest(res, 27)
    @test z ≈ 0.6666666666667
    @test phi ≈ 1.767145867644
    (z, phi) = Healpix.pix2zphiNest(res, 28)
    @test z ≈ 0.8125
    @test phi ≈ 1.832595714594
    (z, phi) = Healpix.pix2zphiNest(res, 29)
    @test z ≈ 0.8125
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 30)
    @test z ≈ 0.9166666666667
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiNest(res, 31)
    @test z ≈ 0.9166666666667
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 32)
    @test z ≈ 0.9791666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 33)
    @test z ≈ 0.1666666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 34)
    @test z ≈ 0.3333333333333
    @test phi ≈ 4.123340357837
    (z, phi) = Healpix.pix2zphiNest(res, 35)
    @test z ≈ 0.3333333333333
    @test phi ≈ 3.730641276138
    (z, phi) = Healpix.pix2zphiNest(res, 36)
    @test z ≈ 0.5
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 37)
    @test z ≈ 0.5
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiNest(res, 38)
    @test z ≈ 0.6666666666667
    @test phi ≈ 4.516039439535
    (z, phi) = Healpix.pix2zphiNest(res, 39)
    @test z ≈ 0.6666666666667
    @test phi ≈ 4.123340357837
    (z, phi) = Healpix.pix2zphiNest(res, 40)
    @test z ≈ 0.8125
    @test phi ≈ 4.450589592586
    (z, phi) = Healpix.pix2zphiNest(res, 41)
    @test z ≈ 0.5
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiNest(res, 42)
    @test z ≈ 0.6666666666667
    @test phi ≈ 3.730641276138
    (z, phi) = Healpix.pix2zphiNest(res, 43)
    @test z ≈ 0.6666666666667
    @test phi ≈ 3.337942194439
    (z, phi) = Healpix.pix2zphiNest(res, 44)
    @test z ≈ 0.8125
    @test phi ≈ 3.403392041389
    (z, phi) = Healpix.pix2zphiNest(res, 45)
    @test z ≈ 0.8125
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 46)
    @test z ≈ 0.9166666666667
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiNest(res, 47)
    @test z ≈ 0.9166666666667
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiNest(res, 48)
    @test z ≈ 0.9791666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 49)
    @test z ≈ 0.1666666666667
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiNest(res, 50)
    @test z ≈ 0.3333333333333
    @test phi ≈ 5.694136684632
    (z, phi) = Healpix.pix2zphiNest(res, 51)
    @test z ≈ 0.3333333333333
    @test phi ≈ 5.301437602933
    (z, phi) = Healpix.pix2zphiNest(res, 52)
    @test z ≈ 0.5
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiNest(res, 53)
    @test z ≈ 0.5
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 54)
    @test z ≈ 0.6666666666667
    @test phi ≈ 6.08683576633
    (z, phi) = Healpix.pix2zphiNest(res, 55)
    @test z ≈ 0.6666666666667
    @test phi ≈ 5.694136684632
    (z, phi) = Healpix.pix2zphiNest(res, 56)
    @test z ≈ 0.8125
    @test phi ≈ 6.02138591938
    (z, phi) = Healpix.pix2zphiNest(res, 57)
    @test z ≈ 0.5
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiNest(res, 58)
    @test z ≈ 0.6666666666667
    @test phi ≈ 5.301437602933
    (z, phi) = Healpix.pix2zphiNest(res, 59)
    @test z ≈ 0.6666666666667
    @test phi ≈ 4.908738521234
    (z, phi) = Healpix.pix2zphiNest(res, 60)
    @test z ≈ 0.8125
    @test phi ≈ 4.974188368184
    (z, phi) = Healpix.pix2zphiNest(res, 61)
    @test z ≈ 0.8125
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiNest(res, 62)
    @test z ≈ 0.9166666666667
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 63)
    @test z ≈ 0.9166666666667
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiNest(res, 64)
    @test z ≈ 0.9791666666667
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiNest(res, 65)
    @test z ≈ -0.5
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiNest(res, 66)
    @test z ≈ -0.3333333333333
    @test phi ≈ 0.1963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 67)
    @test z ≈ -0.3333333333333
    @test phi ≈ 6.08683576633
    (z, phi) = Healpix.pix2zphiNest(res, 68)
    @test z ≈ -0.1666666666667
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiNest(res, 69)
    @test z ≈ -0.1666666666667
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 70)
    @test z ≈ 0
    @test phi ≈ 0.5890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 71)
    @test z ≈ 0
    @test phi ≈ 0.1963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 72)
    @test z ≈ 0.1666666666667
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 73)
    @test z ≈ -0.1666666666667
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 74)
    @test z ≈ 0
    @test phi ≈ 6.08683576633
    (z, phi) = Healpix.pix2zphiNest(res, 75)
    @test z ≈ 0
    @test phi ≈ 5.694136684632
    (z, phi) = Healpix.pix2zphiNest(res, 76)
    @test z ≈ 0.1666666666667
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 77)
    @test z ≈ 0.1666666666667
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiNest(res, 78)
    @test z ≈ 0.3333333333333
    @test phi ≈ 0.1963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 79)
    @test z ≈ 0.3333333333333
    @test phi ≈ 6.08683576633
    (z, phi) = Healpix.pix2zphiNest(res, 80)
    @test z ≈ 0.5
    @test phi ≈ 0
    (z, phi) = Healpix.pix2zphiNest(res, 81)
    @test z ≈ -0.5
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiNest(res, 82)
    @test z ≈ -0.3333333333333
    @test phi ≈ 1.767145867644
    (z, phi) = Healpix.pix2zphiNest(res, 83)
    @test z ≈ -0.3333333333333
    @test phi ≈ 1.374446785946
    (z, phi) = Healpix.pix2zphiNest(res, 84)
    @test z ≈ -0.1666666666667
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiNest(res, 85)
    @test z ≈ -0.1666666666667
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 86)
    @test z ≈ 0
    @test phi ≈ 2.159844949343
    (z, phi) = Healpix.pix2zphiNest(res, 87)
    @test z ≈ 0
    @test phi ≈ 1.767145867644
    (z, phi) = Healpix.pix2zphiNest(res, 88)
    @test z ≈ 0.1666666666667
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 89)
    @test z ≈ -0.1666666666667
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiNest(res, 90)
    @test z ≈ 0
    @test phi ≈ 1.374446785946
    (z, phi) = Healpix.pix2zphiNest(res, 91)
    @test z ≈ 0
    @test phi ≈ 0.9817477042468
    (z, phi) = Healpix.pix2zphiNest(res, 92)
    @test z ≈ 0.1666666666667
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiNest(res, 93)
    @test z ≈ 0.1666666666667
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiNest(res, 94)
    @test z ≈ 0.3333333333333
    @test phi ≈ 1.767145867644
    (z, phi) = Healpix.pix2zphiNest(res, 95)
    @test z ≈ 0.3333333333333
    @test phi ≈ 1.374446785946
    (z, phi) = Healpix.pix2zphiNest(res, 96)
    @test z ≈ 0.5
    @test phi ≈ 1.570796326795
    (z, phi) = Healpix.pix2zphiNest(res, 97)
    @test z ≈ -0.5
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiNest(res, 98)
    @test z ≈ -0.3333333333333
    @test phi ≈ 3.337942194439
    (z, phi) = Healpix.pix2zphiNest(res, 99)
    @test z ≈ -0.3333333333333
    @test phi ≈ 2.94524311274
    (z, phi) = Healpix.pix2zphiNest(res, 100)
    @test z ≈ -0.1666666666667
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiNest(res, 101)
    @test z ≈ -0.1666666666667
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiNest(res, 102)
    @test z ≈ 0
    @test phi ≈ 3.730641276138
    (z, phi) = Healpix.pix2zphiNest(res, 103)
    @test z ≈ 0
    @test phi ≈ 3.337942194439
    (z, phi) = Healpix.pix2zphiNest(res, 104)
    @test z ≈ 0.1666666666667
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiNest(res, 105)
    @test z ≈ -0.1666666666667
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiNest(res, 106)
    @test z ≈ 0
    @test phi ≈ 2.94524311274
    (z, phi) = Healpix.pix2zphiNest(res, 107)
    @test z ≈ 0
    @test phi ≈ 2.552544031042
    (z, phi) = Healpix.pix2zphiNest(res, 108)
    @test z ≈ 0.1666666666667
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiNest(res, 109)
    @test z ≈ 0.1666666666667
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiNest(res, 110)
    @test z ≈ 0.3333333333333
    @test phi ≈ 3.337942194439
    (z, phi) = Healpix.pix2zphiNest(res, 111)
    @test z ≈ 0.3333333333333
    @test phi ≈ 2.94524311274
    (z, phi) = Healpix.pix2zphiNest(res, 112)
    @test z ≈ 0.5
    @test phi ≈ 3.14159265359
    (z, phi) = Healpix.pix2zphiNest(res, 113)
    @test z ≈ -0.5
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiNest(res, 114)
    @test z ≈ -0.3333333333333
    @test phi ≈ 4.908738521234
    (z, phi) = Healpix.pix2zphiNest(res, 115)
    @test z ≈ -0.3333333333333
    @test phi ≈ 4.516039439535
    (z, phi) = Healpix.pix2zphiNest(res, 116)
    @test z ≈ -0.1666666666667
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiNest(res, 117)
    @test z ≈ -0.1666666666667
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiNest(res, 118)
    @test z ≈ 0
    @test phi ≈ 5.301437602933
    (z, phi) = Healpix.pix2zphiNest(res, 119)
    @test z ≈ 0
    @test phi ≈ 4.908738521234
    (z, phi) = Healpix.pix2zphiNest(res, 120)
    @test z ≈ 0.1666666666667
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiNest(res, 121)
    @test z ≈ -0.1666666666667
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiNest(res, 122)
    @test z ≈ 0
    @test phi ≈ 4.516039439535
    (z, phi) = Healpix.pix2zphiNest(res, 123)
    @test z ≈ 0
    @test phi ≈ 4.123340357837
    (z, phi) = Healpix.pix2zphiNest(res, 124)
    @test z ≈ 0.1666666666667
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiNest(res, 125)
    @test z ≈ 0.1666666666667
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiNest(res, 126)
    @test z ≈ 0.3333333333333
    @test phi ≈ 4.908738521234
    (z, phi) = Healpix.pix2zphiNest(res, 127)
    @test z ≈ 0.3333333333333
    @test phi ≈ 4.516039439535
    (z, phi) = Healpix.pix2zphiNest(res, 128)
    @test z ≈ 0.5
    @test phi ≈ 4.712388980385
    (z, phi) = Healpix.pix2zphiNest(res, 129)
    @test z ≈ -0.9791666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 130)
    @test z ≈ -0.9166666666667
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiNest(res, 131)
    @test z ≈ -0.9166666666667
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 132)
    @test z ≈ -0.8125
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 133)
    @test z ≈ -0.8125
    @test phi ≈ 1.308996938996
    (z, phi) = Healpix.pix2zphiNest(res, 134)
    @test z ≈ -0.6666666666667
    @test phi ≈ 1.374446785946
    (z, phi) = Healpix.pix2zphiNest(res, 135)
    @test z ≈ -0.6666666666667
    @test phi ≈ 0.9817477042468
    (z, phi) = Healpix.pix2zphiNest(res, 136)
    @test z ≈ -0.5
    @test phi ≈ 1.178097245096
    (z, phi) = Healpix.pix2zphiNest(res, 137)
    @test z ≈ -0.8125
    @test phi ≈ 0.2617993877991
    (z, phi) = Healpix.pix2zphiNest(res, 138)
    @test z ≈ -0.6666666666667
    @test phi ≈ 0.5890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 139)
    @test z ≈ -0.6666666666667
    @test phi ≈ 0.1963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 140)
    @test z ≈ -0.5
    @test phi ≈ 0.3926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 141)
    @test z ≈ -0.5
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 142)
    @test z ≈ -0.3333333333333
    @test phi ≈ 0.9817477042468
    (z, phi) = Healpix.pix2zphiNest(res, 143)
    @test z ≈ -0.3333333333333
    @test phi ≈ 0.5890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 144)
    @test z ≈ -0.1666666666667
    @test phi ≈ 0.7853981633974
    (z, phi) = Healpix.pix2zphiNest(res, 145)
    @test z ≈ -0.9791666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 146)
    @test z ≈ -0.9166666666667
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiNest(res, 147)
    @test z ≈ -0.9166666666667
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 148)
    @test z ≈ -0.8125
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 149)
    @test z ≈ -0.8125
    @test phi ≈ 2.879793265791
    (z, phi) = Healpix.pix2zphiNest(res, 150)
    @test z ≈ -0.6666666666667
    @test phi ≈ 2.94524311274
    (z, phi) = Healpix.pix2zphiNest(res, 151)
    @test z ≈ -0.6666666666667
    @test phi ≈ 2.552544031042
    (z, phi) = Healpix.pix2zphiNest(res, 152)
    @test z ≈ -0.5
    @test phi ≈ 2.748893571891
    (z, phi) = Healpix.pix2zphiNest(res, 153)
    @test z ≈ -0.8125
    @test phi ≈ 1.832595714594
    (z, phi) = Healpix.pix2zphiNest(res, 154)
    @test z ≈ -0.6666666666667
    @test phi ≈ 2.159844949343
    (z, phi) = Healpix.pix2zphiNest(res, 155)
    @test z ≈ -0.6666666666667
    @test phi ≈ 1.767145867644
    (z, phi) = Healpix.pix2zphiNest(res, 156)
    @test z ≈ -0.5
    @test phi ≈ 1.963495408494
    (z, phi) = Healpix.pix2zphiNest(res, 157)
    @test z ≈ -0.5
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 158)
    @test z ≈ -0.3333333333333
    @test phi ≈ 2.552544031042
    (z, phi) = Healpix.pix2zphiNest(res, 159)
    @test z ≈ -0.3333333333333
    @test phi ≈ 2.159844949343
    (z, phi) = Healpix.pix2zphiNest(res, 160)
    @test z ≈ -0.1666666666667
    @test phi ≈ 2.356194490192
    (z, phi) = Healpix.pix2zphiNest(res, 161)
    @test z ≈ -0.9791666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 162)
    @test z ≈ -0.9166666666667
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiNest(res, 163)
    @test z ≈ -0.9166666666667
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiNest(res, 164)
    @test z ≈ -0.8125
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 165)
    @test z ≈ -0.8125
    @test phi ≈ 4.450589592586
    (z, phi) = Healpix.pix2zphiNest(res, 166)
    @test z ≈ -0.6666666666667
    @test phi ≈ 4.516039439535
    (z, phi) = Healpix.pix2zphiNest(res, 167)
    @test z ≈ -0.6666666666667
    @test phi ≈ 4.123340357837
    (z, phi) = Healpix.pix2zphiNest(res, 168)
    @test z ≈ -0.5
    @test phi ≈ 4.319689898686
    (z, phi) = Healpix.pix2zphiNest(res, 169)
    @test z ≈ -0.8125
    @test phi ≈ 3.403392041389
    (z, phi) = Healpix.pix2zphiNest(res, 170)
    @test z ≈ -0.6666666666667
    @test phi ≈ 3.730641276138
    (z, phi) = Healpix.pix2zphiNest(res, 171)
    @test z ≈ -0.6666666666667
    @test phi ≈ 3.337942194439
    (z, phi) = Healpix.pix2zphiNest(res, 172)
    @test z ≈ -0.5
    @test phi ≈ 3.534291735289
    (z, phi) = Healpix.pix2zphiNest(res, 173)
    @test z ≈ -0.5
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 174)
    @test z ≈ -0.3333333333333
    @test phi ≈ 4.123340357837
    (z, phi) = Healpix.pix2zphiNest(res, 175)
    @test z ≈ -0.3333333333333
    @test phi ≈ 3.730641276138
    (z, phi) = Healpix.pix2zphiNest(res, 176)
    @test z ≈ -0.1666666666667
    @test phi ≈ 3.926990816987
    (z, phi) = Healpix.pix2zphiNest(res, 177)
    @test z ≈ -0.9791666666667
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiNest(res, 178)
    @test z ≈ -0.9166666666667
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 179)
    @test z ≈ -0.9166666666667
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiNest(res, 180)
    @test z ≈ -0.8125
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiNest(res, 181)
    @test z ≈ -0.8125
    @test phi ≈ 6.02138591938
    (z, phi) = Healpix.pix2zphiNest(res, 182)
    @test z ≈ -0.6666666666667
    @test phi ≈ 6.08683576633
    (z, phi) = Healpix.pix2zphiNest(res, 183)
    @test z ≈ -0.6666666666667
    @test phi ≈ 5.694136684632
    (z, phi) = Healpix.pix2zphiNest(res, 184)
    @test z ≈ -0.5
    @test phi ≈ 5.890486225481
    (z, phi) = Healpix.pix2zphiNest(res, 185)
    @test z ≈ -0.8125
    @test phi ≈ 4.974188368184
    (z, phi) = Healpix.pix2zphiNest(res, 186)
    @test z ≈ -0.6666666666667
    @test phi ≈ 5.301437602933
    (z, phi) = Healpix.pix2zphiNest(res, 187)
    @test z ≈ -0.6666666666667
    @test phi ≈ 4.908738521234
    (z, phi) = Healpix.pix2zphiNest(res, 188)
    @test z ≈ -0.5
    @test phi ≈ 5.105088062083
    (z, phi) = Healpix.pix2zphiNest(res, 189)
    @test z ≈ -0.5
    @test phi ≈ 5.497787143782
    (z, phi) = Healpix.pix2zphiNest(res, 190)
    @test z ≈ -0.3333333333333
    @test phi ≈ 5.694136684632
    (z, phi) = Healpix.pix2zphiNest(res, 191)
    @test z ≈ -0.3333333333333
    @test phi ≈ 5.301437602933
    (z, phi) = Healpix.pix2zphiNest(res, 192)
    @test z ≈ -0.1666666666667
    @test phi ≈ 5.497787143782
end


## Test udgrade
A = Healpix.HealpixMap{Float64, Healpix.NestedOrder}(1.0:Healpix.nside2npix(4))
nest_ref = [8.5,24.5,40.5,56.5,72.5,88.5,104.5,120.5,136.5,152.5,168.5,184.5]
@test Healpix.udgrade(A, 1).pixels ≈ nest_ref
@test Healpix.udgrade(A, 8).pixels ≈ repeat(A.pixels, inner=4)

ring_ref = [30.0625,33.4375,36.8125,40.1875,94.75,92.75,96.75,100.75,153.0625,156.4375,
    159.8125,163.1875]
A = Healpix.HealpixMap{Float64, Healpix.RingOrder}(1.0:Healpix.nside2npix(4))
@test Healpix.udgrade(A, 1).pixels ≈ ring_ref


unseen_ref = [9.0,24.5,40.5,56.5,72.5,88.5,104.5,120.5,136.5,152.5,168.5,184.5]
A = Healpix.HealpixMap{Float64, Healpix.NestedOrder}(1.0:Healpix.nside2npix(4))
A[1] = Healpix.UNSEEN
Healpix.udgrade(A, 1).pixels
@test Healpix.udgrade(A, 1).pixels ≈ unseen_ref
@test Healpix.udgrade(A, 1; pess=true) ≈ [Healpix.UNSEEN, unseen_ref[2:end]...]

@test A ≈ Healpix.udgrade(A, 4)
