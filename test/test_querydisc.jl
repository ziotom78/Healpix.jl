# NSIDE = 1
resol = Healpix.Resolution(1)
(theta, phi) = (0.689126661018239, 1.0582462744482535)
radius = 0.058619561676582874
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 2 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1 ]

(theta, phi) = (1.4140877251282506, 0.6744301121442584)
radius = 0.01695673391711479
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 5 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1, 5 ]

(theta, phi) = (0.8608403810587308, 1.1599355989295899)
radius = 0.04460724018855071
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 2, 6 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1 ]

(theta, phi) = (0.41161737458473074, 1.3067332065505721)
radius = 0.14661027372966898
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 2 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1, 2 ]

(theta, phi) = (1.0083251189919613, 1.5282558373734476)
radius = 0.14507605046902744
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 2, 6 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1, 2, 6 ]

(theta, phi) = (1.24579812272374, 1.2393799362503723)
radius = 0.15056798001463167
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 6 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1, 6 ]

(theta, phi) = (1.0595618910606672, 0.7796227404700381)
radius = 0.12662589083197673
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 5, 6 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1 ]

(theta, phi) = (0.9947752927045007, 1.3916937464463668)
radius = 0.06181036072925033
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 2, 6 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1, 6 ]

(theta, phi) = (1.1052388022518542, 0.5024882388747463)
radius = 0.023733482672056272
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 5 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1 ]


# NSIDE = 2
resol = Healpix.Resolution(2)
(theta, phi) = (0.8274206694161718, 0.6778863939924552)
radius = 0.06415569523479893
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 5, 6, 14 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1, 5, 6, 14 ]

(theta, phi) = (0.5999273204030076, 1.1322040673838987)
radius = 0.05424135575858182
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 5, 6 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1, 6 ]

(theta, phi) = (1.2756980126463742, 0.2589674949713131)
radius = 0.08330004470187188
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 5, 13, 14, 21 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 5, 13, 21 ]

(theta, phi) = (1.0408783167453979, 0.962528002658057)
radius = 0.08859347139983494
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 5, 6, 14, 22 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 6, 14 ]

(theta, phi) = (0.7963089271915055, 0.003904968819053311)
radius = 0.05071562818698281
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 5, 12, 13 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 5, 12, 13 ]

(theta, phi) = (1.0322098664193138, 0.13309315730222127)
radius = 0.00011432024368608222
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 5, 12, 13 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 5, 13 ]

(theta, phi) = (0.9685697720544718, 1.405628356841054)
radius = 0.0017668454106800381
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 6, 7, 15 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 6, 15 ]

(theta, phi) = (0.6017922080874089, 0.6836348548676232)
radius = 0.07305780571992057
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 5, 6 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1, 5 ]

(theta, phi) = (0.9022067535596878, 0.1035934739894827)
radius = 0.047596326377552685
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 5, 12, 13 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 5, 12, 13 ]


# NSIDE = 4
resol = Healpix.Resolution(4)
(theta, phi) = (1.515165746791338, 0.635024185382709)
radius = 0.01162590354453882
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 74, 75, 90 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 90 ]

(theta, phi) = (1.362037560020116, 0.2161573215303783)
radius = 0.026052162498527282
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 57, 73, 74, 89 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 57, 73, 74, 89 ]

(theta, phi) = (0.4763613735907769, 1.117758874772761)
radius = 0.030679666541001314
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 6, 14, 15 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 6 ]

(theta, phi) = (0.9459783755784011, 0.7942160339012698)
radius = 0.0353334764739088
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 26, 27, 43 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 26, 27, 43 ]

(theta, phi) = (1.2912679807917973, 0.8440209515955615)
radius = 0.017244528032354524
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 43, 58, 59, 75 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 59, 75 ]

(theta, phi) = (1.079464814945391, 0.035490466577281006)
radius = 0.002631117290719266
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 41, 57 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 41 ]

(theta, phi) = (0.6739597536634145, 1.0511209675498052)
radius = 0.007346847528606954
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 14, 15, 27 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 14, 15, 27 ]

(theta, phi) = (1.1314024983092204, 0.8335890264333707)
radius = 0.03583874270454281
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 43, 58, 59, 75 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 43, 59 ]

(theta, phi) = (1.2726758389245316, 0.7783261238489829)
radius = 0.028540674536373753
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 43, 58, 59, 75 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 43, 58, 59, 75 ]

(theta, phi) = (1.3574870394971803, 1.0111065342319239)
radius = 0.006429328420746328
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 59, 75, 76, 91 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 59, 76 ]


# NSIDE = 8
resol = Healpix.Resolution(8)
(theta, phi) = (0.8076806206890037, 0.7952557360370545)
radius = 0.01102721479796156
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 88, 116, 117, 149 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 88, 116, 117 ]

(theta, phi) = (1.3817261212420426, 0.47307073901260505)
radius = 0.018575838226070014
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 275, 276, 307, 339 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 307 ]

(theta, phi) = (1.5600877921358713, 1.3876539929454466)
radius = 0.005063882107503555
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 344, 375, 376, 408 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 344, 376, 408 ]

(theta, phi) = (1.091284114698803, 0.17501925444952274)
radius = 0.009184678411068776
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 177, 178, 210 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 177, 210 ]

(theta, phi) = (0.09701774177400761, 1.1057641283541604)
radius = 0.020168845590830178
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1, 2, 5, 6 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1, 6 ]

(theta, phi) = (0.2691412953817959, 1.1392606263960878)
radius = 9.903440123815941e-05
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 6, 14, 15 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 6, 15 ]

(theta, phi) = (1.0216479231322462, 0.6727245529779553)
radius = 0.01928243290401354
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 148, 149, 180, 212 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 180 ]

(theta, phi) = (1.030238817205236, 1.5057632873503757)
radius = 0.020209816327264818
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 152, 153, 184, 185, 217 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 153, 184 ]

(theta, phi) = (1.0927608639812159, 1.5433716553328818)
radius = 0.020528161007589297
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 153, 184, 185, 217 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 184, 217 ]

(theta, phi) = (0.5638465893912242, 0.7535965734705176)
radius = 0.009240198128134463
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 43, 63, 64 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 43, 63 ]


# NSIDE = 16
resol = Healpix.Resolution(16)
(theta, phi) = (1.5341049101993827, 1.3109663847459547)
radius = 0.005717969122250481
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1390, 1454, 1455, 1518 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1390, 1454, 1518 ]

(theta, phi) = (1.5018211299114723, 0.5014786542516827)
radius = 0.0016190814132463777
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1318, 1381, 1382, 1446 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1382, 1446 ]

(theta, phi) = (1.3038894126099487, 1.043600462743714)
radius = 0.0010117144033290762
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1067, 1068, 1131 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1131 ]

(theta, phi) = (0.7123347304324052, 0.08340174946183303)
radius = 0.00996569896486288
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 313, 314, 365, 366 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 313, 365 ]

(theta, phi) = (0.9606611903817225, 0.18185603658521668)
radius = 0.007701318985326741
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 547, 610, 611, 675 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 610, 675 ]

(theta, phi) = (0.6134272449395771, 1.4953667092786198)
radius = 0.005501346665823176
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 231, 275, 276, 325 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 231, 276 ]

(theta, phi) = (1.5586002974302926, 1.067911363575043)
radius = 0.004599875374474255
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1452, 1515, 1516, 1580 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1452, 1515 ]

(theta, phi) = (0.8169257109714797, 1.0695178008515478)
radius = 0.00021976994020986404
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 430, 431, 491, 492 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 431, 491 ]

(theta, phi) = (1.1226016584627345, 1.3356194220952733)
radius = 0.008107569131018387
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 814, 815, 878, 943 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 815, 878 ]

(theta, phi) = (1.262368100056086, 0.41309822061848983)
radius = 0.010938668106977556
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 996, 997, 1061, 1062, 1125 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 997, 1061 ]


# NSIDE = 32
resol = Healpix.Resolution(32)
(theta, phi) = (1.1041178494634607, 0.7223512801224912)
radius = 0.004230338859433805
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 3152, 3279, 3280, 3407, 3408 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 3279, 3408 ]

(theta, phi) = (0.8640509276017856, 1.3351377045019193)
radius = 0.0015327504888317854
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 2011, 2012, 2140, 2141, 2268 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 2140 ]

(theta, phi) = (0.2305395271389823, 1.169712232120978)
radius = 0.005413378345406392
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 118, 119, 151, 152, 188 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 151, 188 ]

(theta, phi) = (1.5136684013393589, 0.6469295852088414)
radius = 0.004043721926152334
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 5582, 5710, 5837, 5838 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 5710, 5838 ]

(theta, phi) = (0.5683354153554419, 1.3176511461918188)
radius = 0.0010100848617587486
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 858, 943, 1032 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 943 ]

(theta, phi) = (1.4636377121902646, 0.9738226880978805)
radius = 0.004289830651352997
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 5332, 5333, 5461, 5588 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 5461 ]

(theta, phi) = (0.6089685551576646, 0.6870486621978211)
radius = 0.0020971469387025487
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1022, 1023, 1115 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1023, 1115 ]

(theta, phi) = (0.5995210175626202, 0.5276658514410062)
radius = 0.002221077032839156
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 932, 1020, 1021, 1112, 1113 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1020 ]

(theta, phi) = (1.0465022470086962, 0.5805581105305407)
radius = 0.0017530588498535047
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 2893, 3020, 3021, 3149 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 2893, 3020, 3149 ]

(theta, phi) = (0.9929190379217407, 1.1050799003318914)
radius = 0.002694676917431618
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 2647, 2648, 2775 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 2775 ]


# NSIDE = 64
resol = Healpix.Resolution(64)
(theta, phi) = (1.0821121265621954, 1.2559636677618118)
radius = 0.0006050770961973326
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 12724, 12980, 13236 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 12980 ]

(theta, phi) = (0.7832698855515476, 1.2142436658472298)
radius = 0.0011921473384799543
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 6890, 7126, 7127, 7368 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 6890, 7127 ]

(theta, phi) = (0.6120767192768516, 1.2365233621923988)
radius = 0.0009593364958215149
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 4177, 4361, 4362, 4550, 4551 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 4361, 4362, 4550 ]

(theta, phi) = (0.45782657645221075, 0.9258177253178457)
radius = 0.002078863644146202
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 2401, 2402, 2541, 2542 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 2401, 2542 ]

(theta, phi) = (1.248616053388374, 1.2638598718334468)
radius = 0.0008800861192997422
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 16564, 16565, 16820 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 16820 ]

(theta, phi) = (0.9988924082178857, 1.1132824744512053)
radius = 8.478251894957192e-05
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 10926, 11182, 11438 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 11182 ]

(theta, phi) = (1.097667503988652, 1.0897805604713573)
radius = 0.0017126313743913255
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 12973, 13229, 13485, 13486 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 13229 ]

(theta, phi) = (0.3799603230734639, 0.530354349459068)
radius = 0.0009005529055655211
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 1634, 1635, 1750, 1751 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 1634, 1751 ]

(theta, phi) = (1.2110594057802377, 0.2502984766180111)
radius = 0.002694907287006249
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 15499, 15754, 15755, 16011 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 15499, 15755, 16011 ]

(theta, phi) = (0.6194331412908418, 0.11354740398959659)
radius = 0.0027518068203280068
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=0) == [  ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=1) == [ 4328, 4515, 4516, 4708 ]
@test Healpix.queryDiscRing(resol, theta, phi, radius, fact=4) == [ 4328, 4516 ]
