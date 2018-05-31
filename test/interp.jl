resol = Healpix.Resolution(4)

function compare(r1::Healpix.RingInfo, r2::Healpix.RingInfo)
    @test r1.ring == r2.ring
    @test r1.firstpixidx == r2.firstpixidx
    @test r1.numofpixels == r2.numofpixels
    @test r1.colatitude_rad ≈ r2.colatitude_rad
    @test r1.shifted == r2.shifted
end

compare(Healpix.getringinfo(resol,  1), Healpix.RingInfo(1,    1,  4, 2.044801990e-01, true))
compare(Healpix.getringinfo(resol,  2), Healpix.RingInfo(2,    5,  8, 4.111378623e-01, true))
compare(Healpix.getringinfo(resol,  3), Healpix.RingInfo(3,   13, 12, 6.223684886e-01, true))
compare(Healpix.getringinfo(resol,  4), Healpix.RingInfo(4,   25, 16, 8.410686706e-01, true))
compare(Healpix.getringinfo(resol,  5), Healpix.RingInfo(5,   41, 16, 1.047197551e+00, false))
compare(Healpix.getringinfo(resol,  6), Healpix.RingInfo(6,   57, 16, 1.230959417e+00, true))
compare(Healpix.getringinfo(resol,  7), Healpix.RingInfo(7,   73, 16, 1.403348248e+00, false))
compare(Healpix.getringinfo(resol,  8), Healpix.RingInfo(8,   89, 16, 1.570796327e+00, true))
compare(Healpix.getringinfo(resol,  9), Healpix.RingInfo(9,  105, 16, 1.738244406e+00, false))
compare(Healpix.getringinfo(resol, 10), Healpix.RingInfo(10, 121, 16, 1.910633236e+00, true))
compare(Healpix.getringinfo(resol, 11), Healpix.RingInfo(11, 137, 16, 2.094395102e+00, false))
compare(Healpix.getringinfo(resol, 12), Healpix.RingInfo(12, 153, 16, 2.300523983e+00, true))
compare(Healpix.getringinfo(resol, 13), Healpix.RingInfo(13, 169, 12, 2.519224165e+00, true))
compare(Healpix.getringinfo(resol, 14), Healpix.RingInfo(14, 181,  8, 2.730454791e+00, true))
compare(Healpix.getringinfo(resol, 15), Healpix.RingInfo(15, 189,  4, 2.937112455e+00, true))

(idx, w) = Healpix.getinterpolRing(resol, 2.149214549e-01, 7.772156558e-01)
@test idx == [3, 0, 4, 5] + 1
@test w ≈ [4.945957119e-03, 9.445296357e-01, 2.578858162e-02, 2.473582554e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.163955573e-01, 2.378372948e+00)
@test idx == [1, 2, 6, 7] + 1
@test w ≈ [9.290373638e-01, 1.330516466e-02, 2.720057586e-02, 3.045689563e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.615696557e-01, 3.913463913e+00)
@test idx == [3, 0, 1, 2] + 1
@test w ≈ [5.246295665e-02, 5.246295665e-02, 5.926731317e-02, 8.358067735e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.820288832e-01, 5.453327200e+00)
@test idx == [0, 1, 2, 3] + 1
@test w ≈ [2.744925417e-02, 2.744925417e-02, 5.264563016e-02, 8.924558615e-01]

(idx, w) = Healpix.getinterpolRing(resol, 3.801627711e-01, 3.487819350e-01)
@test idx == [3, 0, 11, 4] + 1
@test w ≈ [4.166209221e-02, 1.082239110e-01, 4.753586505e-02, 8.025781317e-01]

(idx, w) = Healpix.getinterpolRing(resol, 3.841583285e-01, 1.205197444e+00)
@test idx == [0, 1, 5, 6] + 1
@test w ≈ [9.566151592e-02, 3.489030289e-02, 8.394478326e-01, 3.000034856e-02]

(idx, w) = Healpix.getinterpolRing(resol, 4.235101843e-01, 1.959967433e+00)
@test idx == [5, 6, 15, 16] + 1
@test w ≈ [4.228852420e-03, 9.371985710e-01, 4.432409073e-02, 1.424848581e-02]

(idx, w) = Healpix.getinterpolRing(resol, 4.049509500e-01, 2.732246699e+00)
@test idx == [1, 2, 6, 7] + 1
@test w ≈ [2.277075717e-02, 7.167219647e-03, 2.056090835e-02, 9.495011148e-01]

(idx, w) = Healpix.getinterpolRing(resol, 4.026864270e-01, 3.523315887e+00)
@test idx == [1, 2, 7, 8] + 1
@test w ≈ [1.050971310e-02, 3.038611055e-02, 1.340336926e-02, 9.457008071e-01]

(idx, w) = Healpix.getinterpolRing(resol, 3.992881129e-01, 4.295575726e+00)
@test idx == [2, 3, 8, 9] + 1
@test w ≈ [4.388525470e-02, 1.345474045e-02, 2.894260210e-02, 9.137174028e-01]

(idx, w) = Healpix.getinterpolRing(resol, 3.911012261e-01, 5.130210349e+00)
@test idx == [2, 3, 10, 11] + 1
@test w ≈ [2.268827654e-02, 7.426741266e-02, 8.741589128e-01, 2.888539796e-02]

(idx, w) = Healpix.getinterpolRing(resol, 3.640174506e-01, 5.883317214e+00)
@test idx == [3, 0, 10, 11] + 1
@test w ≈ [1.720495755e-01, 5.596234990e-02, 7.046605832e-03, 7.649414687e-01]

(idx, w) = Healpix.getinterpolRing(resol, 5.866457293e-01, 2.759300775e-01)
@test idx == [11, 4, 12, 13] + 1
@test w ≈ [2.514350428e-02, 1.439738245e-01, 8.084591176e-01, 2.242355356e-02]

(idx, w) = Healpix.getinterpolRing(resol, 6.448612764e-01, 7.523113308e-01)
@test idx == [12, 13, 25, 26] + 1
@test w ≈ [5.669213305e-02, 8.404602716e-01, 6.008921433e-02, 4.275838098e-02]

(idx, w) = Healpix.getinterpolRing(resol, 6.175323839e-01, 1.335645119e+00)
@test idx == [5, 6, 14, 15] + 1
@test w ≈ [1.830227282e-02, 4.592630150e-03, 9.273760418e-01, 4.972905527e-02]

(idx, w) = Healpix.getinterpolRing(resol, 6.676720908e-01, 1.797850703e+00)
@test idx == [14, 15, 28, 29] + 1
@test w ≈ [5.261204923e-02, 7.402385930e-01, 1.909525107e-01, 1.619684709e-02]

(idx, w) = Healpix.getinterpolRing(resol, 6.297404338e-01, 2.394172483e+00)
@test idx == [16, 17, 29, 30] + 1
@test w ≈ [8.962043091e-01, 7.008769308e-02, 1.359409293e-02, 2.011390489e-02]

(idx, w) = Healpix.getinterpolRing(resol, 6.160658501e-01, 2.842695991e+00)
@test idx == [7, 8, 16, 17] + 1
@test w ≈ [2.627410560e-02, 3.563605882e-03, 6.873655707e-02, 9.014257314e-01]

(idx, w) = Healpix.getinterpolRing(resol, 6.156440114e-01, 3.425196760e+00)
@test idx == [7, 8, 18, 19] + 1
@test w ≈ [4.421977241e-03, 2.741278682e-02, 9.278470164e-01, 4.031821955e-02]

(idx, w) = Healpix.getinterpolRing(resol, 5.865905571e-01, 3.954832255e+00)
@test idx == [8, 9, 19, 20] + 1
@test w ≈ [7.868499277e-02, 9.069353005e-02, 7.864546513e-01, 4.416682592e-02]

(idx, w) = Healpix.getinterpolRing(resol, 6.146918563e-01, 4.474086218e+00)
@test idx == [9, 10, 20, 21] + 1
@test w ≈ [2.919810024e-02, 7.144320205e-03, 9.204132049e-01, 4.324437464e-02]

(idx, w) = Healpix.getinterpolRing(resol, 6.241163390e-01, 4.949489045e+00)
@test idx == [20, 21, 36, 37] + 1
@test w ≈ [4.679523224e-02, 9.452127743e-01, 7.162661404e-03, 8.293320078e-04]

(idx, w) = Healpix.getinterpolRing(resol, 6.193384905e-01, 5.500182190e+00)
@test idx == [10, 11, 22, 23] + 1
@test w ≈ [7.128507683e-03, 7.215993880e-03, 9.811469111e-01, 4.508587337e-03]

(idx, w) = Healpix.getinterpolRing(resol, 6.036427361e-01, 6.020782191e+00)
@test idx == [11, 4, 22, 23] + 1
@test w ≈ [7.394376608e-02, 1.470697917e-02, 1.050819043e-03, 9.102984357e-01]

(idx, w) = Healpix.getinterpolRing(resol, 8.274281073e-01, 2.151641185e-01)
@test idx == [23, 12, 24, 25] + 1
@test w ≈ [5.555191388e-03, 5.681587353e-02, 8.927062622e-01, 4.492267293e-02]

(idx, w) = Healpix.getinterpolRing(resol, 8.479960921e-01, 5.677950591e-01)
@test idx == [24, 25, 41, 42] + 1
@test w ≈ [5.230287239e-02, 9.140898937e-01, 1.862249945e-02, 1.498473450e-02]

(idx, w) = Healpix.getinterpolRing(resol, 8.530366695e-01, 9.813302854e-01)
@test idx == [25, 26, 42, 43] + 1
@test w ≈ [1.001232813e-03, 9.409380099e-01, 2.909209423e-02, 2.896866306e-02]

(idx, w) = Healpix.getinterpolRing(resol, 8.118243739e-01, 1.414234280e+00)
@test idx == [14, 15, 27, 28] + 1
@test w ≈ [1.068427498e-01, 2.687591642e-02, 7.785114220e-01, 8.776991186e-02]

(idx, w) = Healpix.getinterpolRing(resol, 8.510939650e-01, 1.776456182e+00)
@test idx == [28, 29, 44, 45] + 1
@test w ≈ [9.288085188e-01, 2.255543158e-02, 2.316493602e-02, 2.547111356e-02]

(idx, w) = Healpix.getinterpolRing(resol, 8.469525322e-01, 2.111243588e+00)
@test idx == [28, 29, 45, 46] + 1
@test w ≈ [1.202296055e-01, 8.512258180e-01, 1.780503206e-02, 1.073954438e-02]

(idx, w) = Healpix.getinterpolRing(resol, 8.298442146e-01, 2.559754416e+00)
@test idx == [16, 17, 30, 31] + 1
@test w ≈ [3.137041286e-02, 1.995307416e-02, 9.312577744e-01, 1.741873863e-02]

(idx, w) = Healpix.getinterpolRing(resol, 8.351939785e-01, 2.991975062e+00)
@test idx == [17, 18, 31, 32] + 1
@test w ≈ [2.110666172e-02, 5.755190822e-03, 8.573328282e-01, 1.158053193e-01]

(idx, w) = Healpix.getinterpolRing(resol, 8.369234356e-01, 3.312134186e+00)
@test idx == [17, 18, 31, 32] + 1
@test w ≈ [3.303479765e-03, 1.565048254e-02, 6.447390906e-02, 9.165721286e-01]

(idx, w) = Healpix.getinterpolRing(resol, 7.935487577e-01, 3.757499820e+00)
@test idx == [18, 19, 33, 34] + 1
@test w ≈ [7.033548766e-02, 1.469478838e-01, 7.291829434e-01, 5.353368516e-02]

(idx, w) = Healpix.getinterpolRing(resol, 8.059051117e-01, 4.093045825e+00)
@test idx == [19, 20, 33, 34] + 1
@test w ≈ [1.097929030e-01, 5.099141143e-02, 6.474078631e-02, 7.744748992e-01]

(idx, w) = Healpix.getinterpolRing(resol, 8.826328785e-01, 4.496918030e+00)
@test idx == [34, 35, 51, 52] + 1
@test w ≈ [3.887387033e-02, 7.594842839e-01, 1.106393219e-01, 9.100252379e-02]

(idx, w) = Healpix.getinterpolRing(resol, 8.850317142e-01, 4.921077416e+00)
@test idx == [36, 37, 52, 53] + 1
@test w ≈ [7.620012603e-01, 2.471934151e-02, 9.993830310e-02, 1.133410951e-01]

(idx, w) = Healpix.getinterpolRing(resol, 8.049519283e-01, 5.288600692e+00)
@test idx == [21, 22, 36, 37] + 1
@test w ≈ [6.597727244e-02, 9.916544458e-02, 2.729058803e-02, 8.075666949e-01]

(idx, w) = Healpix.getinterpolRing(resol, 8.502428215e-01, 5.663061886e+00)
@test idx == [37, 38, 54, 55] + 1
@test w ≈ [7.560943765e-02, 8.798836942e-01, 2.577532155e-02, 1.873154664e-02]

(idx, w) = Healpix.getinterpolRing(resol, 8.447191377e-01, 6.098214209e+00)
@test idx == [39, 24, 55, 40] + 1
@test w ≈ [9.538285349e-01, 2.846183070e-02, 8.341681147e-03, 9.367953286e-03]

(idx, w) = Healpix.getinterpolRing(resol, 1.014475223e+00, -4.403746725e-02)
@test idx == [39, 24, 55, 40] + 1
@test w ≈ [9.717542604e-02, 6.157150822e-02, 9.433853059e-02, 7.469145352e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.067543701e+00, 3.863265116e-01)
@test idx == [40, 41, 56, 57] + 1
@test w ≈ [1.443089173e-02, 8.748489119e-01, 5.715682308e-02, 5.356337328e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.049067875e+00, 7.563528245e-01)
@test idx == [41, 42, 57, 58] + 1
@test w ≈ [7.321055114e-02, 9.166114741e-01, 5.841784462e-03, 4.336190282e-03]

(idx, w) = Healpix.getinterpolRing(resol, 1.032937195e+00, 1.149483697e+00)
@test idx == [26, 27, 42, 43] + 1
@test w ≈ [3.963171988e-02, 2.955002924e-02, 6.782295736e-02, 8.629952935e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.088084110e+00, 1.562225721e+00)
@test idx == [43, 44, 59, 60] + 1
@test w ≈ [1.696888931e-02, 7.605335927e-01, 1.161047380e-01, 1.063927799e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.054500634e+00, 1.922730895e+00)
@test idx == [44, 45, 60, 61] + 1
@test w ≈ [9.968051352e-02, 8.605773864e-01, 2.399651775e-02, 1.574558235e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.005940898e+00, 2.349579042e+00)
@test idx == [29, 30, 45, 46] + 1
@test w ≈ [1.034466403e-01, 9.670315310e-02, 1.347435722e-02, 7.863758494e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.017296147e+00, 2.703703526e+00)
@test idx == [30, 31, 46, 47] + 1
@test w ≈ [8.922389380e-02, 5.583779831e-02, 9.838245911e-02, 7.565558488e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.072244926e+00, 3.139970677e+00)
@test idx == [47, 48, 63, 64] + 1
@test w ≈ [3.567351803e-03, 8.601292079e-01, 6.871469833e-02, 6.758874195e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.081298280e+00, 3.517281899e+00)
@test idx == [48, 49, 64, 65] + 1
@test w ≈ [3.527718291e-02, 7.791526022e-01, 1.008231170e-01, 8.474709785e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.073788275e+00, 3.943102138e+00)
@test idx == [50, 51, 65, 66] + 1
@test w ≈ [8.202074943e-01, 3.509042904e-02, 6.641432573e-02, 7.828775091e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.034224901e+00, 4.366972270e+00)
@test idx == [34, 35, 51, 52] + 1
@test w ≈ [2.388977002e-02, 3.904488399e-02, 8.242393335e-01, 1.128260125e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.051768672e+00, 4.694258711e+00)
@test idx == [51, 52, 67, 68] + 1
@test w ≈ [4.501990264e-02, 9.301048543e-01, 1.358607052e-02, 1.128917258e-02]

(idx, w) = Healpix.getinterpolRing(resol, 9.997062844e-01, 5.094753809e+00)
@test idx == [36, 37, 52, 53] + 1
@test w ≈ [1.212610847e-01, 1.091349019e-01, 2.025286792e-02, 7.493511455e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.020925057e+00, 5.531339255e+00)
@test idx == [37, 38, 54, 55] + 1
@test w ≈ [5.283845693e-02, 7.461818368e-02, 7.979934738e-01, 7.454988555e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.046254099e+00, 5.882666420e+00)
@test idx == [38, 39, 54, 55] + 1
@test w ≈ [2.379642813e-03, 2.197359378e-03, 1.982183047e-02, 9.756011673e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.254361768e+00, 1.586402026e-01)
@test idx == [71, 56, 72, 73] + 1
@test w ≈ [8.299019033e-02, 7.812565204e-01, 8.091249558e-02, 5.484079366e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.190772965e+00, 5.526066525e-01)
@test idx == [41, 42, 56, 57] + 1
@test w ≈ [1.296377637e-01, 8.904989623e-02, 7.250478088e-02, 7.088075592e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.277322826e+00, 9.989285444e-01)
@test idx == [58, 59, 74, 75] + 1
@test w ≈ [6.990692359e-01, 3.198405694e-02, 1.227067602e-01, 1.462399469e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.242583797e+00, 1.376450097e+00)
@test idx == [59, 60, 75, 76] + 1
@test w ≈ [9.278114413e-01, 4.757396209e-03, 3.337158867e-02, 3.405957382e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.269247941e+00, 1.768885308e+00)
@test idx == [60, 61, 76, 77] + 1
@test w ≈ [7.744487652e-01, 3.445642649e-03, 1.100689910e-01, 1.120366012e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.183915065e+00, 2.143062129e+00)
@test idx == [45, 46, 60, 61] + 1
@test w ≈ [1.389445801e-01, 1.170625734e-01, 3.179609672e-02, 7.121967498e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.263822021e+00, 2.532194129e+00)
@test idx == [61, 62, 78, 79] + 1
@test w ≈ [4.194200239e-02, 7.674272970e-01, 1.051939475e-01, 8.543675309e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.189132118e+00, 2.916904929e+00)
@test idx == [47, 48, 62, 63] + 1
@test w ≈ [1.302338546e-01, 9.738300885e-02, 5.573716928e-02, 7.166459672e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.202692089e+00, 3.307287249e+00)
@test idx == [48, 49, 63, 64] + 1
@test w ≈ [8.892091593e-02, 6.490495144e-02, 6.605419578e-02, 7.801199368e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.277364937e+00, 3.744855187e+00)
@test idx == [65, 66, 81, 82] + 1
@test w ≈ [7.043570642e-01, 2.645194452e-02, 1.248520126e-01, 1.443389787e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.185252504e+00, 4.154705415e+00)
@test idx == [50, 51, 66, 67] + 1
@test w ≈ [1.044984180e-01, 1.442306247e-01, 6.912665995e-01, 6.000435776e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.274974832e+00, 4.525150828e+00)
@test idx == [67, 68, 83, 84] + 1
@test w ≈ [7.273957377e-01, 1.727788873e-02, 1.217391141e-01, 1.335872595e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.264897666e+00, 4.908605903e+00)
@test idx == [67, 68, 84, 85] + 1
@test w ≈ [2.712243336e-04, 8.028584316e-01, 9.850165699e-02, 9.836868711e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.228710819e+00, 5.312356398e+00)
@test idx == [53, 54, 69, 70] + 1
@test w ≈ [5.778010096e-03, 6.458468007e-03, 9.602992672e-01, 2.746425472e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.215572964e+00, 5.669032493e+00)
@test idx == [54, 55, 69, 70] + 1
@test w ≈ [4.721785603e-02, 3.651254000e-02, 5.857464039e-02, 8.576949636e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.234495479e+00, 6.041621309e+00)
@test idx == [70, 71, 87, 72] + 1
@test w ≈ [1.127759508e-01, 8.667119199e-01, 1.261778343e-02, 7.894345891e-03]

(idx, w) = Healpix.getinterpolRing(resol, 1.354725810e+00, 2.660903331e-02)
@test idx == [71, 56, 72, 73] + 1
@test w ≈ [1.219138980e-01, 1.601370779e-01, 6.693012671e-01, 4.864775700e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.410990220e+00, 3.872021369e-01)
@test idx == [72, 73, 88, 89] + 1
@test w ≈ [1.335902273e-02, 9.410031161e-01, 2.345776277e-02, 2.218009845e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.415816606e+00, 7.920680428e-01)
@test idx == [74, 75, 89, 90] + 1
@test w ≈ [9.098189486e-01, 1.572000932e-02, 3.596582195e-02, 3.849522012e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.360892021e+00, 1.184663480e+00)
@test idx == [58, 59, 75, 76] + 1
@test w ≈ [1.190228621e-01, 1.272589081e-01, 7.411154740e-01, 1.260275578e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.356423554e+00, 1.579361629e+00)
@test idx == [59, 60, 76, 77] + 1
@test w ≈ [1.301642071e-01, 1.420384276e-01, 7.119231135e-01, 1.587425189e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.398840917e+00, 1.952878468e+00)
@test idx == [60, 61, 76, 77] + 1
@test w ≈ [1.378003695e-02, 1.236626377e-02, 2.632892997e-02, 9.475247693e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.373892920e+00, 2.369046460e+00)
@test idx == [61, 62, 78, 79] + 1
@test w ≈ [7.984085171e-02, 9.102478410e-02, 8.019990580e-01, 2.713530614e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.427258878e+00, 2.778134829e+00)
@test idx == [79, 80, 94, 95] + 1
@test w ≈ [7.933762387e-01, 6.382946488e-02, 6.076436370e-02, 8.202993272e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.443786262e+00, 3.180727971e+00)
@test idx == [80, 81, 95, 96] + 1
@test w ≈ [6.829137070e-01, 7.559045224e-02, 9.668110468e-02, 1.448147361e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.354413791e+00, 3.574195571e+00)
@test idx == [64, 65, 81, 82] + 1
@test w ≈ [1.130861473e-01, 1.707748017e-01, 6.433690939e-01, 7.276995718e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.364397858e+00, 3.975989465e+00)
@test idx == [65, 66, 82, 83] + 1
@test w ≈ [8.478042199e-02, 1.411645535e-01, 6.774730576e-01, 9.658196692e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.371461931e+00, 4.322328249e+00)
@test idx == [66, 67, 83, 84] + 1
@test w ≈ [9.124100287e-02, 9.372641190e-02, 8.095567847e-01, 5.475800555e-03]

(idx, w) = Healpix.getinterpolRing(resol, 1.384622757e+00, 4.709659901e+00)
@test idx == [67, 68, 83, 84] + 1
@test w ≈ [5.506667184e-02, 5.355690328e-02, 6.194659401e-03, 8.851817655e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.390653496e+00, 5.075035017e+00)
@test idx == [68, 69, 84, 85] + 1
@test w ≈ [4.245575594e-02, 3.118446478e-02, 7.089380431e-02, 8.554659750e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.362350144e+00, 5.473559160e+00)
@test idx == [69, 70, 85, 86] + 1
@test w ≈ [1.335844850e-01, 1.042389502e-01, 4.702328606e-02, 7.151532788e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.360225362e+00, 5.878483663e+00)
@test idx == [70, 71, 86, 87] + 1
@test w ≈ [1.327200989e-01, 1.174288560e-01, 2.291865359e-02, 7.269323916e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.599241123e+00, 1.507301684e-01)
@test idx == [103, 88, 104, 105] + 1
@test w ≈ [9.643491612e-02, 7.336927465e-01, 1.046700305e-01, 6.520230685e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.579444393e+00, 5.400873787e-01)
@test idx == [88, 89, 105, 106] + 1
@test w ≈ [1.182395914e-01, 8.301141541e-01, 3.226231945e-02, 1.938393507e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.615350203e+00, 9.895368230e-01)
@test idx == [90, 91, 106, 107] + 1
@test w ≈ [7.193669693e-01, 1.455726091e-02, 1.277603178e-01, 1.383154520e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.585367886e+00, 1.384915925e+00)
@test idx == [91, 92, 107, 108] + 1
@test w ≈ [8.886391417e-01, 2.433950299e-02, 4.119073669e-02, 4.583061858e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.565702276e+00, 1.762951234e+00)
@test idx == [76, 77, 91, 92] + 1
@test w ≈ [1.553578742e-02, 1.488588638e-02, 1.035659576e-02, 9.592217304e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.584409713e+00, 2.198178287e+00)
@test idx == [93, 94, 109, 110] + 1
@test w ≈ [8.290218281e-01, 8.967902174e-02, 3.271355515e-02, 4.858559501e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.573233048e+00, 2.539436055e+00)
@test idx == [93, 94, 110, 111] + 1
@test w ≈ [3.289345075e-02, 9.525544479e-01, 7.761787977e-03, 6.790313325e-03]

(idx, w) = Healpix.getinterpolRing(resol, 1.556697970e+00, 2.991687932e+00)
@test idx == [79, 80, 95, 96] + 1
@test w ≈ [3.213984335e-02, 5.205554974e-02, 8.074916993e-01, 1.083129076e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.548287295e+00, 3.384326115e+00)
@test idx == [80, 81, 96, 97] + 1
@test w ≈ [5.133439831e-02, 8.308955185e-02, 7.633379430e-01, 1.022381068e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.585334905e+00, 3.714465380e+00)
@test idx == [96, 97, 113, 114] + 1
@test w ≈ [3.761514669e-02, 8.755604603e-01, 4.698863049e-02, 3.983576251e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.554031199e+00, 4.140423585e+00)
@test idx == [82, 83, 98, 99] + 1
@test w ≈ [4.570518921e-02, 5.441616345e-02, 8.607320542e-01, 3.914659313e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.583952510e+00, 4.560816716e+00)
@test idx == [99, 100, 115, 116] + 1
@test w ≈ [8.163656227e-01, 1.050656457e-01, 3.032561344e-02, 4.824311808e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.551037089e+00, 4.895336191e+00)
@test idx == [84, 85, 99, 100] + 1
@test w ≈ [6.302835410e-02, 5.497382001e-02, 3.010148585e-02, 8.518963400e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.558273267e+00, 5.335548510e+00)
@test idx == [85, 86, 101, 102] + 1
@test w ≈ [3.089759571e-02, 4.389012409e-02, 8.448458321e-01, 8.036644813e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.560260922e+00, 5.691659682e+00)
@test idx == [86, 87, 101, 102] + 1
@test w ≈ [3.185558091e-02, 3.106186045e-02, 5.910774370e-03, 9.311717843e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.588580220e+00, 6.085795343e+00)
@test idx == [102, 103, 119, 104] + 1
@test w ≈ [2.368034854e-03, 8.914265474e-01, 5.338409134e-02, 5.282132639e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.751093487e+00, 1.914423495e-02)
@test idx == [104, 105, 135, 120] + 1
@test w ≈ [8.803477800e-01, 4.511676098e-02, 3.363409650e-02, 4.090136249e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.783866814e+00, 3.766054576e-01)
@test idx == [104, 105, 120, 121] + 1
@test w ≈ [3.013624091e-02, 7.052154750e-01, 1.431699785e-01, 1.214783056e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.785396474e+00, 7.359002918e-01)
@test idx == [105, 106, 121, 122] + 1
@test w ≈ [9.156918496e-02, 6.349092209e-01, 1.712369071e-01, 1.022846870e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.773782348e+00, 1.202139487e+00)
@test idx == [107, 108, 122, 123] + 1
@test w ≈ [7.452481558e-01, 4.860193588e-02, 9.045382491e-02, 1.156960834e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.728573392e+00, 1.541054504e+00)
@test idx == [91, 92, 107, 108] + 1
@test w ≈ [3.325185754e-02, 2.450343987e-02, 7.136271976e-02, 8.708819828e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.756993254e+00, 1.918099916e+00)
@test idx == [108, 109, 124, 125] + 1
@test w ≈ [1.030262695e-01, 7.882146653e-01, 6.695193627e-02, 4.180712900e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.777406002e+00, 2.371869229e+00)
@test idx == [110, 111, 125, 126] + 1
@test w ≈ [7.419820463e-01, 3.084780886e-02, 1.045174864e-01, 1.226526584e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.753826594e+00, 2.746049052e+00)
@test idx == [110, 111, 126, 127] + 1
@test w ≈ [6.588770748e-03, 9.030214539e-01, 4.584962692e-02, 4.454014844e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.757740845e+00, 3.156422697e+00)
@test idx == [112, 113, 127, 128] + 1
@test w ≈ [8.534108749e-01, 3.349340379e-02, 5.227686913e-02, 6.081885218e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.730614741e+00, 3.559474066e+00)
@test idx == [96, 97, 113, 114] + 1
@test w ≈ [1.986030858e-02, 2.570405468e-02, 8.932312312e-01, 6.120440551e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.782426275e+00, 3.886387874e+00)
@test idx == [113, 114, 129, 130] + 1
@test w ≈ [7.689535583e-02, 6.668126941e-01, 1.546451647e-01, 1.016467853e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.762506035e+00, 4.322418883e+00)
@test idx == [115, 116, 130, 131] + 1
@test w ≈ [8.532909049e-01, 5.971271725e-03, 6.939088215e-02, 7.134694122e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.779515720e+00, 4.707569946e+00)
@test idx == [115, 116, 131, 132] + 1
@test w ≈ [9.333654487e-03, 7.512580603e-01, 1.226420582e-01, 1.167662271e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.718607381e+00, 5.079958700e+00)
@test idx == [100, 101, 116, 117] + 1
@test w ≈ [6.614058002e-02, 5.113174111e-02, 5.648697579e-02, 8.262407031e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.760590315e+00, 5.473089858e+00)
@test idx == [117, 118, 133, 134] + 1
@test w ≈ [5.473885856e-02, 8.156361027e-01, 7.296478357e-02, 5.666025512e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.743387337e+00, 5.868002046e+00)
@test idx == [118, 119, 134, 135] + 1
@test w ≈ [5.554737124e-02, 9.146193069e-01, 1.662478246e-02, 1.320853937e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.928481425e+00, 2.006267335e-01)
@test idx == [120, 121, 136, 137] + 1
@test w ≈ [8.930393748e-01, 9.833898656e-03, 4.750548018e-02, 4.962124640e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.899373338e+00, 6.305174618e-01)
@test idx == [105, 106, 121, 122] + 1
@test w ≈ [2.576100254e-02, 3.955586352e-02, 8.359810324e-01, 9.870210151e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.876684423e+00, 9.514637558e-01)
@test idx == [106, 107, 121, 122] + 1
@test w ≈ [1.136526786e-01, 8.327895165e-02, 6.193057798e-02, 7.411377918e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.906536951e+00, 1.361809664e+00)
@test idx == [107, 108, 122, 123] + 1
@test w ≈ [1.264561043e-02, 1.111628684e-02, 3.141550400e-02, 9.448225987e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.895621871e+00, 1.784303891e+00)
@test idx == [108, 109, 124, 125] + 1
@test w ≈ [3.973458041e-02, 4.734394571e-02, 8.730336089e-01, 3.988786498e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.951130318e+00, 2.175461836e+00)
@test idx == [125, 126, 141, 142] + 1
@test w ≈ [7.486178835e-01, 3.100406507e-02, 1.014250145e-01, 1.189530368e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.932190104e+00, 2.589998379e+00)
@test idx == [126, 127, 142, 143] + 1
@test w ≈ [7.985030933e-01, 8.418819333e-02, 4.746583722e-02, 6.984287614e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.925886100e+00, 2.950770916e+00)
@test idx == [127, 128, 143, 144] + 1
@test w ≈ [9.040885356e-01, 1.290804293e-02, 4.033331844e-02, 4.267010301e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.914059311e+00, 3.311174358e+00)
@test idx == [127, 128, 144, 145] + 1
@test w ≈ [6.689288470e-02, 9.144630143e-01, 1.059290208e-02, 8.051198916e-03]

(idx, w) = Healpix.getinterpolRing(resol, 1.960513132e+00, 3.690312886e+00)
@test idx == [128, 129, 145, 146] + 1
@test w ≈ [7.482000039e-02, 6.537423119e-01, 1.635942466e-01, 1.078434411e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.878165657e+00, 4.107800291e+00)
@test idx == [114, 115, 129, 130] + 1
@test w ≈ [1.016226607e-01, 8.671656954e-02, 3.211940862e-02, 7.795413612e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.923994634e+00, 4.549404339e+00)
@test idx == [131, 132, 147, 148] + 1
@test w ≈ [8.485042815e-01, 7.878532374e-02, 3.017750278e-02, 4.253289200e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.945528802e+00, 4.861343452e+00)
@test idx == [131, 132, 148, 149] + 1
@test w ≈ [9.777195174e-02, 7.123324725e-01, 1.178663901e-01, 7.202918565e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.879719029e+00, 5.256524395e+00)
@test idx == [117, 118, 132, 133] + 1
@test w ≈ [1.101740658e-01, 6.915429963e-02, 9.386066301e-02, 7.268109716e-01]

(idx, w) = Healpix.getinterpolRing(resol, 1.919161218e+00, 5.664599318e+00)
@test idx == [133, 134, 150, 151] + 1
@test w ≈ [7.172566552e-02, 8.818665489e-01, 2.669451412e-02, 1.971327150e-02]

(idx, w) = Healpix.getinterpolRing(resol, 1.895750284e+00, 6.082619283e+00)
@test idx == [119, 104, 134, 135] + 1
@test w ≈ [4.409379293e-02, 4.223983240e-02, 9.810206540e-03, 9.038561681e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.052081277e+00, 4.127474467e-02)
@test idx == [135, 120, 136, 137] + 1
@test w ≈ [9.093021785e-02, 1.393342347e-01, 6.888322816e-01, 8.090326580e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.091890932e+00, 4.075109167e-01)
@test idx == [120, 121, 137, 138] + 1
@test w ≈ [6.299634942e-03, 7.327621415e-03, 9.491687089e-01, 3.720403476e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.087951781e+00, 8.280528712e-01)
@test idx == [121, 122, 138, 139] + 1
@test w ≈ [1.372314997e-02, 2.134028244e-02, 8.601258123e-01, 1.048107553e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.059010906e+00, 1.166517057e+00)
@test idx == [122, 123, 138, 139] + 1
@test w ≈ [1.019554984e-01, 9.059912440e-02, 2.381052084e-02, 7.836348563e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.085650517e+00, 1.532296997e+00)
@test idx == [123, 124, 139, 140] + 1
@test w ≈ [2.845852568e-02, 1.912797906e-02, 9.337246528e-02, 8.590410300e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.130752169e+00, 1.986466754e+00)
@test idx == [141, 142, 156, 157] + 1
@test w ≈ [7.754412320e-01, 4.817850260e-02, 7.787258356e-02, 9.850768181e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.100018824e+00, 2.312082224e+00)
@test idx == [141, 142, 157, 158] + 1
@test w ≈ [1.092662877e-01, 8.634511606e-01, 1.670595116e-02, 1.057660055e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.105695556e+00, 2.742342278e+00)
@test idx == [142, 143, 158, 159] + 1
@test w ≈ [1.576814797e-02, 9.294095792e-01, 2.832572177e-02, 2.649655105e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.052799231e+00, 3.177930195e+00)
@test idx == [127, 128, 144, 145] + 1
@test w ≈ [9.223324726e-02, 1.341242238e-01, 7.020552275e-01, 7.158730143e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.072951409e+00, 3.486452912e+00)
@test idx == [128, 129, 144, 145] + 1
@test w ≈ [7.256200572e-02, 4.413083161e-02, 1.076049752e-01, 7.757021875e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.137101401e+00, 3.878478789e+00)
@test idx == [145, 146, 161, 162] + 1
@test w ≈ [9.794060269e-02, 6.948768907e-01, 1.291855164e-01, 7.799699013e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.088561403e+00, 4.293045209e+00)
@test idx == [130, 131, 146, 147] + 1
@test w ≈ [1.802695460e-02, 1.371901701e-02, 6.569617686e-02, 9.025578515e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.074352685e+00, 4.705936370e+00)
@test idx == [131, 132, 147, 148] + 1
@test w ≈ [5.632580050e-02, 5.274153456e-02, 1.463930336e-02, 8.762933616e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.094719000e+00, 5.114944728e+00)
@test idx == [149, 150, 164, 165] + 1
@test w ≈ [9.733683111e-01, 2.506035359e-02, 7.462274571e-04, 8.251078410e-04]

(idx, w) = Healpix.getinterpolRing(resol, 2.079917751e+00, 5.481143315e+00)
@test idx == [133, 134, 149, 150] + 1
@test w ≈ [4.273069270e-02, 3.605252878e-02, 3.904407903e-02, 8.821726995e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.115407077e+00, 5.885408296e+00)
@test idx == [150, 151, 166, 167] + 1
@test w ≈ [1.161272165e-02, 8.864511760e-01, 5.228617071e-02, 4.964993161e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.258871960e+00, 1.803457322e-01)
@test idx == [136, 137, 167, 152] + 1
@test w ≈ [1.092688788e-01, 9.279898808e-02, 3.251841880e-02, 7.654137143e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.289640532e+00, 5.967327945e-01)
@test idx == [137, 138, 153, 154] + 1
@test w ≈ [2.536647291e-02, 2.743278049e-02, 9.286663171e-01, 1.853442952e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.325980904e+00, 9.936909561e-01)
@test idx == [154, 155, 169, 170] + 1
@test w ≈ [8.567258866e-01, 2.687310934e-02, 7.009552866e-02, 4.630547536e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.273799548e+00, 1.390628538e+00)
@test idx == [139, 140, 155, 156] + 1
@test w ≈ [5.948218963e-02, 7.016696349e-02, 8.344867408e-01, 3.586410612e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.319259993e+00, 1.787519828e+00)
@test idx == [156, 157, 170, 171] + 1
@test w ≈ [8.668930048e-01, 4.743715209e-02, 7.375197083e-03, 7.829464605e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.309119850e+00, 2.182197401e+00)
@test idx == [157, 158, 171, 172] + 1
@test w ≈ [9.060128129e-01, 5.468284522e-02, 1.306122438e-02, 2.624311754e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.265683831e+00, 2.574990418e+00)
@test idx == [142, 143, 158, 159] + 1
@test w ≈ [7.484947689e-02, 9.417172969e-02, 7.834806639e-01, 4.749812951e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.294276700e+00, 2.942604403e+00)
@test idx == [143, 144, 158, 159] + 1
@test w ≈ [1.535747737e-02, 1.495017771e-02, 6.515768992e-03, 9.631765759e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.312479071e+00, 3.361013294e+00)
@test idx == [160, 161, 173, 174] + 1
@test w ≈ [8.897971795e-01, 5.553854254e-02, 4.424387007e-03, 5.023989091e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.319902347e+00, 3.750943468e+00)
@test idx == [161, 162, 174, 175] + 1
@test w ≈ [8.642748091e-01, 4.711820527e-02, 2.979194306e-02, 5.881504256e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.323472566e+00, 4.130766233e+00)
@test idx == [162, 163, 175, 176] + 1
@test w ≈ [8.781426996e-01, 1.692559507e-02, 6.409413480e-02, 4.083757054e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.277479067e+00, 4.559244024e+00)
@test idx == [147, 148, 163, 164] + 1
@test w ≈ [4.359925865e-02, 6.819932214e-02, 7.904818813e-01, 9.771953793e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.293845635e+00, 4.918323088e+00)
@test idx == [148, 149, 164, 165] + 1
@test w ≈ [1.540869028e-02, 1.699020331e-02, 9.439849643e-01, 2.361614212e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.254698294e+00, 5.255534892e+00)
@test idx == [149, 150, 164, 165] + 1
@test w ≈ [1.371444054e-01, 8.517130785e-02, 9.090374388e-02, 6.867805428e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.335132718e+00, 5.705440592e+00)
@test idx == [166, 167, 178, 179] + 1
@test w ≈ [8.175226144e-01, 2.422998539e-02, 9.548824229e-02, 6.275915793e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.269487300e+00, 6.077566661e+00)
@test idx == [151, 136, 166, 167] + 1
@test w ≈ [7.883863050e-02, 7.173068026e-02, 2.004960718e-02, 8.293810821e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.520579872e+00, 2.648638089e-01)
@test idx == [168, 169, 187, 180] + 1
@test w ≈ [9.877668141e-01, 5.815050238e-03, 1.044647365e-03, 5.373488285e-03]

(idx, w) = Healpix.getinterpolRing(resol, 2.541450067e+00, 7.756938929e-01)
@test idx == [168, 169, 180, 181] + 1
@test w ≈ [1.658364704e-02, 8.781953334e-01, 5.391060606e-02, 5.131041348e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.514913813e+00, 1.331701892e+00)
@test idx == [154, 155, 170, 171] + 1
@test w ≈ [2.145299695e-03, 1.756365586e-02, 9.377824236e-01, 4.250862084e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.543821469e+00, 1.796433428e+00)
@test idx == [170, 171, 181, 182] + 1
@test w ≈ [6.102243982e-02, 8.225299396e-01, 2.476956413e-02, 9.167805646e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.482294580e+00, 2.361156715e+00)
@test idx == [157, 158, 172, 173] + 1
@test w ≈ [8.229596618e-02, 8.656344913e-02, 8.232637402e-01, 7.876844484e-03]

(idx, w) = Healpix.getinterpolRing(resol, 2.559458941e+00, 2.864864270e+00)
@test idx == [172, 173, 183, 184] + 1
@test w ≈ [2.308132040e-02, 7.864407326e-01, 1.623522630e-01, 2.812568397e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.562069696e+00, 3.443195446e+00)
@test idx == [174, 175, 183, 184] + 1
@test w ≈ [7.365629126e-01, 6.059940420e-02, 2.352661520e-02, 1.793110680e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.566192249e+00, 3.889603536e+00)
@test idx == [174, 175, 184, 185] + 1
@test w ≈ [5.552734546e-02, 7.221181356e-01, 1.217619937e-01, 1.005925252e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.506202743e+00, 4.478563018e+00)
@test idx == [162, 163, 176, 177] + 1
@test w ≈ [5.682081792e-03, 5.385797669e-02, 8.902155794e-01, 5.024436208e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.487185853e+00, 4.961375377e+00)
@test idx == [164, 165, 176, 177] + 1
@test w ≈ [1.268583310e-01, 1.963588761e-02, 2.088614928e-02, 8.326196321e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.486833278e+00, 5.453248632e+00)
@test idx == [165, 166, 177, 178] + 1
@test w ≈ [9.085086684e-02, 5.725549043e-02, 7.246402548e-02, 7.794296172e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.507106131e+00, 5.974733165e+00)
@test idx == [166, 167, 178, 179] + 1
@test w ≈ [1.581753547e-02, 3.959180913e-02, 8.416321367e-02, 8.604274417e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.751429923e+00, 4.135303144e-01)
@test idx == [180, 181, 191, 188] + 1
@test w ≈ [8.746718774e-01, 2.383112936e-02, 2.402823835e-02, 7.746875487e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.689999731e+00, 1.182371359e+00)
@test idx == [169, 170, 181, 182] + 1
@test w ≈ [4.631682676e-02, 1.452039836e-01, 8.040794700e-01, 4.399719682e-03]

(idx, w) = Healpix.getinterpolRing(resol, 2.770597294e+00, 2.006543154e+00)
@test idx == [182, 183, 188, 189] + 1
@test w ≈ [7.615901930e-01, 4.416343002e-02, 4.323826337e-02, 1.510081136e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.731264459e+00, 2.720661427e+00)
@test idx == [182, 183, 189, 190] + 1
@test w ≈ [3.580544847e-02, 9.602766329e-01, 3.008856255e-03, 9.090623442e-04]

(idx, w) = Healpix.getinterpolRing(resol, 2.706636709e+00, 3.565204154e+00)
@test idx == [174, 175, 184, 185] + 1
@test w ≈ [7.791190867e-02, 3.484675314e-02, 8.523204828e-01, 3.492085540e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.704369805e+00, 4.359710081e+00)
@test idx == [175, 176, 185, 186] + 1
@test w ≈ [2.143389449e-02, 1.020566566e-01, 8.318466651e-01, 4.466278377e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.721431698e+00, 5.089689476e+00)
@test idx == [177, 178, 185, 186] + 1
@test w ≈ [3.329385225e-02, 9.422932613e-03, 1.876857953e-02, 9.385146356e-01]

(idx, w) = Healpix.getinterpolRing(resol, 2.757524566e+00, 5.850901197e+00)
@test idx == [186, 187, 191, 188] + 1
@test w ≈ [4.379924360e-02, 8.252122683e-01, 1.015423562e-01, 2.944613199e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.967997516e+00, 7.479336534e-01)
@test idx == [191, 188, 189, 190] + 1
@test w ≈ [5.800865797e-02, 8.664704311e-01, 3.776045547e-02, 3.776045547e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.960633628e+00, 2.366832475e+00)
@test idx == [189, 190, 191, 188] + 1
@test w ≈ [9.077348403e-01, 3.475060927e-02, 2.875727520e-02, 2.875727520e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.892099009e+00, 3.953399824e+00)
@test idx == [184, 185, 190, 191] + 1
@test w ≈ [1.015841648e-01, 1.162323202e-01, 7.690330575e-01, 1.315045749e-02]

(idx, w) = Healpix.getinterpolRing(resol, 2.949545649e+00, 5.494536152e+00)
@test idx == [190, 191, 188, 189] + 1
@test w ≈ [1.714477856e-02, 9.524532703e-01, 1.520097559e-02, 1.520097559e-02]

