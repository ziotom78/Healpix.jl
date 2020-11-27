@test Healpix.numberOfAlms(10, 5) == 51
@test Healpix.numberOfAlms(10, 7) == 60
@test Healpix.numberOfAlms(12, 7) == 76
@test Healpix.numberOfAlms(12, 12) == 91
@test_throws DomainError(-1, "`lmax` is not positive or zero") Healpix.numberOfAlms(-1, 1)
@test_throws DomainError(-1, "`mmax` is not positive or zero") Healpix.numberOfAlms(4, -1)
@test_throws DomainError((5, 7), "`lmax` and `mmax` are inconsistent") Healpix.numberOfAlms(
    5,
    7,
)

alm = Healpix.Alm(10, 8)
@test Healpix.almIndex(alm, 4, 2) == 24
@test Healpix.almIndex(alm, 5, 2) == 25
@test Healpix.almIndex(alm, 5, 3) == 33
@test Healpix.almIndex(alm, [4, 6, 5], [3, 4, 5]) == [32, 41, 46]

alm = Healpix.Alm{ComplexF32}(10, 8)
@test Healpix.almIndex(alm, 4, 2) == 24
@test Healpix.almIndex(alm, 5, 2) == 25
@test Healpix.almIndex(alm, 5, 3) == 33
@test Healpix.almIndex(alm, [4, 6, 5], [3, 4, 5]) == [32, 41, 46]

alm = Healpix.readAlmFromFITS("alm.fits", ComplexF64)
@test alm[1] ≈ (5.443205775735e+03 + 0.000000000000e+00im) atol = eps
@test alm[2] ≈ (-3.143659646589e+03 + 0.000000000000e+00im) atol = eps
@test alm[3] ≈ (-8.445976910202e-07 + 0.000000000000e+00im) atol = eps
@test alm[4] ≈ (3.003475555079e-07 + 0.000000000000e+00im) atol = eps
@test alm[5] ≈ (-1.094164444296e-06 + 0.000000000000e+00im) atol = eps
@test alm[6] ≈ (3.745732939005e-07 + 0.000000000000e+00im) atol = eps
@test alm[7] ≈ (-1.344818023454e-06 + 0.000000000000e+00im) atol = eps
@test alm[8] ≈ (6.658742467775e-01 + -3.280201017201e+01im) atol = eps
@test alm[9] ≈ (1.200156696497e-15 + -1.483545355539e-15im) atol = eps
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

## test alm2cl
testalm = Healpix.Alm(2, 2, ComplexF64.(1:6))
@test isapprox(Healpix.alm2cl(testalm), [1.0, 12.0, 26.2])
