import CFITSIO

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

alm = Healpix.readAlmFromFITS("alm.fits", ComplexF64).alm
@test alm[1] ≈ (5.443205775735018e+03 + 0.000000000000e+00im) atol = eps
@test alm[2] ≈ (-3.143659646588944e+03 + 0.000000000000e+00im) atol = eps
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
@test alm[14] ≈ (-6.895464505702813e-01 + 1.5972357769351044e+01im) atol = eps
@test alm[15] ≈ (4.344648310121e-15 + -1.892253965354e-15im) atol = eps
@test alm[16] ≈ (6.488192068017e-02 + 3.967268430232e+00im) atol = eps
@test alm[17] ≈ (1.340722726247e-15 + -9.822172672993e-15im) atol = eps
@test alm[18] ≈ (1.485579846001e-01 + 7.942515768792e-01im) atol = eps
@test alm[19] ≈ (6.938894269949477e-01 + -1.0288216337491448e+01im) atol = eps
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


## test gaussbeam

# spin-0
bl_0 = Healpix.gaussbeam(deg2rad(1), 512; pol=false)
@test bl_0[50+1] ≈ 0.9323560235080715
@test bl_0[500+1] ≈ 0.0010276785784086898

# spin-2
bl_2 = Healpix.gaussbeam(deg2rad(1), 512; pol=true)
@test bl_2[50+1] ≈ 0.9324584647703964
@test bl_2[500+1] ≈ 0.001027791493098124

## test writeAlmToFITS & readAlmFromFITS:
testalm = Healpix.Alm(2, 2, ComplexF64.(1:6))
file_name = tempname() * "testalm.fits"
Healpix.writeAlmToFITS(file_name, testalm)
almread = Healpix.readAlmFromFITS(file_name, ComplexF64)

@test almread.alm == testalm.alm

## test almxfl & almxfl! through '*'

alm = Healpix.Alm(3,3)
#let's fill the alms as a_ℓm = ℓ + m
for ℓ in 0:3
    for m in 0:ℓ
        alm.alm[Healpix.almIndex(alm,ℓ,m)] = ℓ + m
    end
end

alm_test = Healpix.Alm(3,3)
#let's fill the alms as a_ℓm = (ℓ + m)*ℓ
for ℓ in 0:3
    for m in 0:ℓ
        alm_test.alm[Healpix.almIndex(alm_test,ℓ,m)] = (ℓ + m)*ℓ
    end
end

C_l = Vector{Float64}(0:5)
almnew = alm*C_l
@test alm_test.alm == almnew.alm
@test alm_test.alm != alm.alm

## test almExplicitIndex
#let's see if the function reproduces the indexing shown on the file "alm.fits"
alm = Healpix.readAlmFromFITS("alm.fits", ComplexF64)
f = CFITSIO.fits_open_table("alm.fits")
numOfRows = CFITSIO.fits_get_num_rows(f)
idx = Array{Int64}(undef, numOfRows)
CFITSIO.fits_read_col(f, 1, 1, 1, idx)

idx_test_1 = Healpix.almExplicitIndex(alm)
idx_test_2 = Healpix.almExplicitIndex(alm.lmax, alm.mmax)

@test idx == idx_test_1
@test idx == idx_test_2

## test each_ell

alm = Healpix.Alm(5,5)
ell = Healpix.each_ell(alm, [0,1,2])
test_ell = [0, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 2, 3, 4, 5]

@test ell == test_ell

## test each_ell_idx

idx = Healpix.each_ell_idx(alm, [0,2])
test_idx = [1, 2, 3, 4, 5, 6, 12, 13, 14, 15]

@test idx == test_idx

## test each_m

ms = Healpix.each_m(alm, [2, 1, 0])
test_ms = [0, 1, 2, 0, 1, 0]

@test ms == test_ms

## test each_m_idx

idx = Healpix.each_m_idx(alm, [0, 3])
test_idx = [1, 4, 9, 13, 16]

@test idx == test_idx

## test each_ell_m

ellm = Healpix.each_ell_m(alm)
test_ellm = [(0, 0), (1, 0), (2, 0), (3, 0), (4, 0), (5, 0), (1, 1), (2, 1),
            (3, 1), (4, 1), (5, 1), (2, 2), (3, 2), (4, 2), (5, 2), (3, 3),
            (4, 3), (5, 3), (4, 4), (5, 4), (5, 5)]

@test ellm == test_ellm

## test alm algebra

alm1 = Healpix.Alm(3,3, ones(ComplexF64, Healpix.numberOfAlms(3)) .+ 1.00im)
alm2 = Healpix.Alm(3,3, 2*ones(ComplexF64, Healpix.numberOfAlms(3)) .+ 1.00im)

#+
test_alm = [3.0 + 2.0im, 3.0 + 2.0im, 3.0 + 2.0im,
            3.0 + 2.0im, 3.0 + 2.0im, 3.0 + 2.0im,
            3.0 + 2.0im, 3.0 + 2.0im, 3.0 + 2.0im,
            3.0 + 2.0im]

@test isapprox(test_alm, (alm1 + alm2).alm)

#-
test_alm = [1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im,
            1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im,
            1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im,
            1.0 + 0.0im]

@test isapprox(test_alm, (alm2 - alm1).alm)

#alm*alm
test_alm = [2.0 + 0.0im, 2.0 + 0.0im, 2.0 + 0.0im,
            2.0 + 0.0im, 1.0 + 3.0im, 1.0 + 3.0im,
            1.0 + 3.0im, 1.0 + 3.0im, 1.0 + 3.0im,
            1.0 + 3.0im]

@test isapprox(test_alm, (alm1 * alm2).alm)

#alm*c
c = 2.0

test_alm = [2.0 + 2.0im, 2.0 + 2.0im, 2.0 + 2.0im,
            2.0 + 2.0im, 2.0 + 2.0im, 2.0 + 2.0im,
            2.0 + 2.0im, 2.0 + 2.0im, 2.0 + 2.0im,
            2.0 + 2.0im]

@test isapprox(test_alm, (alm1 * c).alm)


#alm/alm
test_alm = [0.5 + 0.0im, 0.5 + 0.0im, 0.5 + 0.0im,
            0.5 + 0.0im, 0.6 + 0.2im, 0.6 + 0.2im,
            0.6 + 0.2im, 0.6 + 0.2im, 0.6 + 0.2im,
            0.6 + 0.2im]

@test isapprox(test_alm, (alm1/alm2).alm)

#alm/c
test_alm = [0.5 + 0.5im, 0.5 + 0.5im, 0.5 + 0.5im,
            0.5 + 0.5im, 0.5 + 0.5im, 0.5 + 0.5im,
            0.5 + 0.5im, 0.5 + 0.5im, 0.5 + 0.5im,
            0.5 + 0.5im]

@test isapprox(test_alm, (alm1/c).alm)

#dot
test_res = 44.0
@test isapprox(test_res, LinearAlgebra.:dot(alm1, alm2))
