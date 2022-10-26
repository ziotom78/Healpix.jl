##test readClFromFITS/writeClToFITS

C_l_test = Vector{Float64}(0:20)
file_name = tempname() * ".fits"
Healpix.writeClToFITS(file_name, C_l_test)
C_l_read = Healpix.readClFromFITS(file_name, Float64)

@test C_l_test == C_l_read

## test cl2dl/dl2cl

### with monopole
C_l_test = Vector{Float64}(randn(10))
D_l_test = Healpix.cl2dl(C_l_test, 0)
C_l_conv = Healpix.dl2cl(D_l_test, 0)

#we expect the monopole component to be 0
@test isapprox(C_l_conv[1], 0.0)
@test isapprox(C_l_conv[2:end], C_l_test[2:end])

### from dipole on
C_l_test_2 = Vector{Float64}(randn(10))
D_l_test_2 = Healpix.cl2dl(C_l_test_2, 1)
C_l_conv_2 = Healpix.dl2cl(D_l_test_2, 0)

#we expect the monopole component to be 0
@test isapprox(C_l_conv_2[1], 0.0)
@test isapprox(C_l_conv_2[2:end], C_l_test_2)

### with no monopole and no dipole
C_l_test_3 = Vector{Float64}(randn(12))
D_l_test_3 = Healpix.cl2dl(C_l_test_3, 2)
C_l_conv_3 = Healpix.dl2cl(D_l_test_3, 0) #here the PS has been 'completed' by cl2dl

#we expect the monopole component to be 0
@test isapprox(C_l_conv_3[1], 0.0)
@test isapprox(C_l_conv_3[2], 0.0)
@test isapprox(C_l_conv_3[3:end], C_l_test_3)

### with no monopole and no dipole - STARTING FROM Dl's
D_l_test_4 = Vector{Float64}(randn(12))
C_l_test_4 = Healpix.dl2cl(D_l_test_4, 2)
D_l_conv_4 = Healpix.cl2dl(C_l_test_4, 0) #here the PS has been 'completed' by cl2dl

#we expect the monopole component to be 0
@test isapprox(D_l_conv_4[1], 0.0)
@test isapprox(D_l_conv_4[2], 0.0)
@test isapprox(D_l_conv_4[3:end], D_l_test_4)

## test synalm

cl = Vector(0:4)
alm = Healpix.synalm(cl, 4)

@test alm.alm[1] == 0.0 + 0.0im #\ell=0 should be 0 because of our cl's
@test imag.(alm.alm[1:5]) == zeros(Float64, 5) #field should be always real
@test count(i->(i==0.00000000000e+00), real.(alm.alm[2:end])) == 0 #none of the other alms should be exactly 0.
@test count(i->(i==0.00000000000e+00), imag.(alm.alm[6:end])) == 0 # "

#other methods of the function
alm = Healpix.synalm(cl)

@test alm.alm[1] == 0.0 + 0.0im #\ell=0 should be 0 because of our cl's
@test imag.(alm.alm[1:5]) == zeros(Float64, 5) #field should be always real
@test count(i->(i==0.00000000000e+00), real.(alm.alm[2:end])) == 0 #none of the other alms should be exactly 0.
@test count(i->(i==0.00000000000e+00), imag.(alm.alm[6:end])) == 0

import Random
rng = Random.seed!(1234)
alm = Healpix.synalm(cl, rng)

@test alm.alm[1] == 0.0 + 0.0im #\ell=0 should be 0 because of our cl's
@test imag.(alm.alm[1:5]) == zeros(Float64, 5) #field should be always real
@test count(i->(i==0.00000000000e+00), real.(alm.alm[2:end])) == 0 #none of the other alms should be exactly 0.
@test count(i->(i==0.00000000000e+00), imag.(alm.alm[6:end])) == 0

alm = Healpix.synalm(cl, 4, rng)

@test alm.alm[1] == 0.0 + 0.0im #\ell=0 should be 0 because of our cl's
@test imag.(alm.alm[1:5]) == zeros(Float64, 5) #field should be always real
@test count(i->(i==0.00000000000e+00), real.(alm.alm[2:end])) == 0 #none of the other alms should be exactly 0.
@test count(i->(i==0.00000000000e+00), imag.(alm.alm[6:end])) == 0

alm = Healpix.synalm(cl, 4, 4, rng)

@test alm.alm[1] == 0.0 + 0.0im #\ell=0 should be 0 because of our cl's
@test imag.(alm.alm[1:5]) == zeros(Float64, 5) #field should be always real
@test count(i->(i==0.00000000000e+00), real.(alm.alm[2:end])) == 0 #none of the other alms should be exactly 0.
@test count(i->(i==0.00000000000e+00), imag.(alm.alm[6:end])) == 0

## test spin-0 anafast

nside = 4
map = Healpix.HealpixMap{Float64,Healpix.RingOrder}(2 .* ones(Healpix.nside2npix(nside)))
test_cl = Healpix.anafast(map)

test_cl_ref = [
    50.25609318985963,
    2.083026679225843e-34,
    4.7018489094583055e-7,
    4.066000172958629e-34,
    4.606958743985832e-7,
    7.050822195517214e-34,
    4.302032807806699e-7,
    1.0729590299934532e-33,
    5.800538451198812e-7,
    5.8544687237680685e-34,
    8.960148278384837e-7,
    1.763468998685261e-33,
]

@test isapprox(test_cl, test_cl_ref)

## test anafast cross spectrum

map2 = Healpix.HealpixMap{Float64,Healpix.RingOrder}(3 .* ones(Healpix.nside2npix(nside)))
test_ccl = Healpix.anafast(map, map2)

test_ccl_ref = [
    75.38413978478945
    3.890609325797566e-34
    7.052773364188918e-7
    6.508566039307384e-34
    6.910438115978033e-7
    4.127231154504506e-34
    6.453049211709552e-7
    1.8179780948330372e-34
    8.700807676800265e-7
    -5.07503235735685e-35
    1.344022241757911e-6
    1.0083830113062526e-33
]

@test isapprox(test_ccl, test_ccl_ref)
