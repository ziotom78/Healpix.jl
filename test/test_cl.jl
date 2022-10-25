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

import Random
rng = Random.seed!(1234) #we set the seed to guarantee reproducibility
cl = Vector(0:4)
alm = Healpix.synalm(cl, 4, rng)

ref_alm = [
    -0.0 + 0.0im,
    1.0872084924285859 + 0.0im,
    0.5943192583849796 + 0.0im,
    0.7996791530503793 + 0.0im,
    -1.8728302820268359 + 0.0im,
    -0.29669466345292583 + 0.5083460918445625im,
    -0.6856709022761192 + 2.054763056064037im,
    -1.034857629674229 + 1.0859956426035118im,
    -0.8736960340524973 + 0.21500322847968445im,
    0.3248927294469158 - 0.3049012551964323im,
    0.016624443835621037 + 0.1293962624603848im,
    -1.2719064982797272 + 1.151859436947903im,
    -0.662762592300875 - 0.8453680298536129im,
    0.4221137801332629 - 2.790353436891747im,
    0.4550092269865841 + 0.38238373487493016im
]

@test isapprox(alm.alm, ref_alm)
