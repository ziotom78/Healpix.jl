##test readClFromFITS/writeClToFITS

C_l_test = Vector{Float64}(0:20)
Healpix.writeClToFITS("Cltest.fits", C_l_test)
C_l_read = Healpix.readClFromFITS("Cltest.fits", Float64)

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
