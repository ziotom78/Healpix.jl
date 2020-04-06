import Healpix
using Test

## spin 0 map2alm
nside = 4
lmax = 4
map = Healpix.Map{Float64, Healpix.RingOrder}(2 .* ones(Healpix.nside2npix(nside)))
alm = Healpix.map2alm(map, lmax=lmax, mmax=lmax, niter=0)

test_alm_spin0 = [
    7.08981540e+00,  6.26606889e-17, -4.64452418e-02, -1.43573675e-16,
    -1.10275269e-01,  0.00000000e+00, 0.00000000e+00,  0.00000000e+00,
    0.00000000e+00,  0.00000000e+00, 0.00000000e+00,  0.00000000e+00,
    0.00000000e+00,  0.00000000e+00, -7.87873039e-04-9.64866195e-20im]
@test isapprox(alm.alm , test_alm_spin0)

# test kwarg defaults
@test isapprox(Healpix.map2alm(map).alm[1:2], 
    [7.08915321e+00+0.0im, -3.33069318e-17+0.0im])

## test type convert
map_int = Healpix.Map{Int64, Healpix.RingOrder}(
    2 .* ones(Int64,Healpix.nside2npix(nside)))
alm = Healpix.map2alm(map_int, lmax=lmax, mmax=lmax, niter=0)
@test isapprox(alm.alm , test_alm_spin0)

## spin 0 map2alm niter=3
map = Healpix.Map{Float64, Healpix.RingOrder}(2 .* ones(Healpix.nside2npix(nside)))
alm_niter_3 = Healpix.map2alm(map, lmax=lmax, mmax=lmax, niter=3).alm
alm_niter_0 = Healpix.map2alm(map, lmax=lmax, mmax=lmax, niter=0).alm
reference_Δ = [
    -4.72790441e-06+0.00000000e+00im,  6.60402068e-19+0.00000000e+00im,
    4.64356457e-02+0.00000000e+00im, -1.51317123e-18+0.00000000e+00im,
    1.10263459e-01+0.00000000e+00im, -3.88659647e-22-6.99765171e-23im,
    0.00000000e+00+0.00000000e+00im, -8.36466992e-22-1.50602326e-22im,
    0.00000000e+00+0.00000000e+00im,  3.21098183e-23+5.21377677e-22im,
    0.00000000e+00+0.00000000e+00im,  1.00695080e-22+1.63501911e-21im,
    0.00000000e+00+0.00000000e+00im,  0.00000000e+00+0.00000000e+00im,
    7.87761619e-04+3.09456747e-20im]
@test isapprox(alm_niter_3 .- alm_niter_0, reference_Δ)

## spin 0 alm2map
nside = 2
lmax = 2
nalms = Healpix.numberOfAlms(lmax, lmax)
alm = Healpix.Alm{ComplexF64}(lmax, lmax, 2 .* ones(ComplexF64, nalms))
map = Healpix.alm2map(alm, nside)

test_map_spin0 = [
        1.22822792,  3.61032581,  3.61032581,  1.22822792, -0.33740788,
       -0.16286106,  1.80075965,  4.40319186,  4.40319186,  1.80075965,
       -0.16286106, -0.33740788, -0.4312723 , -1.13862492, -0.90401688,
        2.07742993,  4.11691608,  2.07742993, -0.90401688, -1.13862492,
       -0.25082501, -1.68800153, -0.63028243,  2.30273478,  2.30273478,
       -0.63028243, -1.68800153, -0.25082501,  0.859566  , -0.41667555,
       -1.5554869 ,  0.05254053,  1.52313774,  0.05254053, -1.5554869 ,
       -0.41667555,  1.19694074, -0.29055765, -0.67742383,  0.26296318,
        0.26296318, -0.67742383, -0.29055765,  1.19694074,  1.03769816,
        0.21777049,  0.21777049,  1.03769816]
@test isapprox(map.pixels, test_map_spin0)

## test type convert
alm_bf = Healpix.Alm{BigFloat}(lmax, lmax, 2 .* ones(BigFloat, nalms))
map_bf = Healpix.alm2map(alm, nside)
@test isapprox(map_bf.pixels, test_map_spin0)

## spin 2 map2alm
nside = 4
lmax = 4
map = Healpix.PolarizedMap{Float64, Healpix.RingOrder}(
    2 .* ones(Healpix.nside2npix(nside)),
    2 .* ones(Healpix.nside2npix(nside)),
    2 .* ones(Healpix.nside2npix(nside)))
alms = Healpix.map2alm(map, lmax=lmax, mmax=lmax, niter=0)
test_alm_spin2 = [ 
    0.00000000e+00,  0.00000000e+00,
    -6.49104757e+00, -1.31064234e-16,
    -2.27889895e+00,  0.00000000e+00,
    0.00000000e+00,  0.00000000e+00,
    0.00000000e+00,  0.00000000e+00,
    0.00000000e+00,  0.00000000e+00,
    0.00000000e+00,  0.00000000e+00,
    2.36716190e-02+2.89893725e-18im]
@test isapprox(alms[1].alm, test_alm_spin0)
@test isapprox(alms[2].alm, test_alm_spin2)
@test isapprox(alms[3].alm, test_alm_spin2)

## test type convert
map_int = Healpix.PolarizedMap{BigFloat, Healpix.RingOrder}(
    2 .* ones(BigFloat, Healpix.nside2npix(nside)),
    2 .* ones(BigFloat, Healpix.nside2npix(nside)),
    2 .* ones(BigFloat, Healpix.nside2npix(nside)))
alms_int = Healpix.map2alm(map_int, lmax=lmax, mmax=lmax, niter=0)
@test isapprox(alms_int[1].alm, test_alm_spin0)
@test isapprox(alms_int[2].alm, test_alm_spin2)
@test isapprox(alms_int[3].alm, test_alm_spin2)

## spin 2 map2alm niter=3
alm_niter_3 = Healpix.map2alm(map, lmax=lmax, mmax=lmax, niter=3)[2].alm
alm_niter_0 = Healpix.map2alm(map, lmax=lmax, mmax=lmax, niter=0)[2].alm
reference_Δ = [ 
    0.00000000e+00+0.00000000e+00im,  0.00000000e+00+0.00000000e+00im,
    2.16700759e-02+0.00000000e+00im, -4.93711681e-17+0.00000000e+00im,
    4.23910198e-02+0.00000000e+00im,  0.00000000e+00+0.00000000e+00im,
    -7.22390428e-18-4.13077515e-18im,  7.75707020e-18-1.35791630e-17im,
    -4.18018942e-18-2.39951192e-18im,  3.66933589e-18+6.33078775e-18im,
    3.66577068e-18-2.13953406e-18im, -9.62096654e-19-1.59215462e-18im,
    1.79413778e-17+4.01876168e-17im,  4.51299239e-17-1.96930366e-17im,
    -3.38942495e-03+5.32175148e-18im]
@test isapprox(alm_niter_3 .- alm_niter_0, reference_Δ)

## spin 2 alm2map
nside = 2
lmax = 2
nalms = Healpix.numberOfAlms(lmax, lmax)
alm_t = Healpix.Alm{ComplexF64}(lmax, lmax, 2 .* ones(ComplexF64, nalms))
alm_e = Healpix.Alm{ComplexF64}(lmax, lmax, 2 .* ones(ComplexF64, nalms))
alm_b = Healpix.Alm{ComplexF64}(lmax, lmax, 2 .* ones(ComplexF64, nalms))

maps = Healpix.alm2map([alm_t, alm_e, alm_b], nside)

@test isapprox(maps.i , test_map_spin0)

test_map_spin2_q = [
    -1.9631492 ,  1.00333301, -0.59650858,  1.06275217, -2.6071711 ,
    -1.4882688 ,  0.1809384 , -0.25943678, -0.72916617,  0.72899969,
    1.43862464, -0.69806834, -1.78405186, -2.22862401, -1.17525562,
    -0.82688372, -0.99110781,  0.01416045,  1.20357653,  0.29450851,
    -1.70135994, -1.49205262, -1.49205262, -1.70135994, -0.73579893,
    0.83901787,  0.83901787, -0.73579893, -0.99110781, -0.82688372,
    -1.17525562, -2.22862401, -1.78405186,  0.29450851,  1.20357653,
    0.01416045, -0.25943678,  0.1809384 , -1.4882688 , -2.6071711 ,
    -0.69806834,  1.43862464,  0.72899969, -0.72916617,  1.00333301,
    -1.9631492 ,  1.06275217, -0.59650858]
@test isapprox(maps.q , test_map_spin2_q)

test_map_spin2_u = [ 
    1.06275217, -0.59650858,  1.00333301, -1.9631492 , -0.69806834,
     1.43862464,  0.72899969, -0.72916617, -0.25943678,  0.1809384 ,
    -1.4882688 , -2.6071711 , -1.78405186,  0.29450851,  1.20357653,
     0.01416045, -0.99110781, -0.82688372, -1.17525562, -2.22862401,
    -0.73579893,  0.83901787,  0.83901787, -0.73579893, -1.70135994,
    -1.49205262, -1.49205262, -1.70135994, -0.99110781,  0.01416045,
     1.20357653,  0.29450851, -1.78405186, -2.22862401, -1.17525562,
    -0.82688372, -0.72916617,  0.72899969,  1.43862464, -0.69806834,
    -2.6071711 , -1.4882688 ,  0.1809384 , -0.25943678, -0.59650858,
     1.06275217, -1.9631492 ,  1.00333301]
@test isapprox(maps.u , test_map_spin2_u)

## test type conversion
alm_t = Healpix.Alm{Complex{Float16}}(lmax, lmax, 2 .* ones(Complex{Float16}, nalms))
alm_e = Healpix.Alm{Complex{Float16}}(lmax, lmax, 2 .* ones(Complex{Float16}, nalms))
alm_b = Healpix.Alm{Complex{Float16}}(lmax, lmax, 2 .* ones(Complex{Float16}, nalms))
maps_float = Healpix.alm2map([alm_t, alm_e, alm_b], nside)

@test isapprox(maps_float.i , test_map_spin0)
@test isapprox(maps_float.q , test_map_spin2_q)
@test isapprox(maps_float.u , test_map_spin2_u)
