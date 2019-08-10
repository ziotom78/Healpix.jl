m = Healpix.Map{Float64,Healpix.RingOrder}(1)
m.pixels = 1.0:12.0

# Do not run @test here, just check that the function can be ran
bmp = Healpix.equirectangular(m, Dict(:width => 50))
bmp = Healpix.mollweide(m, Dict(:width => 50))
bmp = Healpix.orthographic(m, Dict(:width => 50))
bmp = Healpix.gnomonic(m, Dict(:width => 50))
bmp = Healpix.orthographic2(m, Dict(:width => 50))
