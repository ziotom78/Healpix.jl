m = Healpix.HealpixMap{Float64,Healpix.RingOrder}(1)
m.pixels = 1.0:12.0

let projections = [
    (Healpix.equiproj, Healpix.equiprojinv, "Equirectangular"),
    (Healpix.mollweideproj, Healpix.mollweideprojinv, "Mollweide"),
    ]
    @testset "$name projection" for (idx, name) in enumerate((x[3] for x in projections))
        @testset "Point ($x_in, $y_in)" for (x_in, y_in) in [
            (0.0, 0.0),
            (0.1, 0.1),
            (0.9, 0.9),
            (0.1, 0.9),
            (0.9, 0.1),
            (0.5, 0.1),
            (0.5, 0.9),
            (0.1, 0.5),
            (0.9, 0.5),
            ]

            projfn = projections[idx][1]
            invprojfn = projections[idx][2]
            
            visible_inv, lat, long = invprojfn(x_in, y_in)
            (!visible_inv) && continue
            visible, x, y = projfn(lat, long)
            
            @test visible_inv == visible

            @test x ≈ x_in
            @test y ≈ y_in
        end
    end
end

# Do not run @test here, just check that the function can be ran
bmp = Healpix.equirectangular(m, Dict(:width => 50))
bmp = Healpix.mollweide(m, Dict(:width => 50))
bmp = Healpix.orthographic(m, Dict(:width => 50))
bmp = Healpix.gnomonic(m, Dict(:width => 50))
bmp = Healpix.orthographic2(m, Dict(:width => 50))
