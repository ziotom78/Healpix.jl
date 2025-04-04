################################################################################

function pixelWindowPath(datadir, nside)
    joinpath(datadir, @sprintf("pixel_window_n%04d.fits", nside))
end

function readPixelWindowT(fileName, nside)
    f = CFITSIO.fits_open_table(fileName)
    try
        pixwin = Vector(Float64, CFITSIO.fits_get_num_rows(f))
        CFITSIO.fits_read_col(f, 1, 1, 1, pixwin)

        return pixwin
    finally
        CFITSIO.fits_close_file(f)
    end
end

function readPixelWindowP(fileName, nside)
    local pixwinT
    local pixwinP

    f = CFITSIO.fits_open_table(fileName)
    try
        pixwinT = Vector(Float64, CFITSIO.fits_get_num_rows(f))
        pixwinP = Vector(Float64, CFITSIO.fits_get_num_rows(f))
        CFITSIO.fits_read_col(f, 1, 1, 1, pixwinT)
        CFITSIO.fits_read_col(f, 2, 1, 1, pixwinP)

        return (pixwinT, pixwinP)
    finally
        CFITSIO.fits_close_file(f)
    end
end
