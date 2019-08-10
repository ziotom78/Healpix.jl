################################################################################

function pixelWindowPath(datadir, nside)
    joinpath(datadir, @sprintf("pixel_window_n%04d.fits", nside))
end

function readPixelWindowT(fileName, nside)
    f = FITSIO.fits_open_table(fileName)
    try
        pixwin = Array(Float64, FITSIO.fits_get_num_rows(f))
        FITSIO.fits_read_col(f, 1, 1, 1, pixwin)

        return pixwin
    finally
        FITSIO.fits_close_file(f)
    end
end

function readPixelWindowP(fileName, nside)
    local pixwinT
    local pixwinP

    f = FITSIO.fits_open_table(fileName)
    try
        pixwinT = Array(Float64, FITSIO.fits_get_num_rows(f))
        pixwinP = Array(Float64, FITSIO.fits_get_num_rows(f))
        FITSIO.fits_read_col(f, 1, 1, 1, pixwinT)
        FITSIO.fits_read_col(f, 2, 1, 1, pixwinP)

        return (pixwinT, pixwinP)
    finally
        FITSIO.fits_close_file(f)
    end
end
