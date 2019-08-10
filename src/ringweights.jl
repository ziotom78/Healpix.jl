################################################################################

function ringWeightPath(datadir, nside)
    joinpath(datadir, @sprintf("weight_ring_n%05d.fits", nside))
end

function readWeightRing(fileName, nside)
    f = FITSIO.fits_open_table(fileName)
    try
        weights = Array(Float64, 2 * nside)
        FITSIO.fits_read_col(f, 1, 1, 1, weights)

        return weights
    finally
        FITSIO.fits_close_file(f)
    end
end
