"""
    readClFromFITS{T <: Real}(f::CFITSIO.FITSFile, t::Type{T}) -> Vector{T}
    readClFromFITS{T <: Real}(fileName::String, t::Type{T}) -> Vector{T}
Read a set of C_ℓ coefficients from a FITS file.
"""
function readClFromFITS(f::CFITSIO.FITSFile, t::Type{T}) where {T <: Real}
    numOfRows = CFITSIO.fits_get_num_rows(f)

    Cl = Vector{T}(undef, numOfRows)

    CFITSIO.fits_read_col(f, 2, 1, 1, Cl)

    return Cl
end

function readClFromFITS(fileName, t::Type{T}) where {T <: Real}
    f = CFITSIO.fits_open_table(fileName)
    try
        result = readClFromFITS(f, t)
        return result
    finally
        CFITSIO.fits_close_file(f)
end
end


##############################################################################

"""
    writeClToFITS{T <: Real}(f::CFITSIO.FITSFile, Cl::Vector{T})
    writeClFromFITS{T <: Real}(fileName::String, Cl::Vector{T})
Write a set of C_ℓ coefficients to a FITS file.
"""

function writeClToFITS(f::CFITSIO.FITSFile, Cl::Vector{T}) where {T <: Real}

    idx = Vector{Int64}(range(1, length(Cl)))

    CFITSIO.fits_write_col(f, 1, 1, 1, idx)
    CFITSIO.fits_write_col(f, 2, 1, 1, Cl)

end

function writeClToFITS(fileName, Cl::Vector{T}, overwrite = true) where {T <: Real}
    if overwrite
        f = CFITSIO.fits_clobber_file(fileName)
    else
        f = CFITSIO.fits_create_file(fileName)
    end
    try
        CFITSIO.fits_create_binary_tbl(f, 0, [("Index", "1I", ""),("C_l", "1D", "")], "Cl")
        writeClToFITS(f, Cl)
    finally
    end
    CFITSIO.fits_close_file(f)
end
