"""
    readClFromFITS{T <: Real}(f::CFITSIO.FITSFile, t::Type{T}) -> Vector{T}
    readClFromFITS{T <: Real}(fileName::String, t::Type{T}) -> Vector{T}
Read a set of C_ℓ coefficients from a FITS file.
"""
function readClFromFITS(f::CFITSIO.FITSFile, t::Type{T}; col_num = 2) where {T <: Real}
    numOfRows = CFITSIO.fits_get_num_rows(f)

    Cl = Vector{T}(undef, numOfRows)

    CFITSIO.fits_read_col(f, col_num, 1, 1, Cl)

    return Cl
end

function readClFromFITS(fileName, t::Type{T}; col_num = 2) where {T <: Real}
    f = CFITSIO.fits_open_table(fileName)
    try
        result = readClFromFITS(f, t; col_num = 2)
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

function writeClToFITS(fileName, Cl::Vector{T}; overwrite = true) where {T <: Real}
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


#########################################################################


"""
    dl2cl(dl::AbstractVector{T}, lmin::Integer) where {T <: Real}

Convert a set of ``D_{\\ell}`` to ``C_{\\ell}`` power spectrum, where
``C_{\\ell} = 2\\pi D_{\\ell} / \\ell (\\ell + 1)``. The first components are
set to zero if not present. The monopole component is set to zero in any case to avoid Inf values.

#ARGUMENTS:
- `dl::AbstractVector{T}` : Array of Dℓ components
- `lmin::Integer` : minimum l in the representation of the Dℓ power spectrum

#RETURNS:
- `Vector{T}` : Array of Cℓ power spectrum components
"""

function dl2cl(dl::AbstractVector{T}, lmin::Integer) where {T <: Real}
    lmax = length(dl)+lmin-1
    l_s = Vector{Int}(lmin:lmax)
    cl =  dl .* 2π ./ (l_s .* (l_s .+ 1))
    #fill the missing initial components (monopole and/or dipole) with zeros
    head = zeros(lmin)
    cl = append!(head, cl)
    if lmin == 0
        cl[1] = 0
    end
    return cl
end


"""
    cl2dl(cl::AbstractVector{T}, lmin::Integer) where {T <: Real}

Convert a set of ``C_{\\ell}`` to ``D_{\\ell}`` power spectrum, where
``D_{\\ell} = \\ell (\\ell + 1) C_{\\ell} / 2\\pi``.
The first components are set to zero if not present.

#ARGUMENTS:
- `cl::AbstractVector{T}` : Array of Cℓ components
- `lmin::Integer` : minimum l in the representation of the Dℓ power spectrum

#RETURNS:
- `Vector{T}` : Array of Dℓ power spectrum components
"""

function cl2dl(cl::AbstractVector{T}, lmin::Integer) where {T <: Real}
    lmax = length(cl)+lmin-1
    l_s = Vector{Int}(lmin:lmax)
    head = zeros(lmin)
    #fill the missing initial components (monopole and/or dipole) with zeros
    cl = append!(head, cl)
    dl =  cl ./ 2π .* (l_s .* (l_s .+ 1))
    return dl
end
