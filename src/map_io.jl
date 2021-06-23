################################################################################
# Function to read and write maps from/to files

function readpixelsfromfits(
    f::CFITSIO.FITSFile,
    nside,
    column,
    orderType::Union{Type{RingOrder}, Type{NestedOrder}},
    t::Type{T},
) where {T <: Number}
    
    result = HealpixMap{T, orderType}(
        Array{T}(undef, nside2npix(nside)),
    )

    CFITSIO.fits_read_col(f, column, 1, 1, result.pixels)

    result
end

function readpixelsfromfits(
    f::CFITSIO.FITSFile,
    nside,
    column,
    orderType::Union{Type{RingOrder}, Type{NestedOrder}},
    t::Type{Union{Nothing, T}}
) where {T <: Number}
    
    result = HealpixMap{Union{Nothing, T}, orderType}(nside)
    tmparray = Vector{T}(undef, nside2npix(nside))

    CFITSIO.fits_read_col(f, column, 1, 1, tmparray)

    for i in 1:length(tmparray)
        result.pixels[i] = (tmparray[i] ≈ UNSEEN) ? nothing : tmparray[i]
    end
    
    result
end

"""
    readMapFromFITS{T}(f::CFITSIO.FITSFILE, column, t::Type{T})
    readMapFromFITS{T}(fileName::String, column, t::Type{T})

Read a Healpix map from the specified (1-base indexed) column in a
FITS file. The values will be read as numbers of type T. If the code
fails, CFITSIO will raise an exception. (Refer to the CFITSIO library
for more information.)
"""
function readMapFromFITS(f::CFITSIO.FITSFile, column, t::Type{T}) where {T}
    value, comment = CFITSIO.fits_read_keyword(f, "NSIDE")
    nside = parse(Int, value)

    value, comment = CFITSIO.fits_read_keyword(f, "ORDERING")
    ringOrdering = uppercase(strip(value[2:end - 1])) == "RING"

    repeat = (CFITSIO.fits_get_coltype(f, column))[2]
    nrows = CFITSIO.fits_get_num_rows(f)

    if repeat * nrows != nside2npix(nside)
        error("Wrong number of pixels in column $column of FITS file (NSIDE=$nside)")
    end

    orderType = ringOrdering ? RingOrder : NestedOrder
    readpixelsfromfits(f, nside, column, orderType, t)
end

function readMapFromFITS(fileName::AbstractString, column, t::Type{T}) where {T}
    f = CFITSIO.fits_open_table(fileName)
    result = readMapFromFITS(f, column, t)
    CFITSIO.fits_close_file(f)

    result
end

"""
    readPolarizedMapFromFITS{T}(fileName::AbstractString, column, t::Type{T})

Read a polarized map (I/Q/U) from a FITS file and return a
[`PolarizedHealpixMap`](@ref) object.

The parameter `column` can be either a number or a 3-element tuple. In the
first case, three consecutive columns will be read starting from `column`
(1-based index), otherwise the three column indices will be used.

"""
function readPolarizedMapFromFITS(
    fileName::AbstractString,
    column,
    t::Type{T},
) where {T}

    if length(column) == 1
        column_i = 1
        column_q = 2
        column_u = 3
    elseif length(column) == 3
        column_i, column_q, column_u = column
    else
        throw(DomainError(
            column,
            "The column must be either a number or three numbers (I/Q/U)",
        ))
    end

    f = CFITSIO.fits_open_table(fileName)
    i, q, u = (readMapFromFITS(f, colidx, t) for colidx in (column_i, column_q, column_u))
    CFITSIO.fits_close_file(f)

    PolarizedHealpixMap(i, q, u)
end

################################################################################

function writepixelstofits(f::CFITSIO.FITSFile, column, arr::Vector{T}) where {T <: Number}
    CFITSIO.fits_write_col(f, column, 1, 1, arr)
end

function writepixelstofits(f::CFITSIO.FITSFile, column, arr::Vector{Union{T, Nothing}}) where {T <: Number}
    converted = T[isnothing(x) ? T(UNSEEN) : x for x in arr]
    CFITSIO.fits_write_col(f, column, 1, 1, converted)
end

"""
    savePixelsToFITS(map::HealpixMap{T}, f::CFITSIO.FITSFile, column) where {T <: Number}

Save the pixels of `map` into the column with index/name `column` in the FITS
file, which must have been already opened.
"""
function savePixelsToFITS(
    map::HealpixMap{T},
    f::CFITSIO.FITSFile,
    column;
    write_keywords=true,
) where {T}

    if write_keywords
        CFITSIO.fits_update_key(f, "PIXTYPE", "HEALPIX", "HEALPIX pixelisation")
        CFITSIO.fits_update_key(f, "NSIDE", map.resolution.nside, "Value of NSIDE")
        CFITSIO.fits_update_key(f, "FIRSTPIX", 1, "First pixel (1 based)")
        CFITSIO.fits_update_key(
            f,
            "LASTPIX",
            map.resolution.numOfPixels,
            "Last pixel (1 based)",
        )
        CFITSIO.fits_update_key(f, "INDXSCHM", "IMPLICIT", "Indexing: IMPLICIT or EXPLICIT")
    end

    writepixelstofits(f, column, map.pixels)
end

"""
    saveToFITS{T, O <: Order}(map::HealpixMap{T, O},
                              f::CFITSIO.FITSFile,
                              column)
    saveToFITS{T, O <: Order}(map::HealpixMap{T, O},
                              fileName::String,
                              typechar="D",
                              unit="",
                              extname="MAP")

Save a Healpix map in the specified (1-based index) column in a FITS
file. If the code fails, CFITSIO will raise an exception. (Refer to the
CFITSIO library for more information.)
"""
function saveToFITS(
    map::HealpixMap{T,RingOrder},
    f::CFITSIO.FITSFile,
    column;
    write_keywords=true,
) where {T}

    CFITSIO.fits_update_key(f, "ORDERING", "RING")
    savePixelsToFITS(map, f, column, write_keywords=write_keywords)

end

function saveToFITS(
    map::HealpixMap{T,NestedOrder},
    f::CFITSIO.FITSFile,
    column;
    write_keywords=true,
) where {T}

    CFITSIO.fits_update_key(f, "ORDERING", "NEST")
    savePixelsToFITS(map, f, column, write_keywords=write_keywords)

end

function saveToFITS(
    map::HealpixMap{T,O},
    fileName::AbstractString;
    typechar="D",
    unit="",
    extname="MAP",
) where {T, O <: Order}

    f = CFITSIO.fits_create_file(fileName)
    try
        CFITSIO.fits_create_binary_tbl(f, 0, [("PIXELS", "1$typechar", unit)], extname)
        saveToFITS(map, f, 1)
    finally
        CFITSIO.fits_close_file(f)
    end

end

function saveToFITS(
    map::PolarizedHealpixMap{T,O},
    fileName::AbstractString;
    typechar="D",
    unit="",
    extname="MAP",
) where {T, O <: Order}

    f = CFITSIO.fits_create_file(fileName)
    try
        CFITSIO.fits_create_binary_tbl(
            f,
            0,
            [
                ("I", "1$typechar", unit),
                ("Q", "1$typechar", unit),
                ("U", "1$typechar", unit),
            ],
            extname,
        )
        saveToFITS(map.i, f, 1)
        saveToFITS(map.q, f, 2, write_keywords=false)
        saveToFITS(map.u, f, 3, write_keywords=false)
    finally
        CFITSIO.fits_close_file(f)
    end

end

@doc raw"""
    saveToFITS(map::HealpixMap{T, O}, filename::AbstractString, typechar="D", unit="", extname="MAP") where {T <: Number, O <: Order}
    saveToFITS(map::PolarizedHealpixMap{T, O}, filename::AbstractString, typechar="D", unit="", extname="MAP") where {T <: Number, O <: Order}

Save a map into a FITS file. The name of the file is specified in
`filename`; if it begins with `!`, existing files will be overwritten
without warning. The parameter `typechar` specifies the data type to
be used in the FITS file: the default (`D`) will save 64-bit
floating-point values. See the CCFITSIO documentation for other
values. The keyword `unit` specifies the measure unit used for the
pixels in the map. The keyword `extname` specifies the name of the HDU
where the map pixels will be written.

"""
saveToFITS
