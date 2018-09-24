module Healpix
###
export nsideok, nside2pixarea, nside2resol
export Resolution, nside2npix, npix2nside
export ang2pixNest, ang2pixRing, pix2angNest, pix2angRing
export vec2pixNest, vec2pixRing, pix2vecNest, pix2vecRing
export Order, RingOrder, NestedOrder, Map
export ang2vec, vec2ang, ang2pix, pix2ang
export readMapFromFITS, savePixelsToFITS, saveToFITS, conformables, ringWeightPath, readWeightRing
export pixelWindowPath, readPixelWindowT, readPixelWindowP
export Alm, numberOfAlms, almIndexL0, almIndex, readAlmFromFITS
export getringinfo!, getringinfo, getinterpolRing
export pix2xyfRing, xyf2pixRing, pix2xyfNest, xyf2pixNest
export ring2nest, nest2ring

import FITSIO

const NSIDE_MAX = 8192

########################################################################

"""
    nsideok(nside::Integer) -> Bool

Check whether `nside` is a valid `NSIDE` parameter.
"""
nsideok(nside::Integer) = (nside > 0) && ((nside) & (nside - 1) == 0)

########################################################################

"""
    nside2npix(nside::Integer) -> Integer

Return the number of pixels for a Healpix map with the specified
`NSIDE` value. If `NSIDE` is not an integer power of two, the function
throws a `DomainError` exception.
"""
function nside2npix(nside::Integer)
    nsidelog2 = round(Int, log2(nside))
    (2^nsidelog2 == nside) || throw(DomainError())

    12(nside^2)
end

########################################################################

"""
    npix2nside(npix::Integer) -> Integer

Given the number of pixels in a Healpix map, return the `NSIDE`
resolution parameter. If the number is invalid, throw a `DomainError`
exception.
"""
function npix2nside(npix::Integer)
    (npix % 12 == 0) || throw(DomainError())

    square_root = sqrt(npix / 12)
    (square_root^2 == npix / 12) || throw(DomainError())

    convert(Int, round(square_root))
end

################################################################################

"""
    nside2pixarea(nside::Integer) -> Real

Return the solid angle of a pixel in a map with the specified `NSIDE` parameter.
The result is expressed in steradians.
"""
nside2pixarea(nside::Integer) = 4π / nside2npix(nside)

################################################################################

"""
    nside2resol(nside::Integer) -> Real

Return the approximate resolution of a map with the specified `NSIDE`. The
resolution is expressed in radians, and it is the square root of the pixel
size.
"""
nside2resol(nside::Integer) = sqrt(nside2pixarea(nside))

################################################################################

include("math.jl")
include("datatables.jl")
include("resolution.jl")

################################################################################

function calcNestPosForEquator(z, z_abs, scaled_phi)
    jp = floor(Integer, NSIDE_MAX * (0.5 + scaled_phi - z * 0.75))
    jm = floor(Integer, NSIDE_MAX * (0.5 + scaled_phi + z * 0.75))

    idfp = div(jp, NSIDE_MAX) # in {0,4}
    idfm = div(jm, NSIDE_MAX)

    local face_num
    if idfp == idfm
        face_num = (idfp % 4) + 4
    elseif idfp < idfm
        face_num = (idfp % 4)
    else
        face_num = (idfm % 4) + 8
    end

    ix = mod(jm, NSIDE_MAX)
    iy = NSIDE_MAX - mod(jp, NSIDE_MAX) - 1

    (ix, iy, face_num)
end

################################################################################

function calcNestPosForPole(z, z_abs, scaled_phi)
    ntt = floor(Integer, scaled_phi)
    if ntt >= 4
        ntt = 3
    end

    tp = scaled_phi - ntt
    tmp = sqrt(3 * (1 - z_abs)) # in ]0,1]

    jp = floor(Integer, NSIDE_MAX * tp * tmp)
    jm = floor(Integer, NSIDE_MAX * (1 - tp) * tmp)

    # Clip jp and jm
    jp = jp < NSIDE_MAX - 1 ? jp : NSIDE_MAX - 1
    jm = jm < NSIDE_MAX - 1 ? jm : NSIDE_MAX - 1

    local ix, iy, face_num
    if z >= 0
        face_num = ntt # in {0,3}
        ix = NSIDE_MAX - jm - 1
        iy = NSIDE_MAX - jp - 1
    else
        face_num = ntt + 8 # in {8,11} */
        ix =  jp
        iy =  jm
    end

    (ix, iy, face_num)
end

################################################################################

"""
    ang2pixNest(resol::Resolution, theta, phi) -> Integer

Return the index of the pixel which contains the point with
coordinates (`theta`, the colatitude, and `phi`, the longitude), in
radians, for a Healpix map with pixels in nested order. Note that
pixel indexes are 1-based (this is Julia)!
"""
function ang2pixNest(resol::Resolution, theta, phi)

    const nside = resol.nside
    local ix, iy, face_num

    z          = cos(theta)
    z_abs      = abs(z)
    scaled_phi = mod2pi(phi) / (π / 2) # in [0,4[

    if z_abs ≤ 2//3
        (ix, iy, face_num) = calcNestPosForEquator(z, z_abs, scaled_phi)
    else
        (ix, iy, face_num) = calcNestPosForPole(z, z_abs, scaled_phi)
    end

    (ix_hi, ix_low) = divrem(ix, 128)
    (iy_hi, iy_low) = divrem(iy, 128)

    ipf = ((x2pix[ix_hi + 1] + y2pix[iy_hi + 1]) * (128^2)
           + (x2pix[ix_low + 1] + y2pix[iy_low + 1]))
    ipf = floor(Integer, ipf / ((NSIDE_MAX / nside) ^ 2))

    # Add 1 to have a 1-based index
    ipf + face_num * nside * nside + 1

end

################################################################################

function calcRingPosForEquator(resol::Resolution, z, z_abs, tt)

    const jp = floor(Integer, resol.nside * (0.5 + tt - z * 0.75))
    const jm = floor(Integer, resol.nside * (0.5 + tt + z * 0.75))

    const ir = resol.nside + 1 + jp - jm
    const kshift = (mod(ir, 2) == 0) ? 1 : 0

    const nl4 = resol.nsideTimesFour;

    local ip = div(jp + jm - resol.nside + kshift + 1, 2) + 1
    if ip > nl4
	    ip = ip - nl4
    end

    resol.ncap + nl4 * (ir - 1) + ip
end

################################################################################

function calcRingPosForPole(resol::Resolution, z, z_abs, tt)

    const tp = tt - floor(tt)
    const tmp = sqrt(3. * (1. - z_abs))

    jp = floor(Integer, resol.nside * tp * tmp )
    jm = floor(Integer, resol.nside * (1 - tp) * tmp)

    ir = jp + jm + 1
    ip = floor(Integer, tt * ir) + 1
    if ip > 4ir
	    ip -= 4ir
    end

    if z ≤ 0
	    resol.numOfPixels - 2ir * (ir + 1) + ip
    else
        2ir * (ir - 1) + ip
    end
end

################################################################################

"""
    ang2vec(theta, phi) -> Array{Float64}

Given a direction in the sky with colatitude `theta` and longitude
`phi` (in radians), return an array of 3 elements containing the
`x`, `y`, and `z` components of the one-length vector pointing to
that direction.
"""
function ang2vec(theta, phi)
    (0 ≤ theta ≤ π) || throw(DomainError())

    sintheta = sin(theta)
    return [sintheta * cos(phi), sintheta * sin(phi), cos(theta)]
end

################################################################################

"""
    vec2ang(x, y, z) -> (Number, Number)

Given a vector (not necessarily normalized) whose Cartesian components are `x`,
`y`, and `z`, return a pair (`theta`, `phi`) containing the colatitude `theta`
and the longitude `phi` (in radians) of the direction in the sky the vector is
pointing at.
"""
function vec2ang(x, y, z)
    norm = sqrt(x^2 + y^2 + z^2)
    theta = acos(z / norm)
    phi = atan2(y, x)
    if phi < 0
        phi += 2π
    end

    return (theta, phi)
end

################################################################################

"""
    ang2pixRing(resol::Resolution, theta, phi) -> Integer

Return the index of the pixel which contains the point with
coordinates (`theta`, the colatitude, and `phi`, the longitude), in
radians, for a Healpix map with pixels in ring order. Note that pixel
indexes are 1-based (this is Julia)!
"""
function ang2pixRing(resol::Resolution, theta, phi)

    const z = cos(theta)
    const z_abs = abs(z)
    const scaled_phi = mod2pi(phi) / (π / 2) # in [0,4[

    if z_abs ≤ 2//3
        calcRingPosForEquator(resol, z, z_abs, scaled_phi)
    else
        calcRingPosForPole(resol, z, z_abs, scaled_phi)
    end
end

################################################################################

"""
    pix2angNest(resol::Resolution, pixel) -> (Float64, Float64)

Given the (1-based) index of a pixel in a Healpix map in nested
order, return a pair containing the (`colatitude`, `longitude`) angles
corresponding to its center, both expressed in radians.
"""
function pix2angNest(resol::Resolution, pixel)

    const jrll = [ 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4 ]
    const jpll = [ 1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7 ]

    const floatNside = Float64(resol.nside)
    const fact1 = 1. / (3.0 * floatNside^2)
    const fact2 = 2. / (3.0 * floatNside)

    # face number in {0,11} and pixel number within the face
    const faceNum, ipf = divrem(pixel - 1, resol.pixelsPerFace)

    const ip_trunc, ip_low = divrem(ipf, 1024)
    const ip_hi, ip_med = divrem(ip_trunc, 1024)

    const ix = 1024 * pix2x[ip_hi+1] + 32 * pix2x[ip_med+1] + pix2x[ip_low+1]
    const iy = 1024 * pix2y[ip_hi+1] + 32 * pix2y[ip_med+1] + pix2y[ip_low+1]

    # Transforms this in (horizontal, vertical) coordinates
    const jrt = ix + iy # 'vertical' in {0,2*(nside-1)}
    const jpt = ix - iy # 'horizontal' in {-nside+1,nside-1}

    const jr = jrll[faceNum+1] * resol.nside - jrt - 1
    local nr = resol.nside # Equatorial region (the most frequent)
    const z = (2resol.nside - jr) * fact2
    local kshift = Int(mod(jr - resol.nside, 2))
    if jr < resol.nside
        # North polar cap
	nr = jr
	z = 1. - nr^2 * fact1
	kshift = 0
    elseif jr > 3 * resol.nside
        # South polar cap
	nr = resol.nsideTimesFour - jr
	z = - 1. + nr^2 * fact1
	kshift = 0
    end

    local jp = div(jpll[faceNum+1] * nr + jpt + 1 + kshift, 2)
    if jp > resol.nsideTimesFour
        jp = jp - resol.nsideTimesFour
    end
    if jp < 1
        jp = jp + resol.nsideTimesFour
    end

    (acos(z), (jp - (kshift + 1) * 0.5) * (π / (2nr)))
end

################################################################################

"""
    pix2angRing(resol::Resolution, pixel) -> (Float64, Float64)

Given the (1-based) index of a pixel in a Healpix map in ring
order, return a pair containing the (`colatitude`, `longitude`) angles
corresponding to its center, both expressed in radians.
"""
function pix2angRing(resol::Resolution, pixel)
    const fact1 = 1.5 * resol.nside
    const fact2 = 3.0 * resol.pixelsPerFace

    # Any reference to equations in this routine refers to Gorski et al. (2005)

    if pixel ≤ resol.ncap
        # North Polar cap

        const p_h   = pixel / 2 # Defined in Gorsky et al. (2005) before Eq. (2)
        const floor_p_h = floor(p_h)
            # Eq. (2)
        const i = floor(Integer, sqrt(p_h - sqrt(floor_p_h))) + 1 # counted from N. pole
            # Eq. (3)
        const j  = pixel - 2i*(i - 1)

            # Colatitude: Eq. (4); longitude: Eq. (5)
        return (acos(1 - i^2 / fact2),
                    (Float64(j) - 0.5) * π / (2i))
    elseif pixel ≤ resol.nsideTimesTwo * (5resol.nside + 1)
        # Equatorial belt

        # Zero-based index for pixels in the equatorial region
        const ip    = pixel - resol.ncap - 1
        # Eq. (6) - ring counts from the North pole; resol.nside is
        # the number of pixels in the North Polar ring
        const i = floor(Integer, ip / resol.nsideTimesFour) + resol.nside
        # Eq. (7) - zero-based index of the pixel within this ring
        const j  = Int(mod(ip, resol.nsideTimesFour)) + 1

        # Eq. (9) - this equals 1 if i + resol.nside is odd, 1/2
        # otherwise. It is used to convert j into a longitude (since
        # pixel centers in odd rings are shifted with respect to
        # centers in even rings)
        const s_half = 0.5 * (1 + mod(Float64(i + resol.nside), 2))

        # Colatitude: Eq. (8) in disguise, latitude: Eq. (9)
        return (acos((resol.nsideTimesTwo - i) / fact1),
                    (Float64(j) - s_half) * π / (2resol.nside))
    else
        # South Polar cap

        # The pixels in this cap are handled like the ones in the
        # North Polar cap, except that we must flip the value of "ip".

        const ip = resol.numOfPixels - pixel + 1
        const p_h = ip / 2
        const floor_p_h = floor(p_h)
        const i = floor(Integer, sqrt(p_h - sqrt(floor_p_h))) + 1 # counted from S. pole
        const j = Int(4 * i + 1 - (ip - 2i * (i - 1)))

        return (acos(-1 + i^2 / fact2),
                    (float(j) - 0.5) * π / (2i))
    end
end

################################################################################

vec2pixNest(res::Resolution, x, y, z) = ang2pixNest(res, vec2ang(x, y, z)...)
vec2pixRing(res::Resolution, x, y, z) = ang2pixRing(res, vec2ang(x, y, z)...)
pix2vecNest(res::Resolution, pixel) = ang2vec(pix2angNest(res, pixel)...)
pix2vecRing(res::Resolution, pixel) = ang2vec(pix2angRing(res, pixel)...)

################################################################################

include("interp.jl")
include("xyf.jl")

################################################################################

"""Abstract type representing the ordering of pixels in a Healpix map.
See also `RingOrder` and `NestedOrder`.
"""
abstract type Order end

"""The `RingOrder` type should be used when creating `Map` types in
order to specify that the pixels in the map are sorted in ring
ordering. (See also `NestedOrder`.)
"""
abstract type RingOrder <: Order end

"""The `NestedOrder` type should be used when creating `Map` types in
order to specify that the pixels in the map are sorted in ring
ordering. (See also `RingOrder`.)
"""
abstract type NestedOrder <: Order end

"""A Healpix map. The type `T` is used for the value of the pixels in
a map, and it can be anything (even a string!). The type `O` is used
to specify the ordering of the pixels, and it can either be
`RingOrder` or `NestedOrder`.
"""
mutable struct Map{T, O <: Order}
    pixels::Array{T}
    resolution::Resolution

    """
        Map{T, O <: Order}(nside) -> Map{T, O}

    Create an empty map with the specified NSIDE.
    """
    Map{T, O}(nside::Number) where {T, O <: Order} = new(zeros(T, nside2npix(nside)),
                                                 Resolution(nside))

    """
    Create a map with the specified array of pixels.
    """
    function Map{T, O}(arr::Array{T}) where {T, O <: Order}
        nside = npix2nside(length(arr))    
        new(arr, Resolution(nside))
    end
end

import Base: +, -, *, /

+(a::Map{T,O}, b::Map{T,O}) where {T <: Number, O} = Map{T, O}(a.pixels .+ b.pixels)
-(a::Map{T,O}, b::Map{T,O}) where {T <: Number, O} = Map{T, O}(a.pixels .- b.pixels)
*(a::Map{T,O}, b::Map{T,O}) where {T <: Number, O} = Map{T, O}(a.pixels .* b.pixels)
/(a::Map{T,O}, b::Map{T,O}) where {T <: Number, O} = Map{T, O}(a.pixels ./ b.pixels)

+(a::Map{T,O}, b::Number) where {T <: Number, O} = Map{T, O}(a.pixels .+ b)
-(a::Map{T,O}, b::Number) where {T <: Number, O} = a + (-b)
*(a::Map{T,O}, b::Number) where {T <: Number, O} = Map{T, O}(a.pixels .* b)
/(a::Map{T,O}, b::Number) where {T <: Number, O} = Map{T, O}(a.pixels ./ b)

+(a::Number, b::Map{T,O}) where {T <: Number, O} = b + a
-(a::Number, b::Map{T,O}) where {T <: Number, O} = b + (-a)
*(a::Number, b::Map{T,O}) where {T <: Number, O} = b * a
/(a::Number, b::Map{T,O}) where {T <: Number, O} = Map{T, O}(a ./ b.pixels)

################################################################################

"""
    readMapFromFITS{T <: Number}(f::FITSIO.FITSFILE, column, t::Type{T})
    readMapFromFITS{T <: Number}(fileName::String, column, t::Type{T})

Read a Healpix map from the specified (1-base indexed) column in a
FITS file. The values will be read as numbers of type T. If the code
fails, FITSIO will raise an exception. (Refer to the FITSIO library
for more information.)
"""
function readMapFromFITS(f::FITSIO.FITSFile,
                         column,
                         t::Type{T}) where {T <: Number}
    value, comment = FITSIO.fits_read_keyword(f, "NSIDE")
    const nside = parse(Int, value)

    value, comment = FITSIO.fits_read_keyword(f, "ORDERING")
    const ringOrdering = uppercase(strip(value[2:end-1])) == "RING"

    const repeat = (FITSIO.fits_get_coltype(f, column))[2]
    const nrows = FITSIO.fits_get_num_rows(f)

    if repeat * nrows != nside2npix(nside)
        error("Wrong number of pixels in column $column of FITS file (NSIDE=$nside)")
    end

    if ringOrdering
        result = Map{T, RingOrder}(Array{T}(nside2npix(nside)))
    else
        result = Map{T, NestedOrder}(Array{T}(nside2npix(nside)))
    end
    FITSIO.fits_read_col(f, column, 1, 1, result.pixels)

    result
end

function readMapFromFITS(fileName::AbstractString,
                         column,
                         t::Type{T}) where {T <: Number}
    f = FITSIO.fits_open_table(fileName)
    result = readMapFromFITS(f, column, t)
    FITSIO.fits_close_file(f)

    result
end

################################################################################

"""
    savePixelsToFITS(map::Map{T}, f::FITSIO.FITSFile, column) where {T <: Number}

Save the pixels of `map` into the column with index/name `column` in the FITS
file, which must have been already opened.
"""
function savePixelsToFITS(map::Map{T},
                          f::FITSIO.FITSFile,
                          column) where {T <: Number}

    FITSIO.fits_update_key(f, "PIXTYPE", "HEALPIX",
                           "HEALPIX pixelisation")
    FITSIO.fits_update_key(f, "NSIDE", map.resolution.nside,
                           "Value of NSIDE")
    FITSIO.fits_update_key(f, "FIRSTPIX", 1,
                           "First pixel (1 based)")
    FITSIO.fits_update_key(f, "LASTPIX", map.resolution.numOfPixels,
                           "Last pixel (1 based)")
    FITSIO.fits_update_key(f, "INDXSCHM", "IMPLICIT",
                           "Indexing: IMPLICIT or EXPLICIT")
    FITSIO.fits_write_col(f, column, 1, 1, map.pixels)

end

"""
    saveToFITS{T <: Number, O <: Order}(map::Map{T, O},
                                        f::FITSIO.FITSFile,
                                        column)
    saveToFITS{T <: Number, O <: Order}(map::Map{T, O},
                                        fileName::String,
                                        typechar="D",
                                        unit="",
                                        extname="MAP")

Save a Healpix map in the specified (1-based index) column in a FITS
file. If the code fails, FITSIO will raise an exception. (Refer to the
FITSIO library for more information.)
"""
function saveToFITS(map::Map{T, RingOrder},
                    f::FITSIO.FITSFile,
                    column) where {T <: Number}

    FITSIO.fits_update_key(f, "ORDERING", "RING")
    savePixelsToFITS(map, f, column)

end

function saveToFITS(map::Map{T, NestedOrder},
                    f::FITSIO.FITSFile,
                    column) where {T <: Number}

    FITSIO.fits_update_key(f, "ORDERING", "NEST")
    savePixelsToFITS(map, f, column)

end

"""
    saveToFITS(map::Map{T, O}, filename::AbstractString, typechar="D", unit="", extname="MAP") where {T <: Number, O <: Order}

Save a map into a FITS file. The name of the file is specified in `filename`; if it begins with `!`,
existing files will be overwritten without warning. The parameter `typechar` specifies the data type
to be used in the FITS file: the default (`D`) will save 64-bit floating-point values. See the
CFITSIO documentation for other values. The keyword `unit` specifies the measure unit used for the
pixels in the map. The keyword `extname` specifies the name of the HDU where the map pixels will
be written.
"""
function saveToFITS(map::Map{T, O},
                    fileName::AbstractString;
                    typechar="D",
                    unit="",
                    extname="MAP") where {T <: Number, O <: Order}

    f = FITSIO.fits_create_file(fileName)
    try
        FITSIO.fits_create_binary_tbl(f, 0, [("PIXELS", "1$typechar", unit)], extname)
        saveToFITS(map, f, 1)
    finally
        FITSIO.fits_close_file(f)
    end

end

################################################################################

"""
    conformables{T, S, O1 <: Order, O2 <: Order}(map1::Map{T, O1},
                                                 map2::Map{S, O2}) -> Bool

Determine if two Healpix maps are "conformables", i.e., if their
shape and ordering are the same.
"""
conformables(map1::Map{T, RingOrder}, map2::Map{S, RingOrder}) where {T, S} =
    map1.resolution.nside == map2.resolution.nside

conformables(map1::Map{T, NestedOrder}, map2::Map{S, NestedOrder}) where {T, S} =
    map1.resolution.nside == map2.resolution.nside

conformables(map1::Map{T, O1},
             map2::Map{S, O2}) where {T, S, O1 <: Order, O2 <: Order} = false

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

################################################################################

"""
    ang2pix{T, O <: Order}(map::Map{T, O}, theta::Real, phi::Real)

Convert the direction specified by the colatitude `theta` (∈ [0, π])
and the longitude `phi` (∈ [0, 2π]) into the index of the pixel in the
Healpix map `map`.
"""
function ang2pix(map::Map{T, RingOrder}, theta, phi) where {T}
    ang2pixRing(map.resolution, Float64(theta), Float64(phi))
end

function ang2pix(map::Map{T, NestedOrder}, theta, phi) where {T}
    ang2pixNest(map.resolution, Float64(theta), Float64(phi))
end

################################################################################

"""
    pix2ang{T, O <: Order}(map::Map{T, O}, ipix) -> (Float64, Float64)

Return the pair (`theta`, `phi`), where `theta` is the colatitude and
`phi` the longitude of the direction of the pixel center with index
`ipix`.
"""
function pix2ang(map::Map{T, RingOrder}, ipix) where {T}
    pix2angRing(map.resolution, ipix)
end

function pix2ang(map::Map{T, NestedOrder}, ipix) where {T}
    pix2angNest(map.resolution, ipix)
end

################################################################################

include("projections.jl")
include("alm.jl")
include("mapmaking.jl")

end
