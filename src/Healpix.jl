module Healpix

export Resolution, nside2npix, npix2nside, normalizeAngle, ang2pixNest
export Ordering, Map, conformables, ringWeightPath, readRingWeights
export pixelWindowPath, readPixelWindowT, readPixelWindowP
export Alm, numberOfAlms

import FITSIO

const NSIDE_MAX = 8192

########################################################################

function nside2npix(nside::Integer)
    nsidelog2 = int(log2(nside))
    if 2^nsidelog2 != nside
        throw(DomainError())
    end

    12 * nside * nside
end

########################################################################

function npix2nside(npix::Integer)
    if npix % 12 != 0
        throw(DomainError())
    end

    square_root::Float64 = sqrt(npix / 12)
    if square_root*square_root != npix/12
        throw(DomainError())
    end

    convert(Int, round(square_root))
end

################################################################################

const x2pix = [
       0,     1,     4,     5,    16,    17,    20,    21,    64,    65,
      68,    69,    80,    81,    84,    85,   256,   257,   260,   261,
     272,   273,   276,   277,   320,   321,   324,   325,   336,   337,
     340,   341,  1024,  1025,  1028,  1029,  1040,  1041,  1044,  1045,
    1088,  1089,  1092,  1093,  1104,  1105,  1108,  1109,  1280,  1281,
    1284,  1285,  1296,  1297,  1300,  1301,  1344,  1345,  1348,  1349,
    1360,  1361,  1364,  1365,  4096,  4097,  4100,  4101,  4112,  4113,
    4116,  4117,  4160,  4161,  4164,  4165,  4176,  4177,  4180,  4181,
    4352,  4353,  4356,  4357,  4368,  4369,  4372,  4373,  4416,  4417,
    4420,  4421,  4432,  4433,  4436,  4437,  5120,  5121,  5124,  5125,
    5136,  5137,  5140,  5141,  5184,  5185,  5188,  5189,  5200,  5201,
    5204,  5205,  5376,  5377,  5380,  5381,  5392,  5393,  5396,  5397,
    5440,  5441,  5444,  5445,  5456,  5457,  5460,  5461 ]

################################################################################

const y2pix = [
        0,     2,     8,    10,    32,    34,    40,    42,   128,   130,
      136,   138,   160,   162,   168,   170,   512,   514,   520,   522,
      544,   546,   552,   554,   640,   642,   648,   650,   672,   674,
      680,   682,  2048,  2050,  2056,  2058,  2080,  2082,  2088,  2090,
     2176,  2178,  2184,  2186,  2208,  2210,  2216,  2218,  2560,  2562,
     2568,  2570,  2592,  2594,  2600,  2602,  2688,  2690,  2696,  2698,
     2720,  2722,  2728,  2730,  8192,  8194,  8200,  8202,  8224,  8226,
     8232,  8234,  8320,  8322,  8328,  8330,  8352,  8354,  8360,  8362,
     8704,  8706,  8712,  8714,  8736,  8738,  8744,  8746,  8832,  8834,
     8840,  8842,  8864,  8866,  8872,  8874, 10240, 10242, 10248, 10250,
    10272, 10274, 10280, 10282, 10368, 10370, 10376, 10378, 10400, 10402,
    10408, 10410, 10752, 10754, 10760, 10762, 10784, 10786, 10792, 10794,
    10880, 10882, 10888, 10890, 10912, 10914, 10920, 10922 ]

################################################################################

function ilog2(argument::Uint32)

    local result = 0
    local shifted_argument = argument

    while shifted_argument > 0x0000FFFF
	result += 16
	shifted_argument >>= 16
    end

    if shifted_argument > 0x000000FF
	result |= 8
	shifted_argument >>= 8
    end

    if shifted_argument > 0x0000000F
	result |= 4
	shifted_argument >>= 4
    end

    if shifted_argument > 0x00000003
	result |= 2
	shifted_argument>>=2
    end

    if shifted_argument > 0x00000001
	result |= 1
    end

    result
end

################################################################################

type Resolution
    nside          :: Uint32
    nsideTimesTwo  :: Uint32
    nsideTimesFour :: Uint32
    numOfPixels    :: Uint32

    order          :: Uint32
    pixelsPerFace  :: Uint32
    ncap           :: Uint32
    fact2          :: Float64
    fact1          :: Float64
end

################################################################################

function Resolution(nside::Uint32)
    if nside > NSIDE_MAX
        error("The maximum allowed value for NSIDE is $NSIDE_MAX")
    end

    result = Resolution(zero(Uint32),
                        zero(Uint32),
                        zero(Uint32),
                        zero(Uint32),
                        zero(Uint32),
                        zero(Uint32),
                        zero(Uint32),
                        zero(Float64),
                        zero(Float64))

    result.nside          =   nside
    result.nsideTimesTwo  = 2 * nside
    result.nsideTimesFour = 4 * nside

    # The expression (nside & (nside - 1)) != 0 is a quick check for
    # detecting if nside is a power of two or not
    result.order          = ((nside & (nside - 1)) != 0) ? -1 : ilog2(nside)
    result.pixelsPerFace  = nside * nside
    result.numOfPixels    = 12 * result.pixelsPerFace
    result.ncap           = 2 * (result.pixelsPerFace - nside)
    result.fact2          = 4.0 / result.numOfPixels
    result.fact1          = 2 * nside * result.fact2

    result
end

################################################################################

function Resolution(nside::Int)
    Resolution(convert(Uint32, nside))
end

################################################################################

# Return the same angle as the argument, but in the range [0, 2π)
function normalizeAngle(x::Float64)
    while x >= 2π
        x -= 2π
    end

    while x < 0
        x += 2π
    end

    x
end

################################################################################

function calcNestPosForEquator(z, z_abs, scaled_phi)
    local jp::Int = ifloor(NSIDE_MAX * (0.5 + scaled_phi - z*0.75))
    local jm::Int = ifloor(NSIDE_MAX * (0.5 + scaled_phi + z*0.75))

    local idfp::Int = div(jp, NSIDE_MAX) # in {0,4}
    local idfm::Int = div(jm, NSIDE_MAX)

    local face_num::Int
    if idfp == idfm
        face_num = (idfp % 4) + 4
    elseif idfp < idfm
        face_num = (idfp % 4)
    else
        face_num = (idfm % 4) + 8
    end

    local ix::Int = mod(jm, NSIDE_MAX)
    local iy::Int = NSIDE_MAX - mod(jp, NSIDE_MAX) - 1

    (ix, iy, face_num)
end

################################################################################

function calcNestPosForPole(z, z_abs, scaled_phi)
    local ntt::Int = ifloor(scaled_phi)
    if ntt >= 4
        ntt = 3
    end

    local tp::Float64 = scaled_phi - ntt
    local tmp::Float64 = sqrt (3. * (1. - z_abs)) # in ]0,1]

    local jp::Int = ifloor(NSIDE_MAX * tp * tmp)
    local jm::Int = ifloor(NSIDE_MAX * (1. - tp) * tmp)

    # Clip jp and jm
    jp = jp < NSIDE_MAX-1 ? jp : NSIDE_MAX-1
    jm = jm < NSIDE_MAX-1 ? jm : NSIDE_MAX-1

    local ix::Int, iy::Int, face_num::Int
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

# Return the index of the pixel which contains the point with
# coordinates (theta, phi). Note that indexes are 1-based (this is
# Julia)!
function ang2pixNest(resol::Resolution,
                     theta::Float64,
                     phi::Float64)

    const nside = resol.nside
    local ix::Int, iy::Int, face_num::Int

    local z::Float64          = cos(theta)
    local z_abs::Float64      = abs(z)
    local scaled_phi::Float64 = normalizeAngle(phi) / (0.5 * π) # in [0,4[

    if z_abs <= 2./3.
        (ix, iy, face_num) = calcNestPosForEquator(z, z_abs, scaled_phi)
    else
        (ix, iy, face_num) = calcNestPosForPole(z, z_abs, scaled_phi)
    end

    (ix_hi, ix_low) = divrem(ix, 128)
    (iy_hi, iy_low) = divrem(iy, 128)

    local ipf = ((x2pix[ix_hi + 1] + y2pix[iy_hi + 1]) * (128 * 128)
                 + (x2pix[ix_low + 1] + y2pix[iy_low + 1]))
    ipf = ifloor(ipf / ((NSIDE_MAX / nside) ^ 2))

    # Add 1 to have a 1-based index
    ipf + face_num * nside * nside + 1

end

################################################################################

immutable Ordering
    num :: Int
    Ordering(n :: Integer) = new(n)
end

const Ring = Ordering(0)
const Nested = Ordering(1)

type Map{T <: Number}
    pixels :: Array{T}
    resolution :: Resolution
    ordering :: Ordering

    Map(nside :: Integer, ordering :: Ordering) = new(Array(T, nside2npix(nside)), 
                                                      Resolution(nside), 
                                                      ordering)

    function Map(nside :: Integer, ordering :: Ordering, arr :: Array{T})
        if nside2npix(nside) != length(arr)
            throw(DomainError())
        end

        new(arr, Resolution(nside), ordering)
    end
end

function Map{T <: Number}(f :: FITSIO.FITSFile, column :: Integer, t :: Type{T})
    value, comment = FITSIO.fits_read_keyword(f, "NSIDE")
    const nside = int(value)

    value, comment = FITSIO.fits_read_keyword(f, "ORDERING")
    const ordering = uppercase(strip(value[2:end-1])) == "RING" ? Ring : Nested

    const repeat, width = FITSIO.fits_get_col_repeat(f, column)
    const nrows = FITSIO.fits_get_num_rows(f)

    if repeat * nrows != nside2npix(nside)
        error("Wrong number of pixels in column $column of FITS file (NSIDE=$nside)")
    end

    result = Map{T}(nside, ordering, Array(T, nside2npix(nside)))
    FITSIO.fits_read_col(f, T, column, 1, 1, result.pixels)

    result
end

function Map{T <: Number}(fileName :: ASCIIString, column :: Integer, t :: Type{T})
    f = FITSIO.fits_open_table(fileName)
    result = Map(f, column, t)
    FITSIO.fits_close_file(f)

    result
end

conformables{T}(map1::Map{T}, map2::Map{T}) = 
    map1.resolution.nside == map2.resolution.nside

################################################################################

function ringWeightPath(datadir :: UTF8String, nside)
    @sprintf("%s/weight_ring_n%05d.fits", datadir, nside)
end

function readWeightRing(fileName :: UTF8String, nside)
    f = FITSIO.fits_open_table(fileName)
    try
        weights = Array(Float64, 2 * nside)
        FITSIO.fits_read_col(f, Float64, 1, 1, 1, weights)
    finally
        FITSIO.fits_close_file(f)
    end

    weights
end

################################################################################

function pixelWindowPath(datadir :: UTF8String, nside)
    @sprintf("%s/pixel_window_n%04d.fits", datadir, nside)
end

function readPixelWindowT(fileName :: UTF8String, nside)
    f = FITSIO.fits_open_table(fileName)
    try
        pixwin = Array(Float64, FITSIO.fits_get_num_rows(f))
        FITSIO.fits_read_col(f, Float64, 1, 1, 1, pixwin)
    finally
        FITSIO.fits_close_file(f)
    end

    pixwin
end

function readPixelWindowP(fileName :: UTF8String, nside)
    f = FITSIO.fits_open_table(fileName)
    try
        pixwinT = Array(Float64, FITSIO.fits_get_num_rows(f))
        pixwinP = Array(Float64, FITSIO.fits_get_num_rows(f))
        FITSIO.fits_read_col(f, Float64, 1, 1, 1, pixwinT)
        FITSIO.fits_read_col(f, Float64, 2, 1, 1, pixwinP)
    finally
        FITSIO.fits_close_file(f)
    end

    pixwinT, pixwinP
end

################################################################################

type Alm{T <: Real}
    alm :: Array{T}
    lmax :: Int
    mmax :: Int

    Alm(lmax :: Integer, mmax :: Integer) = new(Array(T, numberOfAlms(lmax, mmax)), 
                                                lmax, mmax)
end

numberOfAlms(lmax :: Integer, mmax :: Integer) = 
    ((mmax + 1) * (mmax + 2)) / 2 + (mmax + 1) * (lmax - mmax)

end
