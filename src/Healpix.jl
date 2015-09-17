module Healpix

export Resolution, nside2npix, npix2nside, normalizeAngle
export ang2pixNest, ang2pixRing, pix2angNest, pix2angRing
export Order, RingOrder, NestedOrder, Map
export ang2vec, vec2ang, ang2pix, pix2ang
export readMapFromFITS, saveToFITS, conformables, ringWeightPath, readRingWeights
export pixelWindowPath, readPixelWindowT, readPixelWindowP
export Alm, numberOfAlms, almIndexL0, almIndex, readAlmFromFITS

import FITSIO

const NSIDE_MAX = 8192

########################################################################

"""
    nside2npix(nside::Integer) -> Integer

Return the number of pixels for a Healpix map with the specified
`NSIDE` value. If `NSIDE` is not an integer power of two, the function
throws a `DomainError` exception."""

function nside2npix(nside::Integer)
    nsidelog2 = round(Int, log2(nside))
    if 2^nsidelog2 != nside
        throw(DomainError())
    end

    12 * nside * nside
end

########################################################################

"""
    npix2nside(nside::Integer) -> Integer

Given the number of pixels in a Healpix map, return the `NSIDE`
resolution parameter. If the number is invalid, throw a `DomainError`
exception."""

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

const pix2x = [
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  0,  1,  0,  1,  2,  3,  2,  3,  0,  1,  0,  1,  2,  3,  2,  3,
  4,  5,  4,  5,  6,  7,  6,  7,  4,  5,  4,  5,  6,  7,  6,  7,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
  8,  9,  8,  9, 10, 11, 10, 11,  8,  9,  8,  9, 10, 11, 10, 11,
 12, 13, 12, 13, 14, 15, 14, 15, 12, 13, 12, 13, 14, 15, 14, 15,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 16, 17, 16, 17, 18, 19, 18, 19, 16, 17, 16, 17, 18, 19, 18, 19,
 20, 21, 20, 21, 22, 23, 22, 23, 20, 21, 20, 21, 22, 23, 22, 23,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31,
 24, 25, 24, 25, 26, 27, 26, 27, 24, 25, 24, 25, 26, 27, 26, 27,
 28, 29, 28, 29, 30, 31, 30, 31, 28, 29, 28, 29, 30, 31, 30, 31 ]

const pix2y = [
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  0,  0,  1,  1,  0,  0,  1,  1,  2,  2,  3,  3,  2,  2,  3,  3,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  4,  4,  5,  5,  4,  4,  5,  5,  6,  6,  7,  7,  6,  6,  7,  7,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
  8,  8,  9,  9,  8,  8,  9,  9, 10, 10, 11, 11, 10, 10, 11, 11,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
 12, 12, 13, 13, 12, 12, 13, 13, 14, 14, 15, 15, 14, 14, 15, 15,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 16, 16, 17, 17, 16, 16, 17, 17, 18, 18, 19, 19, 18, 18, 19, 19,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 20, 20, 21, 21, 20, 20, 21, 21, 22, 22, 23, 23, 22, 22, 23, 23,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 24, 24, 25, 25, 24, 24, 25, 25, 26, 26, 27, 27, 26, 26, 27, 27,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31,
 28, 28, 29, 29, 28, 28, 29, 29, 30, 30, 31, 31, 30, 30, 31, 31 ]

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

"""`Resolution` objects are needed to perform a number of
pixel-related functions, e.g., convert a direction into a pixel number
and vice versa."""

immutable Resolution
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

"""
    Resolution(nside::Uint32) -> Resolution
    Resolution(nside::Int) -> Resolution

Create a `Resolution` object, given a value for `NSIDE`."""

function Resolution(nside::Uint32)
    if nside > NSIDE_MAX || nside < 1
        throw(DomainError())
    end

    # The expression (nside & (nside - 1)) != 0 is a quick check for
    # detecting if nside is a power of two or not
    order          = ((nside & (nside - 1)) != 0) ? -1 : ilog2(nside)
    pixelsPerFace  = nside * nside
    numOfPixels    = 12 * pixelsPerFace
    ncap           = 2 * nside * (nside - 1)
    fact2          = 4.0 / numOfPixels
    fact1          = 2 * nside * fact2

    result = Resolution(nside,
                        2 * nside,
                        4 * nside,
                        numOfPixels,
                        order,
                        pixelsPerFace,
                        ncap,
                        fact2,
                        fact1)

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
    local jp::Int = floor(Integer, NSIDE_MAX * (0.5 + scaled_phi - z*0.75))
    local jm::Int = floor(Integer, NSIDE_MAX * (0.5 + scaled_phi + z*0.75))

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
    local ntt::Int = floor(Integer, scaled_phi)
    if ntt >= 4
        ntt = 3
    end

    local tp::Float64 = scaled_phi - ntt
    local tmp::Float64 = sqrt(3. * (1. - z_abs)) # in ]0,1]

    local jp::Int = floor(Integer, NSIDE_MAX * tp * tmp)
    local jm::Int = floor(Integer, NSIDE_MAX * (1. - tp) * tmp)

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

"""
    ang2pixNest(resol::Resolution, theta::Float64, phi::Float64) -> Integer

Return the index of the pixel which contains the point with
coordinates (`theta`, the colatitude, and `phi`, the longitude), in
radians, for a Healpix map with pixels in nested order. Note that
pixel indexes are 1-based (this is Julia)!"""

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
    ipf = floor(Integer, ipf / ((NSIDE_MAX / nside) ^ 2))

    # Add 1 to have a 1-based index
    ipf + face_num * nside * nside + 1

end

################################################################################

function calcRingPosForEquator(resol::Resolution,
                               z::Float64,
                               z_abs::Float64,
                               tt::Float64)

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

function calcRingPosForPole(resol::Resolution,
                            z::Float64,
                            z_abs::Float64,
                            tt::Float64)

    const tp = tt - floor(tt)
    const tmp = sqrt(3. * (1. - z_abs))

    jp = floor(Integer, resol.nside * tp * tmp )
    jm = floor(Integer, resol.nside * (1. - tp) * tmp)

    ir = jp + jm + 1
    ip = floor(Integer, tt * ir) + 1
    if ip > 4 * ir
	ip -= 4 * ir
    end

    if z <= 0.
	ipix1 = resol.numOfPixels - 2 * ir * (ir + 1) + ip
    else
        ipix1 = 2 * ir * (ir - 1) + ip
    end

    ipix1
end

################################################################################

"""
    ang2vec(theta::Real, phi::Real) -> (Float64, Float64, Float64)

Given a direction in the sky with colatitude `theta` and longitude
`phi` (in radians), return a tuple containing the `x`, `y`, and `z`
components of the one-length vector pointing to that direction. """

function ang2vec(theta::Real, phi::Real)
    if theta < 0 || theta > π
        throw(DomainError())
    end

    sintheta = sin(theta)
    return (sintheta * cos(phi), sintheta * sin(phi), cos(theta))
end

################################################################################

"""
    vec2ang(x::Real, y::Real, z::Real) -> (Float64, Float64)

Given a vector (with any length) whose Cartesian components are `x`,
`y`, and `z`, return a pair (`theta`, `phi`) containing the colatitude
`theta` and the longitude `phi` (in radians) of the direction in the
sky the vector is pointing at. """

function vec2ang(x::Real, y::Real, z::Real)
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
    ang2pixRing(resol::Resolution, theta::Float64, phi::Float64) -> Integer

Return the index of the pixel which contains the point with
coordinates (`theta`, the colatitude, and `phi`, the longitude), in
radians, for a Healpix map with pixels in ring order. Note that pixel
indexes are 1-based (this is Julia)!"""

function ang2pixRing(resol::Resolution,
                     theta::Float64,
                     phi::Float64)

    const z = cos(theta);
    const z_abs = abs(z);
    const scaled_phi = normalizeAngle(phi) / (0.5 * π) # in [0,4[

    local ipix1 :: Int

    if z_abs <= 2./3.
        ipix1 = calcRingPosForEquator(resol, z, z_abs, scaled_phi)
    else
        ipix1 = calcRingPosForPole(resol, z, z_abs, scaled_phi)
    end

    ipix1
end

################################################################################

"""
    pix2angNest(resol::Resolution, pixel::Int) -> (Float64, Float64)

Given the (1-based) index of a pixel in a Healpix map in nested
order, return a pair containing the (`colatitude`, `longitude`) angles
corresponding to its center, both expressed in radians."""

function pix2angNest(resol::Resolution, pixel::Int)

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
    pix2angRing(resol::Resolution, pixel::Int) -> (Float64, Float64)

Given the (1-based) index of a pixel in a Healpix map in ring
order, return a pair containing the (`colatitude`, `longitude`) angles
corresponding to its center, both expressed in radians."""

function pix2angRing(resol::Resolution, pixel::Int)
    const fact1 = 1.5 * resol.nside
    const fact2 = 3.0 * resol.pixelsPerFace

    # Any reference to equations in this routine refers to Gorski et al. (2005)

    if pixel <= resol.ncap
        # North Polar cap

	const p_h   = pixel / 2 # Defined in Gorsky et al. (2005) before Eq. (2)
	const floor_p_h = floor(p_h)
        # Eq. (2)
	const i = floor(Integer, sqrt(p_h - sqrt(floor_p_h))) + 1 # counted from N. pole
        # Eq. (3)
	const j  = pixel - 2i*(i - 1)

        # Colatitude: Eq. (4); longitude: Eq. (5)
	return (acos(1.0 - i^2 / fact2),
                (Float64(j) - 0.5) * π / (2i))
    elseif pixel <= resol.nsideTimesTwo * (5resol.nside + 1)
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
	const j = Int(4. * i + 1 - (ip - 2. * i * (i - 1)))

	return (acos(-1. + i^2 / fact2),
                (float(j) - 0.5) * π / (2i))
    end
end

################################################################################

"Abstract type representing the ordering of pixels in a Healpix map.
See also `RingOrder` and `NestedOrder`."

abstract Order

"""The `RingOrder` type should be used when creating `Map` types in
order to specify that the pixels in the map are sorted in ring
ordering. (See also `NestedOrder`.)"""

type RingOrder <: Order end

"""The `NestedOrder` type should be used when creating `Map` types in
order to specify that the pixels in the map are sorted in ring
ordering. (See also `RingOrder`.)"""

type NestedOrder <: Order end

"""A Healpix map. The type `T` is used for the value of the pixels in
a map, and it can be anything (even a string!). The type `O` is used
to specify the ordering of the pixels, and it can either be
`RingOrder` or `NestedOrder`."""

type Map{T, O <: Order}
    pixels :: Array{T}
    resolution :: Resolution

    """
    Map{T, O <: Order}(nside::Integer) -> Map{T, O}

Create an empty map with the specified NSIDE."""
    Map(nside::Integer) = new(Array(T, nside2npix(nside)),
                                Resolution(nside))

    "Create a map with the specified NSIDE and initialize the value of
    its pixels using the values in arr."
    function Map{T}(nside :: Integer, arr :: Array{T})
        if nside2npix(nside) != length(arr)
            throw(DomainError())
        end

        new(arr, Resolution(nside))
    end
end

################################################################################

"""
    readMapFromFITS{T <: Number}(f::FITSIO.FITSFILE, column::Integer, t::Type{T})
    readMapFromFITS{T <: Number}(fileName::String, column::Integer, t::Type{T})

Read a Healpix map from the specified (1-base indexed) column in a
FITS file. The values will be read as numbers of type T. If the code
fails, FITSIO will raise an exception. (Refer to the FITSIO library
for more information.)"""

function readMapFromFITS{T <: Number}(f :: FITSIO.FITSFile,
                                      column :: Integer,
                                      t :: Type{T})
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
        result = Map{T, RingOrder}(nside, Array(T, nside2npix(nside)))
    else
        result = Map{T, NestedOrder}(nside, Array(T, nside2npix(nside)))
    end
    FITSIO.fits_read_col(f, column, 1, 1, result.pixels)

    result
end

function readMapFromFITS{T <: Number}(fileName :: String,
                                      column :: Integer,
                                      t :: Type{T})
    f = FITSIO.fits_open_table(fileName)
    result = readMapFromFITS(f, column, t)
    FITSIO.fits_close_file(f)

    result
end

################################################################################

function savePixelsToFITS{T <: Number}(map :: Map{T},
                                       f :: FITSIO.FITSFile,
                                       column :: Integer)

    FITSIO.fits_update_key(f, "NSIDE", map.resolution.nside,
                           "Value of NSIDE")
    FITSIO.fits_write_col(f, column, 1, 1, map.pixels)

end

"""
    saveToFITS{T <: Number, O <: Order}(map::Map{T, O},
                                        f::FITSIO.FITSFile,
                                        column::Integer)
    saveToFITS{T <: Number, O <: Order}(map::Map{T, O},
                                        fileName::String,
                                        typechar="D",
                                        unit="",
                                        extname="MAP")

Save a Healpix map in the specified (1-based index) column in a FITS
file. If the code fails, FITSIO will raise an exception. (Refer to the
FITSIO library for more information.)"""

function saveToFITS{T <: Number}(map :: Map{T, RingOrder},
                                 f :: FITSIO.FITSFile,
                                 column :: Integer)

    FITSIO.fits_update_key(f, "ORDERING", "RING")
    savePixelsToFITS(map, f, column)

end

function saveToFITS{T <: Number}(map :: Map{T, NestedOrder},
                                 f :: FITSIO.FITSFile,
                                 column :: Integer)

    FITSIO.fits_update_key(f, "ORDERING", "NEST")
    savePixelsToFITS(map, f, column)

end

function saveToFITS{T <: Number, O <: Order}(map :: Map{T, O},
                                             fileName :: String,
                                             typechar="D",
                                             unit="",
                                             extname="MAP")

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
shape and ordering are the same."""

conformables{T, S}(map1::Map{T, RingOrder}, map2::Map{S, RingOrder}) =
    map1.resolution.nside == map2.resolution.nside

conformables{T, S}(map1::Map{T, NestedOrder}, map2::Map{S, NestedOrder}) =
    map1.resolution.nside == map2.resolution.nside

conformables{T, S, O1 <: Order, O2 <: Order}(map1::Map{T, O1},
                                             map2::Map{S, O2}) = false

################################################################################

function ringWeightPath(datadir :: UTF8String, nside)
    @sprintf("%s/weight_ring_n%05d.fits", datadir, nside)
end

function readWeightRing(fileName :: UTF8String, nside)
    f = FITSIO.fits_open_table(fileName)
    try
        weights = Array(Float64, 2 * nside)
        FITSIO.fits_read_col(f, 1, 1, 1, weights)
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
        FITSIO.fits_read_col(f, 1, 1, 1, pixwin)
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
        FITSIO.fits_read_col(f, 1, 1, 1, pixwinT)
        FITSIO.fits_read_col(f, 2, 1, 1, pixwinP)
    finally
        FITSIO.fits_close_file(f)
    end

    pixwinT, pixwinP
end

################################################################################

"""
    ang2pix{T, O <: Order}(map::Map{T, O}, theta::Real, phi::Real)

Convert the direction specified by the colatitude `theta` (∈ [0, π])
and the longitude `phi` (∈ [0, 2π]) into the index of the pixel in the
Healpix map `map`.
"""

function ang2pix{T}(map::Map{T, RingOrder}, theta::Real, phi::Real)
    ang2pixRing(map.resolution, Float64(theta), Float64(phi))
end

function ang2pix{T}(map::Map{T, NestedOrder}, theta::Real, phi::Real)
    ang2pixNest(map.resolution, Float64(theta), Float64(phi))
end

################################################################################

"""
    ang2pix{T, O <: Order}(map::Map{T, O}, ipix::Integer) -> (Float64, Float64)

Return the pair (`theta`, `phi`), where `theta` is the colatitude and
`phi` the longitude of the direction of the pixel center with index
`ipix`.
"""

function pix2ang{T}(map::Map{T, RingOrder}, ipix::Integer)
    pix2angRing(map.resolution, ipix)
end

function pix2ang{T}(map::Map{T, NestedOrder}, ipix::Integer)
    pix2angNest(map.resolution, ipix)
end

################################################################################

"""An array of a_ℓm numbers."""

type Alm{T <: Number}
    alm :: Array{T}
    lmax :: Int
    mmax :: Int
    tval :: Int

    Alm(lmax, mmax) = new(Array(T, numberOfAlms(lmax, mmax)),
                          lmax, mmax, 2 * lmax + 1)

    function Alm(lmax :: Integer, mmax :: Integer, arr :: Array{T})
        if numberOfAlms(lmax, mmax) != length(arr)
            throw(DomainError())
        end

        new(arr, lmax, mmax, 2 * lmax + 1)
    end
end

################################################################################

"""
    numberOfAlms(lmax::Integer, mmax::Integer) -> Integer
    numberOfAlms(lmax::Integer) -> Integer

Return the size of the array of complex numbers needed to store the
a_lm coefficients in the range of ℓ and m specified by `lmax` and
`mmax`. If `mmax` is not specified, it is assumed to be equal to
`lmax`. If `lmax` and `mmax` are inconsistent or negative, a
`DomainError` exception is thrown."""

function numberOfAlms(lmax :: Integer, mmax :: Integer)
    if mmax < 0 || lmax < 0 || mmax > lmax
        throw(DomainError())
    end

    div((mmax + 1) * (mmax + 2), 2) + (mmax + 1) * (lmax - mmax)
end

numberOfAlms(lmax :: Integer) = numberOfAlms(lmax, lmax)

shr(x :: Integer, y :: Integer) = x >> y
shr{T <: Integer}(x :: Array{T}, y :: Integer) = [a >> y for a in x]

almIndexL0{T}(alm :: Alm{T}, m) = shr((m .* (alm.tval .- m)), 1) + 1
almIndex{T}(alm :: Alm{T}, l, m) = almIndexL0(alm, m) .+ l

################################################################################

"""
    readAlmFromFITS{T <: Complex}(f::FITSIO.FITSFile, t::Type{T}) -> Alm{T}
    readAlmFromFITS{T <: Complex}(fileName::String, t::Type{T}) -> Alm{T}

Read a set of a_ℓm coefficients from a FITS file. If the code fails,
FITSIO will raise an exception. (Refer to the FITSIO library for more
information.)"""

function readAlmFromFITS{T <: Complex}(f :: FITSIO.FITSFile,
                                       t :: Type{T})
    const numOfRows = FITSIO.fits_get_num_rows(f)

    idx = Array(Int64, numOfRows)
    almReal = Array(Float64, numOfRows)
    almImag = Array(Float64, numOfRows)

    FITSIO.fits_read_col(f, 1, 1, 1, idx)
    FITSIO.fits_read_col(f, 2, 1, 1, almReal)
    FITSIO.fits_read_col(f, 3, 1, 1, almImag)

    l = floor(Int64, sqrt(idx - 1))
    m = idx - l.^2 - l - 1
    if count(x -> x < 0, m) > 0
        throw(DomainError())
    end

    result = Alm{T}(maximum(l), maximum(m))
    i = almIndex(result, l, m)
    result.alm = complex(almReal[i], almImag[i])
end

function readAlmFromFITS{T <: Complex}(fileName :: String,
                                       t :: Type{T})
    f = FITSIO.fits_open_table(fileName)
    try
        result = readAlmFromFITS(f, t)
        return result
    finally
        FITSIO.fits_close_file(f)
    end
end

end
