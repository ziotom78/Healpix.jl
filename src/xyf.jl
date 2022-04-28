const JRLL = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
const JPLL = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7]

const UTAB = [
    0, 1, 4, 5, 16, 17,
    20, 21, 64, 65, 68, 69,
    80, 81, 84, 85, 256, 257,
    260, 261, 272, 273, 276, 277,
    320, 321, 324, 325, 336, 337,
    340, 341, 1024, 1025, 1028, 1029,
    1040, 1041, 1044, 1045, 1088, 1089,
    1092, 1093, 1104, 1105, 1108, 1109,
    1280, 1281, 1284, 1285, 1296, 1297,
    1300, 1301, 1344, 1345, 1348, 1349,
    1360, 1361, 1364, 1365, 4096, 4097,
    4100, 4101, 4112, 4113, 4116, 4117,
    4160, 4161, 4164, 4165, 4176, 4177,
    4180, 4181, 4352, 4353, 4356, 4357,
    4368, 4369, 4372, 4373, 4416, 4417,
    4420, 4421, 4432, 4433, 4436, 4437,
    5120, 5121, 5124, 5125, 5136, 5137,
    5140, 5141, 5184, 5185, 5188, 5189,
    5200, 5201, 5204, 5205, 5376, 5377,
    5380, 5381, 5392, 5393, 5396, 5397,
    5440, 5441, 5444, 5445, 5456, 5457,
    5460, 5461, 16384, 16385, 16388, 16389,
    16400, 16401, 16404, 16405, 16448, 16449,
    16452, 16453, 16464, 16465, 16468, 16469,
    16640, 16641, 16644, 16645, 16656, 16657,
    16660, 16661, 16704, 16705, 16708, 16709,
    16720, 16721, 16724, 16725, 17408, 17409,
    17412, 17413, 17424, 17425, 17428, 17429,
    17472, 17473, 17476, 17477, 17488, 17489,
    17492, 17493, 17664, 17665, 17668, 17669,
    17680, 17681, 17684, 17685, 17728, 17729,
    17732, 17733, 17744, 17745, 17748, 17749,
    20480, 20481, 20484, 20485, 20496, 20497,
    20500, 20501, 20544, 20545, 20548, 20549,
    20560, 20561, 20564, 20565, 20736, 20737,
    20740, 20741, 20752, 20753, 20756, 20757,
    20800, 20801, 20804, 20805, 20816, 20817,
    20820, 20821, 21504, 21505, 21508, 21509,
    21520, 21521, 21524, 21525, 21568, 21569,
    21572, 21573, 21584, 21585, 21588, 21589,
    21760, 21761, 21764, 21765, 21776, 21777,
    21780, 21781, 21824, 21825, 21828, 21829,
    21840, 21841, 21844, 21845,
]

const CTAB = [
    0, 1, 256, 257, 2, 3, 258, 259,
    512, 513, 768, 769, 514, 515, 770, 771,
    4, 5, 260, 261, 6, 7, 262, 263,
    516, 517, 772, 773, 518, 519, 774, 775,
    1024, 1025, 1280, 1281, 1026, 1027, 1282, 1283,
    1536, 1537, 1792, 1793, 1538, 1539, 1794, 1795,
    1028, 1029, 1284, 1285, 1030, 1031, 1286, 1287,
    1540, 1541, 1796, 1797, 1542, 1543, 1798, 1799,
    8, 9, 264, 265, 10, 11, 266, 267,
    520, 521, 776, 777, 522, 523, 778, 779,
    12, 13, 268, 269, 14, 15, 270, 271,
    524, 525, 780, 781, 526, 527, 782, 783,
    1032, 1033, 1288, 1289, 1034, 1035, 1290, 1291,
    1544, 1545, 1800, 1801, 1546, 1547, 1802, 1803,
    1036, 1037, 1292, 1293, 1038, 1039, 1294, 1295,
    1548, 1549, 1804, 1805, 1550, 1551, 1806, 1807,
    2048, 2049, 2304, 2305, 2050, 2051, 2306, 2307,
    2560, 2561, 2816, 2817, 2562, 2563, 2818, 2819,
    2052, 2053, 2308, 2309, 2054, 2055, 2310, 2311,
    2564, 2565, 2820, 2821, 2566, 2567, 2822, 2823,
    3072, 3073, 3328, 3329, 3074, 3075, 3330, 3331,
    3584, 3585, 3840, 3841, 3586, 3587, 3842, 3843,
    3076, 3077, 3332, 3333, 3078, 3079, 3334, 3335,
    3588, 3589, 3844, 3845, 3590, 3591, 3846, 3847,
    2056, 2057, 2312, 2313, 2058, 2059, 2314, 2315,
    2568, 2569, 2824, 2825, 2570, 2571, 2826, 2827,
    2060, 2061, 2316, 2317, 2062, 2063, 2318, 2319,
    2572, 2573, 2828, 2829, 2574, 2575, 2830, 2831,
    3080, 3081, 3336, 3337, 3082, 3083, 3338, 3339,
    3592, 3593, 3848, 3849, 3594, 3595, 3850, 3851,
    3084, 3085, 3340, 3341, 3086, 3087, 3342, 3343,
    3596, 3597, 3852, 3853, 3598, 3599, 3854, 3855,
]

@doc raw"""
    pix2xyfRing(resol::Resolution, ipix) :: (Int, Int, Int)

Convert a pixel number into (x, y, face), using RING ordering."""
function pix2xyfRing(resol::Resolution, ipix)
    if ipix < resol.ncap
        # North polar cap
        iring = (1 + isqrt(2 * ipix - 1)) >> 1
        iphi = ipix - 2 * iring * (iring - 1)
        kshift = 0
        nr = iring
        facenum = (iphi - 1) ÷ nr
    elseif ipix < (resol.numOfPixels - resol.ncap)
        # Equatorial region
        ip = ipix - 1 - resol.ncap
        tmp = (resol.order ≥ 0) ? (ip >> (resol.order + 2)) : (ip ÷ resol.nsideTimesFour)
        iring = tmp + resol.nside
        iphi = ip - tmp * resol.nsideTimesFour + 1
        kshift = (iring + resol.nside) & 1
        nr = resol.nside
        ire = iring - resol.nside + 1
        irm = resol.nsideTimesTwo + 2 - ire
        ifm = iphi - ire ÷ 2 + resol.nside - 1
        ifp = iphi - irm ÷ 2 + resol.nside - 1
        if resol.order ≥ 0
            ifm >>= resol.order
            ifp >>= resol.order
        else
            ifm ÷= resol.nside
            ifp ÷= resol.nside
        end

        facenum = if ifp == ifm
            ifp | 4
        else
            if ifp < ifm
                ifp
            else
                ifm + 8
            end
        end
    else
        # South polar cap
        ip = resol.numOfPixels - ipix + 1
        iring = (1 + isqrt(2 * ip - 1)) >> 1
        iphi = 4 * iring + 1 - (ip - 2 * iring * (iring - 1))
        kshift = 0
        nr = iring
        iring = 2 * resol.nsideTimesTwo - iring
        facenum = 8 + (iphi - 1) ÷ nr
    end

    irt = iring - JRLL[facenum+1] * resol.nside + 1
    ipt = 2 * iphi - JPLL[facenum+1] * nr - kshift - 1
    ipt ≥ resol.nsideTimesTwo && (ipt -= 8 * resol.nside)

    ((ipt - irt) >> 1, (-ipt - irt) >> 1, facenum)
end

@doc raw"""
    xyf2pixRing(resol::Resolution, ix, iy, facenum) :: Int

Convert (x, y, face) into a pixel number, using RING ordering."""
function xyf2pixRing(resol::Resolution, ix, iy, facenum)
    jr = JRLL[facenum+1] * resol.nside - ix - iy - 1

    ringinfo = getringinfo(resol, jr, full = false)
    nr = ringinfo.numOfPixels >> 2
    kshift = ringinfo.shifted ? 0 : 1
    jp = (JPLL[facenum+1] * nr + ix - iy + 1 + kshift) ÷ 2
    jp < 1 && (jp += resol.nsideTimesFour)

    (ringinfo.firstPixIdx - 1) + jp
end

function spreadbits(v::Int)
    (
        Int64(UTAB[v&0xff+1]) |
        Int64(UTAB[(v>>8)&0xff+1] << 16) |
        Int64(UTAB[(v>>16)&0xff+1] << 32) |
        Int64(UTAB[(v>>24)&0xff+1] << 48)
    )
end

function compress_bits(v::Int32)
    raw = (v & 0x5555) | ((v & 0x55550000) >> 15)
    CTAB[raw&0xff+1] | (CTAB[raw>>8+1] << 4)
end

function compress_bits(v::Int)
    raw = v & 0x5555555555555555
    raw |= raw >> 15
    return (
        CTAB[raw&0xff+1] |
        (CTAB[(raw>>8)&0xff+1] << 4) |
        (CTAB[(raw>>32)&0xff+1] << 16) |
        (CTAB[(raw>>40)&0xff+1] << 20)
    )
end

@doc raw"""
    pix2xyfNest(resol::Resolution, ipix) :: (Int, Int, Int)

Convert a pixel number into (x, y, face), using NESTED ordering."""
function pix2xyfNest(resol::Resolution, ipix)
    pix = (ipix - 1) & (resol.pixelsPerFace - 1)
    (compress_bits(pix), compress_bits(pix ÷ 2), (ipix - 1) >> (2 * resol.order))
end

@doc raw"""
    xyf2pixNest(resol::Resolution, ix, iy, facenum) :: Int

Convert (x, y, face) into a pixel number, using NESTED ordering."""
function xyf2pixNest(resol::Resolution, ix, iy, facenum)
    facenum << (2 * resol.order) + spreadbits(ix) + 2spreadbits(iy) + 1
end

@doc raw"""
    ring2nest(resol::Resolution, ipix) :: Int

Convert the number of a pixel from RING to NESTED scheme."""
ring2nest(resol::Resolution, ipix) = xyf2pixNest(resol, pix2xyfRing(resol, ipix)...)

@doc raw"""
    nest2ring(resol::Resolution, ipix) :: Int

Convert the number of a pixel from NESTED to RING scheme."""
nest2ring(resol::Resolution, ipix) = xyf2pixRing(resol, pix2xyfNest(resol, ipix)...)
