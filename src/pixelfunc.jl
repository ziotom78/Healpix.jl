################################################################################
# Pixel functions

################################################################################

function calcNestPosForEquator(resol, z, z_abs, scaled_phi, sintheta, havesintheta)
    temp1 = resol.nside * (0.5 + scaled_phi)
    temp2 = resol.nside * (z * 0.75)
    jp = floor(Int, temp1 - temp2)
    jm = floor(Int, temp1 + temp2)

    ifp = jp >> resol.order
    ifm = jm >> resol.order

    local face_num
    face_num = if ifp == ifm
        ifp | 4
    elseif ifp < ifm
        ifp
    else
        ifm + 8
    end

    ix = mod(jm, resol.nside)
    iy = resol.nside - mod(jp, resol.nside) - 1

    xyf2pixNest(resol, ix, iy, face_num)
end

################################################################################

function calcNestPosForPole(resol, z, z_abs, scaled_phi, sintheta, havesintheta)
    ntt = floor(Int, scaled_phi)
    if ntt >= 4
        ntt = 3
    end

    tp = scaled_phi - ntt
    tmp = if (z_abs < 0.99) || (! havesintheta)
        resol.nside * sqrt(3 * (1 - z_abs)) # in ]0,1]
    else
        resol.nside * sintheta / sqrt((1 + z_abs) / 3)
    end

    jp = floor(Int, tp * tmp)
    jm = floor(Int, (1 - tp) * tmp)

    # Clip jp and jm
    jp = jp < resol.nside - 1 ? jp : resol.nside - 1
    jm = jm < resol.nside - 1 ? jm : resol.nside - 1

    if z >= 0
        xyf2pixNest(resol, resol.nside - jm - 1, resol.nside - jp - 1, ntt)
    else
        xyf2pixNest(resol, jp, jm, ntt + 8)
    end
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

    nside = resol.nside
    local ix, iy, face_num

    z = cos(theta)
    z_abs = abs(z)
    scaled_phi = mod2pi(phi) / (π / 2) # in [0,4[

    (sintheta, havesintheta) = if (theta < 0.01) || (theta > 3.14159 - 0.01)
        (sin(theta), true)
    else
        (0.0, false)
    end

    if 3z_abs ≤ 2
        calcNestPosForEquator(resol, z, z_abs, scaled_phi, sintheta, havesintheta)
    else
        calcNestPosForPole(resol, z, z_abs, scaled_phi, sintheta, havesintheta)
    end
end

################################################################################

function calcRingPosForEquator(resol::Resolution, z, z_abs, tt)

    jp = floor(Integer, resol.nside * (0.5 + tt - z * 0.75))
    jm = floor(Integer, resol.nside * (0.5 + tt + z * 0.75))

    ir = resol.nside + 1 + jp - jm
    kshift = (mod(ir, 2) == 0) ? 1 : 0

    nl4 = resol.nsideTimesFour

    local ip = div(jp + jm - resol.nside + kshift + 1, 2) + 1
    if ip > nl4
        ip = ip - nl4
    end

    resol.ncap + nl4 * (ir - 1) + ip
end

################################################################################

function calcRingPosForPole(resol::Resolution, z, z_abs, tt)

    tp = tt - floor(tt)
    tmp = sqrt(3.0 * (1.0 - z_abs))

    jp = floor(Integer, resol.nside * tp * tmp)
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
    (0 ≤ theta ≤ π) || throw(DomainError(theta, "Invalid value of theta"))

    sintheta = sin(theta)
    return (sintheta * cos(phi), sintheta * sin(phi), cos(theta))
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
    phi = atan(y, x)
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

    z = cos(theta)
    z_abs = abs(z)
    scaled_phi = mod2pi(phi) / (π / 2) # in [0,4[

    if z_abs ≤ 2 // 3
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

    jrll = [2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4]
    jpll = [1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7]

    floatNside = Float64(resol.nside)
    fact1 = 1.0 / (3.0 * floatNside^2)
    fact2 = 2.0 / (3.0 * floatNside)

    # face number in {0,11} and pixel number within the face
    faceNum, ipf = divrem(pixel - 1, resol.pixelsPerFace)

    ip_trunc, ip_low = divrem(ipf, 1024)
    ip_hi, ip_med = divrem(ip_trunc, 1024)

    ix = 1024 * pix2x[ip_hi+1] + 32 * pix2x[ip_med+1] + pix2x[ip_low+1]
    iy = 1024 * pix2y[ip_hi+1] + 32 * pix2y[ip_med+1] + pix2y[ip_low+1]

    # Transforms this in (horizontal, vertical) coordinates
    jrt = ix + iy # 'vertical' in {0,2*(nside-1)}
    jpt = ix - iy # 'horizontal' in {-nside+1,nside-1}

    jr = jrll[faceNum+1] * resol.nside - jrt - 1
    local nr = resol.nside # Equatorial region (the most frequent)
    z = (2resol.nside - jr) * fact2
    local kshift = Int(mod(jr - resol.nside, 2))
    if jr < resol.nside
        # North polar cap
        nr = jr
        z = 1.0 - nr^2 * fact1
        kshift = 0
    elseif jr > 3 * resol.nside
        # South polar cap
        nr = resol.nsideTimesFour - jr
        z = -1.0 + nr^2 * fact1
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
    pix2ringpos(resol::Resolution, pixel)

Given the (1-based) index of a pixel in a Healpix map in ring
order, return a pair of numbers (n, i, j) whose meaning is the following:

- `n` can be one of the symbols `:northcap`, `:equator`, or
  `:southcap`, representing the region of the sky
- `i` is the ring index, from 1 to 4NSIDE - 1
- `j` is the pixel-in-ring index
"""
function pix2ringpos(resol::Resolution, pixel)
    # Any reference to equations in this routine refers to Gorski et al. (2005)
    if pixel ≤ resol.ncap
        # North polar cap

        p_h = pixel / 2 # Defined in Gorsky et al. (2005) before Eq. (2)
        floor_p_h = floor(p_h)
        # Eq. (2)
        i = floor(Int, sqrt(p_h - sqrt(floor_p_h))) + 1
        # Eq. (3)
        j = pixel - 2i * (i - 1)

        (:northcap, i, j)
    elseif pixel ≤ resol.nsideTimesTwo * (5resol.nside + 1)
        ip = pixel - resol.ncap - 1
        # Eq. (6) - ring counts from the North pole; resol.nside is
        # the number of pixels in the North Polar ring
        i = floor(Int, ip / resol.nsideTimesFour) + resol.nside
        # Eq. (7) - zero-based index of the pixel within this ring
        j = Int(mod(ip, resol.nsideTimesFour)) + 1
        (:equator, i, j)
    else
        # South polar cap

        ip = resol.numOfPixels - pixel + 1
        p_h = ip / 2
        floor_p_h = floor(p_h)
        i = floor(Int, sqrt(p_h - sqrt(floor_p_h))) + 1 # counted from S. pole
        j = Int(4 * i + 1 - (ip - 2i * (i - 1)))

        (:southcap, i, j)
    end
end

################################################################################

"""
    pix2angRing(resol::Resolution, pixel) -> (Float64, Float64)

Given the (1-based) index of a pixel in a Healpix map in ring
order, return a pair containing the (`colatitude`, `longitude`) angles
corresponding to its center, both expressed in radians.
"""
function pix2angRing(resol::Resolution, pixel)
    fact1 = 1.5 * resol.nside
    fact2 = 3.0 * resol.pixelsPerFace

    # Any reference to equations in this routine refers to Gorski et al. (2005)
    cap, i, j = pix2ringpos(resol, pixel)
    if cap == :northcap
        # Colatitude: Eq. (4); longitude: Eq. (5)
        return (acos(1 - i^2 / fact2), (Float64(j) - 0.5) * π / (2i))
    elseif pixel ≤ resol.nsideTimesTwo * (5resol.nside + 1)
        # Equatorial belt

        # Eq. (9) - this equals 1 if i + resol.nside is odd, 1/2
        # otherwise. It is used to convert j into a longitude (since
        # pixel centers in odd rings are shifted with respect to
        # centers in even rings)
        s_half = 0.5 * (1 + mod(Float64(i + resol.nside), 2))

        # Colatitude: Eq. (8) in disguise, latitude: Eq. (9)
        return (
            acos((resol.nsideTimesTwo - i) / fact1),
            (Float64(j) - s_half) * π / (2resol.nside),
        )
    else
        # South Polar cap

        # The pixels in this cap are handled like the ones in the
        # North Polar cap, except that we must flip the value of "ip".
        return (acos(-1 + i^2 / fact2), (float(j) - 0.5) * π / (2i))
    end
end

################################################################################

vec2pixNest(res::Resolution, x, y, z) = ang2pixNest(res, vec2ang(x, y, z)...)
vec2pixRing(res::Resolution, x, y, z) = ang2pixRing(res, vec2ang(x, y, z)...)
pix2vecNest(res::Resolution, pixel) = ang2vec(pix2angNest(res, pixel)...)
pix2vecRing(res::Resolution, pixel) = ang2vec(pix2angRing(res, pixel)...)

################################################################################

function pix2locRing(res::Resolution, ipix)

    pix = ipix - 1
    if pix < res.ncap
        # North polar cap
        (pix < 0) && println("res = $res, ipix = $ipix, pix = $pix")
        iring = (1 + round(Int, sqrt(1 + 2pix), RoundDown)) >> 1
        iphi = (pix + 1) - 2iring * (iring - 1)

        tmp = iring^2 * res.fact2
        z = 1.0 - tmp
        (sth, have_sth) = if z > 0.99
            (sqrt(tmp * (2 - tmp)), true)
        else
            (0.0, false)
        end

        phi = (iphi - 0.5) * π / 2 / iring

        return (z, phi, sth, have_sth)
    end

    if pix < (res.numOfPixels - res.ncap)
        # Equatorial region
        nl4 = 4 * res.nside
        ip = pix - res.ncap
        tmp = (res.order >= 0) ? (ip >> (res.order + 2)) : (ip / nl4)
        iring = tmp + res.nside
        iphi = ip - nl4 * tmp + 1

        fodd = (((iring + res.nside) & 1) != 0) ? 1.0 : 0.5

        z = (2 * res.nside - iring) * res.fact1
        phi = (iphi - fodd) * π * 0.75 * res.fact1

        return (z, phi, 0.0, false)
    end

    # South polar cap

    ip = res.numOfPixels - pix

    # Counted from south pole
    iring = (1 + round(Int, sqrt(2 * ip - 1), RoundDown)) >> 1
    iphi = 4iring + 1 - (ip - 2iring * (iring - 1))

    tmp = iring^2 * res.fact2
    z = tmp - 1.0

    (sth, have_sth) = if z < -0.99
        (sqrt(tmp * (2 - tmp)), true)
    else
        (0.0, false)
    end

    phi = (iphi - 0.5) * π / 2 / iring

    (z, phi, sth, have_sth)
end

"""
    pix2zphiRing(res::Resolution, pix) -> (z, phi)

Compute the angular coordinates `z = cos(θ), ϕ` of the center of the
pixel with number `pix`, assuming the `RING` numbering scheme for
pixels. *Caution:* this method is inaccurate near the poles at high
resolutions.

"""
function pix2zphiRing(res::Resolution, pix)
    (z, phi, _, _) = pix2locRing(res, pix)
    (z, phi)
end

function pix2locNest(res::Resolution, pix)

    z = 0.0
    phi = 0.0
    sth = 0.0
    have_sth = false

    (ix, iy, face_num) = pix2xyfNest(res, pix)

    jr = (JRLL[face_num+1] << res.order) - ix - iy - 1

    nr = 0

    if res.nside <= jr <= 3res.nside
        nr = res.nside
        z = (2 * res.nside - jr) * res.fact1
    else
        nr = (jr < res.nside) ? jr : (res.nsideTimesFour - jr)

        tmp = nr^2 * res.fact2
        z = (jr < res.nside) ? (1.0 - tmp) : (tmp - 1.0)
        if abs(z) > 0.99
            sth = sqrt(tmp * (2 - tmp))
            have_sth = true
        end
    end

    tmp = round(Int, JPLL[face_num+1]) * nr + ix - iy

    (tmp < 0) && (tmp += 8nr)

    phi = (nr == res.nside) ? (3π / 8 * tmp * res.fact1) : ((π / 4 * tmp) / nr)

    (z, phi, sth, have_sth)
end

"""
    pix2zphiNest(res::Resolution, pix) -> (z, phi)

Compute the angular coordinates `z = cos(θ), ϕ` of the center of the
pixel with number `pix`, assuming the `NEST` numbering scheme for
pixels. *Caution:* this method is inaccurate near the poles at high
resolutions.

"""
function pix2zphiNest(res::Resolution, pix)
    (z, phi, _, _) = pix2locNest(res, pix)
    (z, phi)
end

################################################################################

"""
    ringAbove(res::Resolution, z) -> (ring_number)

Return the number of the next ring to the north of `z = cos(θ)`.
If `z` lies north of all rings, the function returns 0.

"""
function ringAbove(res::Resolution, z)
    az = abs(z)
    3az <= 2 && return round(Int, res.nside * (4 - 3z) * 0.5)
    iring = round(Int, res.nside * sqrt(3(1 - az)))

    (z > 0) ? iring : (4 * res.nside - iring - 1)
end
