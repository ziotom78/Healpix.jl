################################################################################
# Pixel functions

cachesintheta(theta) = ((theta < 0.01) || (theta > 3.14159 - 0.01)) ? (sin(theta), true) : (0.0, false)


function pix2locRing(resol::Resolution, pixel)
    # This is the same as the C++ code (zero-based index)
    pix = pixel - 1
    (sintheta, havesintheta) = (0.0, false)

    if pix < resol.ncap
        iring = (1 + floor(Int, isqrt(1 + 2 * pix))) >> 1
        iphi = (pix + 1) - 2 * iring * (iring - 1)

        tmp = iring * iring * resol.fact2
        z = 1 - tmp
        if z > 0.99
            (sintheta, havesintheta) = (sqrt(tmp * (2 - tmp)), true)
        end

        phi = (iphi - 0.5) * (π/2) / iring
    elseif pix < (resol.numOfPixels - resol.ncap)
        nl4 = resol.nside * 4
        ip = pix - resol.ncap
        tmp = (resol.order >= 0) ? (ip >> (resol.order + 2)) : div(ip, nl4)
        iring = tmp + resol.nside
        iphi = ip - nl4 * tmp + 1
        fodd = (((iring + resol.nside) & 1) != 0) ? 1.0 : 0.5

        z = (2 * resol.nside - iring) * resol.fact1
        phi = (iphi - fodd) * pi * 0.75 * resol.fact1
    else # South polar cap
        ip = resol.numOfPixels - pix
        iring = (1 + floor(Int, isqrt(2ip - 1))) >> 1
        iphi = 4iring + 1 - (ip - 2iring * (iring - 1))

        tmp = (iring * iring) * resol.fact2
        z = tmp - 1
        if z < 0.99
            (sintheta, havesintheta) = (sqrt(tmp * (2 - tmp)), true)
        end

        phi = (iphi - 0.5) * (π/2) / iring
    end

    (z, phi, sintheta, havesintheta)
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

    (0 ≤ theta ≤ π) || throw(DomainError(theta, "Invalid value of theta"))

    nside = resol.nside
    local ix, iy, face_num

    z = cos(theta)
    z_abs = abs(z)
    scaled_phi = mod2pi(phi) / (π / 2) # in [0,4[

    (sintheta, havesintheta) = cachesintheta(theta)

    if 3z_abs ≤ 2
        calcNestPosForEquator(resol, z, z_abs, scaled_phi, sintheta, havesintheta)
    else
        calcNestPosForPole(resol, z, z_abs, scaled_phi, sintheta, havesintheta)
    end
end

################################################################################

function calcRingPosForEquator(resol::Resolution, z, z_abs, tt, sintheta, havesintheta)

    nl4 = 4 * resol.nside
    temp1 = resol.nside * (0.5 + tt)
    temp2 = resol.nside * z * 0.75

    jp = floor(Int, temp1 - temp2)  # Index of ascending edge line
    jm = floor(Int, temp1 + temp2)  # Index of descending edge line

    ir = resol.nside + 1 + jp - jm  # Ring number counted from z=2/3 in {1, 2n+1}
    kshift = 1 - (ir & 1)           # kshift = 1 if ir is even, 0 otherwise

    t1 = jp + jm - resol.nside + kshift + 1 + 2nl4
    ip = (resol.order > 0) ? ((t1 >> 1) & (nl4 - 1)) : mod(t1 >> 1, nl4)

    resol.ncap + (ir - 1) * nl4 + ip + 1
end

################################################################################

function calcRingPosForPole(resol::Resolution, z, z_abs, tt, sintheta, havesintheta)

    tp = tt - floor(Int, tt)

    tmp = if (z_abs < 0.99) || (! havesintheta)
        resol.nside * sqrt(3 * (1 - z_abs))
    else
        resol.nside * sintheta / sqrt((1 + z_abs) / 3)
    end

    jp = floor(Int, tp * tmp)       # Increasing edge line index
    jm = floor(Int, (1 - tp) * tmp) # Decreasing edge line index

    ir = jp + jm + 1
    ip = floor(Int, tt * ir)

    if z > 0
        2 * ir * (ir - 1) + ip + 1
    else
        resol.numOfPixels - 2 * ir * (ir + 1) + ip + 1
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
    theta = atan(sqrt(x^2 + y^2), z)
    phi = atan(y, x)
    if phi < 0
        phi += 2π
    end

    (theta, phi)
end

################################################################################

"""
    zphi2pixRing(resol::Resolution, theta, phi) -> Integer

Return the index of the pixel which contains the point with
coordinates (`theta`, the colatitude, and `phi`, the longitude), in
radians, for a Healpix map with pixels in ring order. Note that pixel
indexes are 1-based (this is Julia)!
"""
function zphi2pixRing(resol::Resolution, z, phi)

    z_abs = abs(z)

    # We do not used mod2pi because we want 1-1 match with C++ code
    scaled_phi = mod(phi * 2 / π, 4.0e0)

    if 3 * z_abs ≤ 2
        calcRingPosForEquator(resol, z, z_abs, scaled_phi, 0.0, false)
    else
        calcRingPosForPole(resol, z, z_abs, scaled_phi, 0.0, false)
    end
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

    (0 ≤ theta ≤ π) || throw(DomainError(theta, "Invalid value of theta"))

    z = cos(theta)
    z_abs = abs(z)

    # We do not use mod2pi because we want 1-1 match with C++ code
    scaled_phi = mod(phi * 2 / π, 4.0e0)

    (sintheta, havesintheta) = cachesintheta(theta)

    if 3 * z_abs ≤ 2
        calcRingPosForEquator(resol, z, z_abs, scaled_phi, sintheta, havesintheta)
    else
        calcRingPosForPole(resol, z, z_abs, scaled_phi, sintheta, havesintheta)
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

    (z, phi, sintheta, havesintheta) = pix2locNest(resol, pixel)

    if havesintheta
        (atan(sintheta, z), phi)
    else
        (acos(z), phi)
    end
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

    (z, phi, sintheta, havesintheta) = pix2locRing(resol, pixel)

    if havesintheta
        (atan(sintheta, z), phi)
    else
        (acos(z), phi)
    end
end

################################################################################

vec2pixNest(res::Resolution, x, y, z) = ang2pixNest(res, vec2ang(x, y, z)...)
vec2pixRing(res::Resolution, x, y, z) = ang2pixRing(res, vec2ang(x, y, z)...)
pix2vecNest(res::Resolution, pixel) = ang2vec(pix2angNest(res, pixel)...)
pix2vecRing(res::Resolution, pixel) = ang2vec(pix2angRing(res, pixel)...)

################################################################################

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
    3az <= 2 && return floor(Int, res.nside * (4 - 3z) * 0.5)
    iring = floor(Int, res.nside * sqrt(3(1 - az)))

    (z > 0) ? iring : (4 * res.nside - iring - 1)
end


@doc raw"""
    ring2z(res::Resolution, ring) -> z

Return the value of `z = \cos(\theta)` for the given ring.

"""
function ring2z(res::Resolution, ring)
    if ring < res.nside
        return 1 - ring^2 * res.fact2
    elseif ring <= 3 * res.nside
        return (res.nsideTimesTwo - ring) * res.fact1
    end

    ring = res.nsideTimesFour - ring
    ring^2 * res.fact2 - 1
end


function set_z_phi(z, phi)
    sintheta = sqrt((1 - z) * (1 + z))
    (sintheta * cos(phi), sintheta * sin(phi), z)
end


v_angle(v1, v2) = atan(norm(SVector(v1) × SVector(v2)), SVector(v1) ⋅ SVector(v2))


function max_pixrad(res::Resolution)
    va = set_z_phi(2/3, π / res.nsideTimesFour)
    t1 = (1 - 1 / res.nside)^2
    vb = set_z_phi(1 - t1/3, 0)
    v_angle(va, vb)
end

function max_pixrad(res::Resolution, ring)
    (ring >= res.nsideTimesTwo) && (ring = res.nsideTimesFour - ring)
    z = ring2z(res, ring)
    z_up = ring2z(res, ring - 1)
    uppos = set_z_phi(z_up, 0)

    if ring <= res.nside
        mypos = set_z_phi(z, π / 4ring)
        v1 = v_angle(mypos, uppos)
        (ring != 1) && (return v1)
        uppos = set_z_phi(ring2z(res, ring + 1), π / (4 * min(res.nside, ring + 1)))
        return max(v1, v_angle(mypos, uppos))
    end

    mypos = set_z_phi(z, 0)
    vdist = v_angle(mypos, uppos)
    hdist = sqrt(1 - z^2) * π / res.nsideTimesFour
    max(hdist, vdist)
end

@doc raw"""
    max_pixrad(res::Resolution, ring)
    max_pixrad(res::Resolution)

Return the maximum angular distance (in radians) between a pixel
center and any of its corners. If `ring` is specified, the result
applies to all the pixels of the given ring; otherwise, all the
pixels on the sphere are considered.

"""
max_pixrad


################################################################################

function loc2vec(z, phi, sintheta, have_sintheta, T::Type{<:Real})
    if ! have_sintheta
        sintheta = sqrt((1 - z) * (1 + z))
    end

    SVector{3, T}(sintheta * cos(phi), sintheta * sin(phi), z)
end

@doc """
    boundariesRing!(resol::Resolution, pix, step, buf::Matrix{T}) where {T <: Real}

See the documentation for [`boundariesRing`](@ref).

"""
function boundariesRing!(resol::Resolution, pix, step, buf::Matrix{T}) where {T <: Real}

    (ix, iy, face) = pix2xyfRing(resol, pix)
    dc = 0.5 / resol.nside
    (xc, yc) = ((ix + 0.5) / resol.nside, (iy + 0.5) / resol.nside)
    d = 1.0 / (step * resol.nside)

    for i in 0:(step - 1)
        for (incr, x, y) in zip(
            (1, 1 + step, 1 + 2step, 1 + 3step),
            (
                xc + dc - i * d,
                xc - dc,
                xc - dc + i * d,
                xc + dc,
            ),
            (
                yc + dc,
                yc + dc - i * d,
                yc - dc,
                yc - dc + i * d,
            ),
        )
            (z, phi, sintheta, have_sintheta) = xyf2loc(x, y, face)
            buf[i + incr, :] = loc2vec(z, phi, sintheta, have_sintheta, T)
        end
    end
end

@doc raw"""
    boundariesRing(resol::Resolution, pix, step, T::Type{<:Real})
    boundariesRing!(resol::Resolution, pix, step, buf::Matrix{T}) where {T <: Real}

Compute a set of directions (3D vectors) along the boundaries of a given
pixel in the RING scheme at some resolution.

The function `boundariesRing` returns a ``N \times 3`` matrix of type
`T` containing ``N`` vectors pointing towards the border of the pixel
with index `pix` in RING scheme. Each edge of the pixel contains
`step` points, and, as every pixel has a diamond-like shape with four
edges, the number ``N`` is equal to `4 * step`.

If you plan to call this function again and again, you should allocate
your own matrix with the results and call `boundariesRing!`, which
accepts the parameter `buf` where the result will be stored. The shape
of this matrix must be `(4step, 3)`, for instance, and its element
can be left undefined:

```julia
step = 10
buf = Matrix{Float64}(undef, 4step, 3)
boundariesRing!(res, pixidx, step, buf)  # This sets `buf`
```

# Examples

Here we show how to use `boundariesRing` and `boundariesRing!`
together to avoid allocating a matrix twice.

```julia-repl
julia> matr = boundariesRing(Resolution(16), 534, 2, Float16)
8×3 Matrix{Float16}:
 0.3535  -0.6123  0.707
 0.3525  -0.6353  0.687
 0.3513  -0.657   0.6665
 0.3762  -0.664   0.646
 0.4014  -0.6694  0.625
 0.4084  -0.645   0.646
 0.414   -0.6196  0.6665
 0.3843  -0.6167  0.687

julia> boundariesRing!(Resolution(16), 535, 2, matr)  # Reuse `matr`

julia> matr
8×3 Matrix{Float16}:
 0.4158  -0.5723  0.707
 0.415   -0.596   0.687
 0.414   -0.6196  0.6665
 0.4397  -0.624   0.646
 0.465   -0.627   0.625
 0.4697  -0.602   0.646
 0.473   -0.576   0.6665
 0.4446  -0.5747  0.687
```

"""
function boundariesRing(resol::Resolution, pix, step, T::Type{<:Real})
    buf = Matrix{T}(undef, 4step, 3)
    boundariesRing!(resol, pix, step, buf)
    buf
end

#################################################################
"""
    getEquatorIdx(nside::Integer)
    getEquatorIdx(res::Resolution)

    Computes the ring index of the equator in a map of `Resolution` `res` or
    NSIDE parameter given by `nside`.

"""
getEquatorIdx(nside::Integer) = 2*nside
getEquatorIdx(res::Resolution) = getEquatorIdx(res.nside)

######################################################

"""
    Create an array of the colatitude in radians (theta) of each ring index in `rings` for a map with resolution `res`.

    If no `rings` array is passed, the computation is performed on all the rings deducted from `res`.

    If an integer `ring` is passed, a single Float value of the colatitude is returned.
"""
function ring2theta(rings::Vector{I}, res::Resolution) where {I<:Integer}
    theta = Vector{Float64}(undef, length(rings))
    ringinfo = RingInfo(0, 0, 0, 0, 0)
    @inbounds for i in 1:length(rings)
        getringinfo!(res, rings[i], ringinfo; full=true)
        theta[i] = ringinfo.colatitude_rad
    end
    theta
end
ring2theta(res::Resolution) = ring2theta(Vector{Int}(1:res.nsideTimesFour-1), res)
#single ring passed as int:
ring2theta(ring::Integer, res::Resolution) = getringinfo(res, ring; full=true).colatitude_rad
    
 
# auxiliary quantites needed to obtain neughbouring pixels, implementation following  https://healpix.jpl.nasa.gov/html/Healpix_cxx/healpix__base2_8cc-source.html    
_xoffset = @SVector [-1, -1, 0, 1, 1, 1, 0, -1]
_yoffset = @SVector [0, 1, 1, 1, 0, -1, -1, -1]
_facearray = @SMatrix [[8, 9, 10, 11, -1, -1, -1, -1, 10, 11, 8, 9];;   # S
    [5, 6, 7, 4, 8, 9, 10, 11, 9, 10, 11, 8];;   # SE
    [-1, -1, -1, -1, 5, 6, 7, 4, -1, -1, -1, -1];;   # E
    [4, 5, 6, 7, 11, 8, 9, 10, 11, 8, 9, 10];;   # SW
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11];;   # center
    [1, 2, 3, 0, 0, 1, 2, 3, 5, 6, 7, 4];;   # NE
    [-1, -1, -1, -1, 7, 4, 5, 6, -1, -1, -1, -1];;   # W
    [3, 0, 1, 2, 3, 0, 1, 2, 4, 5, 6, 7];;   # NW
    [2, 3, 0, 1, -1, -1, -1, -1, 0, 1, 2, 3]]    # N

_swaparray = @SMatrix [[0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3];;   # S
    [0, 0, 0, 0, 0, 0, 0, 0, 6, 6, 6, 6];;   # SE
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;   # E
    [0, 0, 0, 0, 0, 0, 0, 0, 5, 5, 5, 5];;   # SW
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;   # center
    [5, 5, 5, 5, 0, 0, 0, 0, 0, 0, 0, 0];;   # NE
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];;   # W
    [6, 6, 6, 6, 0, 0, 0, 0, 0, 0, 0, 0];;   # NW
    [3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0]]    # N

"""
    getAllNeighbours(pix::Integer, resol::Resolution[; nest::Bool])

Obtain all the neigbouring pixels of `pix`.
"""
function getAllNeighbours(pix::Integer, resol::Resolution; nest::Bool=false)
    result = @MVector zeros(Int, 8)

    if nest
        ix, iy, face_num = pix2xyfNest(resol, pix)
    else
        ix, iy, face_num = pix2xyfRing(resol, pix)
    end

    nside_ = resol.nside
    nsm1 = nside_ - 1


    if (ix > 0) && (ix < nsm1) && (iy > 0) && (iy < nsm1)
        for m in 1:8
            if nest
                result[m] = xyf2pixNest(resol, ix + _xoffset[m], iy + _yoffset[m], face_num)
            else
                result[m] = xyf2pixRing(resol, ix + _xoffset[m], iy + _yoffset[m], face_num)
            end
        end
    else
        for i in 1:8
            x = ix + _xoffset[i]
            y = iy + _yoffset[i]
            nbnum = 4
            if x < 0
                x += nside_
                nbnum -= 1
            elseif x >= nside_
                x -= nside_
                nbnum += 1
            end
            if y < 0
                y += nside_
                nbnum -= 3
            elseif y >= nside_
                y -= nside_
                nbnum += 3
            end

            f = _facearray[face_num+1, nbnum+1]
            if f >= 0
                if _swaparray[face_num+1, nbnum+1] & 1 != 0
                    x = nside_ - x - 1
                end
                if _swaparray[face_num+1, nbnum+1] & 2 != 0
                    y = nside_ - y - 1
                end
                if _swaparray[face_num+1, nbnum+1] & 4 != 0
                    x, y = y, x
                end
                if nest
                    result[i] = xyf2pixNest(resol, x, y, f)
                else
                    result[i] = xyf2pixRing(resol, x, y, f)
                end
            else
                result[i] = -1
            end
        end
    end
    return result
end



