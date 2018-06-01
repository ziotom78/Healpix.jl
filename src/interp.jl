# Ring information and interpolation functions

doc"""
    RingInfo

Information about a ring of pixels, i.e., the set of pixels on a Healpix map
which have the same colatitude. The type is "mutable", so that one object can begin
reused many times without further memory allocations.    
"""
mutable struct RingInfo
    ring::Int
    firstpixidx::Int
    numofpixels::Int
    colatitude_rad::Float64
    shifted::Bool
end

even(x) = (x & 1) == 0

doc"""
    getringinfo!(resol::Resolution, ring, ringinfo::RingInfo; full=true) :: RingInfo

Fill the RingInfo structure with information about the specified ring.
If `full` is `false`, the field `colatitude_rad` (the most expensive in
terms of computation) is set to `NaN`.
"""
function getringinfo!(resol::Resolution, ring, ringinfo::RingInfo; full=true)
    # In the body of this code, firstidx is zero-based (we switch to 1-based
    # index just in the last statement)

    # "northring" is always in the Northern emisphere
    northring = (ring > resol.nsideTimesTwo) ? (resol.nsideTimesFour - ring) : ring
    if northring < resol.nside
        θ = if full
            tmp = northring^2 * resol.fact2
            atan2(sqrt(tmp * (2 - tmp)), 1 - tmp)
        else
            NaN
        end
        numofpixels = 4 * northring
        shifted = true
        firstpixidx = 2 * northring * (northring - 1)
    else
        θ = full ? acos((resol.nsideTimesTwo - northring) * resol.fact1) : NaN
        numofpixels = resol.nsideTimesFour
        shifted = even(northring - resol.nside)
        firstpixidx = resol.ncap + (northring - resol.nside) * numofpixels
    end

    if northring != ring
        # Southern emisphere
        full && (θ = π - θ)
        firstpixidx = resol.numOfPixels - firstpixidx - numofpixels
    end

    ringinfo.ring = ring
    ringinfo.firstpixidx = firstpixidx + 1
    ringinfo.numofpixels = numofpixels
    ringinfo.colatitude_rad = θ
    ringinfo.shifted = shifted
end

doc"""
    getringinfo(resol::Resolution, ring; kwargs...) :: RingInfo

Return a RingInfo structure containing information about
the specified ring. For the list of accepted keyword arguments,
see getringinfo!.
"""
function getringinfo(resol::Resolution, ring; kwargs...)
    ringinfo = RingInfo(0, 0, 0, 0.0, true)
    getringinfo!(resol, ring, ringinfo, kwargs...)
    ringinfo
end

function ringabove(resol::Resolution, z)
    az = abs(z)
    az <= 2 // 3 && return round(Int, resol.nside * (2 - 1.5 * z), RoundToZero)
    
    iring = round(Int, resol.nside * sqrt(3 * (1 - az)), RoundToZero)
    if z > 0
        iring
    else
        resol.nsideTimesFour - iring - 1
    end
end

function ring2idxw(ringinfo::RingInfo, ϕ)
    shift = ringinfo.shifted ? 0.5 : 0
    # Angle subtended by each pixel along this ring
    δϕ = 2π / ringinfo.numofpixels
    tmp = ϕ / δϕ - shift
    i1 = (tmp < 0) ? round(Int, tmp, RoundToZero) - 1 : round(Int, tmp, RoundToZero)
    w1 = (ϕ - (i1 + shift) * δϕ) / δϕ
    i2 = i1 + 1
    i1 < 0 && (i1 += ringinfo.numofpixels)
    i2 >= ringinfo.numofpixels && (i2 -= ringinfo.numofpixels)

    ([ringinfo.firstpixidx + i1, ringinfo.firstpixidx + i2], [1 - w1, w1])
end

doc"""
    getinterpolRing(resol::Resolution, θ, ϕ) :: (Array{Int,1}, Array{Float64, 1})

Return the indices and the weights of the four neighbour pixels for the given
direction (θ, ϕ) in a map with the specified resolution.
"""
function getinterpolRing(resol::Resolution, θ, ϕ)
    @assert (0 ≤ θ ≤ π)

    z = cos(θ)
    ir1 = ringabove(resol, z)
    ir2 = ir1 + 1

    pix = Array{Int}(4) # Zero-based till the end
    wgt = Array{Float64}(4)
    if ir1 > 0
        ring1 = getringinfo(resol, ir1)
        pix1, w1 = ring2idxw(ring1, ϕ)
        pix[1:2] = pix1 - 1 # Switch from 1-based to 0-based
        wgt[1:2] = w1
    end

    if ir2 < resol.nsideTimesFour
        ring2 = getringinfo(resol, ir2)
        pix2, w2 = ring2idxw(ring2, ϕ)
        pix[3:4] = pix2 - 1 # Switch from 1-based to 0-based
        wgt[3:4] = w2
    end

    if ir1 == 0
        wθ = θ / ring2.colatitude_rad
        wgt[3:4] .*= wθ
        fac = (1 - wθ) / 4
        wgt[:] = [fac, fac, wgt[3] + fac, wgt[4] + fac]
        pix[1] = (pix[3] + 2) & 3
        pix[2] = (pix[4] + 2) & 3
    elseif ir2 == resol.nsideTimesFour
        wθ = (θ - ring1.colatitude_rad) / (π - ring1.colatitude_rad)
        wgt[1:2] .*= 1 - wθ
        fac = wθ / 4
        wgt[:] = [wgt[1] + fac, wgt[2] + fac, fac, fac]
        pix[3] = ((pix[1] + 2) & 3) + resol.numOfPixels - 4
        pix[4] = ((pix[2] + 2) & 3) + resol.numOfPixels - 4
    else
        wθ = (θ - ring1.colatitude_rad) / (ring2.colatitude_rad - ring1.colatitude_rad)
        wgt[1:2] .*= (1 - wθ)
        wgt[3:4] .*= wθ
    end

    (pix + 1, wgt)
end