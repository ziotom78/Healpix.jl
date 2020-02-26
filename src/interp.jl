# Ring information and interpolation functions

export RingInfo

@doc raw"""
    RingInfo

Information about a ring of pixels, i.e., the set of pixels on a Healpix map
which have the same colatitude. The type is "mutable", so that one object can begin
reused many times without further memory allocations.

The list of fields defined in this structure is the following:

- `ring`: an integer index, running from 

- `firstPixIdx`: index of the first pixel (using the `RING` scheme)
  belonging to this ring

- `numOfPixels`: number of consecutive pixels within the ring

- `colatitude_rad`: value of the colatitude for this ring (in radians)

- `shifted`: Boolean flag; it is `true` if the longitude of the first
  pixel in the ring is not zero.

# References

See also [`getringinfo!`](@ref) and [`getringinfo`](@ref).

# Example

```julia
import Healpix
res = Healpix.Resolution(256)

# Show information about ring #10
print(getringinfo(res, 10))
```
"""
mutable struct RingInfo
    ring::Int
    firstPixIdx::Int
    numOfPixels::Int
    colatitude_rad::Float64
    shifted::Bool
end

even(x) = (x & 1) == 0

@doc raw"""
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
            atan(sqrt(tmp * (2 - tmp)), 1 - tmp)
        else
            NaN
        end
        numOfPixels = 4 * northring
        shifted = true
        firstPixIdx = 2 * northring * (northring - 1)
    else
        θ = full ? acos((resol.nsideTimesTwo - northring) * resol.fact1) : NaN
        numOfPixels = resol.nsideTimesFour
        shifted = even(northring - resol.nside)
        firstPixIdx = resol.ncap + (northring - resol.nside) * numOfPixels
    end

    if northring != ring
        # Southern emisphere
        full && (θ = π - θ)
        firstPixIdx = resol.numOfPixels - firstPixIdx - numOfPixels
    end

    ringinfo.ring = ring
    ringinfo.firstPixIdx = firstPixIdx + 1
    ringinfo.numOfPixels = numOfPixels
    ringinfo.colatitude_rad = θ
    ringinfo.shifted = shifted
end

@doc raw"""
    getringinfo(resol::Resolution, ring; kwargs...) :: RingInfo

Return a RingInfo structure containing information about
the specified ring. For the list of accepted keyword arguments,
see getringinfo!.
"""
function getringinfo(resol::Resolution, ring; full=true)
    ringinfo = RingInfo(0, 0, 0, 0.0, true)
    getringinfo!(resol, ring, ringinfo, full=full)
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
    δϕ = 2π / ringinfo.numOfPixels
    tmp = ϕ / δϕ - shift
    i1 = (tmp < 0) ? round(Int, tmp, RoundToZero) - 1 : round(Int, tmp, RoundToZero)
    w1 = (ϕ - (i1 + shift) * δϕ) / δϕ
    i2 = i1 + 1
    i1 < 0 && (i1 += ringinfo.numOfPixels)
    i2 >= ringinfo.numOfPixels && (i2 -= ringinfo.numOfPixels)

    ([ringinfo.firstPixIdx + i1, ringinfo.firstPixIdx + i2], [1 - w1, w1])
end

function getinterpolRing(resol::Resolution, θ, ϕ, pix, weights)
    @assert (0 ≤ θ ≤ π)

    z = cos(θ)
    ir1 = ringabove(resol, z)
    ir2 = ir1 + 1

    # Note that "pix" will be zero-based till the end
    if ir1 > 0
        ring1 = getringinfo(resol, ir1)
        pix1, w1 = ring2idxw(ring1, ϕ)
        pix[1:2] = pix1 .- 1 # Switch from 1-based to 0-based
        weights[1:2] = w1
    end

    if ir2 < resol.nsideTimesFour
        ring2 = getringinfo(resol, ir2)
        pix2, w2 = ring2idxw(ring2, ϕ)
        pix[3:4] = pix2 .- 1 # Switch from 1-based to 0-based
        weights[3:4] = w2
    end

    if ir1 == 0
        wθ = θ / ring2.colatitude_rad
        weights[3:4] .*= wθ
        fac = (1 - wθ) / 4
        weights[:] = [fac, fac, weights[3] + fac, weights[4] + fac]
        pix[1] = (pix[3] + 2) & 3
        pix[2] = (pix[4] + 2) & 3
    elseif ir2 == resol.nsideTimesFour
        wθ = (θ - ring1.colatitude_rad) / (π - ring1.colatitude_rad)
        weights[1:2] .*= 1 - wθ
        fac = wθ / 4
        weights[:] = [weights[1] + fac, weights[2] + fac, fac, fac]
        pix[3] = ((pix[1] + 2) & 3) + resol.numOfPixels - 4
        pix[4] = ((pix[2] + 2) & 3) + resol.numOfPixels - 4
    else
        wθ = (θ - ring1.colatitude_rad) / (ring2.colatitude_rad - ring1.colatitude_rad)
        weights[1:2] .*= (1 - wθ)
        weights[3:4] .*= wθ
    end

    # Fix the order from 0-based to 1-based
    pix .+= 1
end

function getinterpolRing(resol::Resolution, θ, ϕ)
    pix = Array{Int}(undef, 4)
    wgt = Array{Float64}(undef, 4)

    getinterpolRing(resol, θ, ϕ, pix, wgt)

    (pix, wgt)
end

@doc raw"""
    getinterpolRing(resol::Resolution, θ, ϕ) -> (Array{Int,1}, Array{Float64, 1})
    getinterpolRing!(resol::Resolution, θ, ϕ, pix, weights) -> (Array{Int,1}, Array{Float64, 1})

Return the indices and the weights of the four neighbour pixels for
the given direction (θ, ϕ) in a map with the specified resolution.

If provided, the parameters `pix` and `weights` should point to two
4-element arrays of integers and floating-points, respectively. They
can be reused in multiple calls to avoid heap allocations and speed up
the code.

"""
getinterpolRing
