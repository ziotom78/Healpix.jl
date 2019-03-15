export tod2map, combinemaps!

@doc raw"""
    tod2map{T,O}(pixidx, tod::Array{T}; nside=128) :: (map, hits)

Create a binned map for a TOD and return a tuple containing the map
itself and the hit map.
"""
function tod2map(pixidx, tod::Array{T}; nside=128, ordering=Healpix.RingOrder) where {T}
    @assert length(pixidx) == length(tod)
    
    binnedmap = Map{T,ordering}(nside)
    hitmap = Map{Int,ordering}(nside)
    
    @inbounds for i in eachindex(pixidx)
        if !ismissing(tod[i])
            binnedmap.pixels[pixidx[i]] += tod[i]
            hitmap.pixels[pixidx[i]] += 1
        end
    end
    
    @inbounds for i in eachindex(binnedmap.pixels)
        if hitmap.pixels[i] > 0
            binnedmap.pixels[i] /= hitmap.pixels[i]
        end
    end
    
    (binnedmap, hitmap)
end

@doc raw"""
    combinemaps{T, O, H}(destmap::Map{T, O}, desthitmap::Map{H, O}, othermap::Map{T, O}, otherhitmap::Map{H, O})

Sum "othermap" to "destmap", assuming that both maps have been
produced by binning TODs. The parameters `desthitmap` and
`otherhitmap` are the two hit maps. At the end of the call, `destmap`
and `desthitmap` are updated.
"""
function combinemaps!(destmap::Map{T,O},
                      desthitmap::Map{H,O},
                      othermap::Map{T,O},
                      otherhitmap::Map{H,O}) where {T,O,H}
    (conformables(destmap, othermap) 
     && conformables(desthitmap, otherhitmap) 
     && conformables(destmap, desthitmap)) || throw(DomainError()) 

    for i in eachindex(destmap.pixels)
        (pix1, pix2) = (destmap.pixels[i], othermap.pixels[i])
        (hit1, hit2) = (desthitmap.pixels[i], otherhitmap.pixels[i])
        hitsum = hit1 + hit2
        if hitsum > 0
            destmap.pixels[i] = (pix1 * hit1 + pix2 * hit2) / hitsum
        end
        desthitmap.pixels[i] = hitsum
    end
end
