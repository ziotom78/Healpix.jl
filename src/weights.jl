################################################################################

function ringWeightPath(datadir, nside)
    joinpath(datadir, @sprintf("weight_ring_n%05d.fits", nside))
end

function readWeightRing(fileName, nside)
    f = CFITSIO.fits_open_table(fileName)
    try
        weights = Array{Float64,1}(undef, 2 * nside)
        CFITSIO.fits_read_col(f, 1, 1, 1, weights)

        return weights
    finally
        CFITSIO.fits_close_file(f)
    end
end


# compressed size of pixel weights
n_fullweights(nside) = ((3 * nside + 1) * (nside + 1)) ÷ 4

"""
    readFullWeights(filename::String)

Read the set of pixel weights used to compute the generalized Fourier
transform of a map.

These weights are usually precomputed; you can download the ones
available in the [Healpy
repository](https://github.com/healpy/healpy-data) using the following
command:

```
git clone --depth 1 https://github.com/healpy/healpy-data
```

# Arguments:
- `filename::String`: filename of the full pixel weights

# Returns:
- `Vector{Float64}`: contains the compressed pixel weights
"""
function readFullWeights(filename::String)
    f = CFITSIO.fits_open_table(filename)
    result = Vector{Float64}(undef, CFITSIO.fits_get_num_rows(f))
    CFITSIO.fits_read_col(f, 1, 1, 1, result)
    result
end

"""
    applyFullWeights!(m::HealpixMap{T, RingOrder}, [wgt::Vector{T}]) where T

Apply a pixel weighting to a map for more accurate SHTs. Note that 
this only helps for `lmax<=1.5*Nside`. If this is not the case, the 
pixel weights may do more harm than good.

Pixel weights are automatically downloaded if not specified. 

# Arguments:
- `m::HealpixMap{T, RingOrder}`: map to modify
- `wgt::Vector{T}` (optional): compressed pixel weights. If not specified, this routine will
    look for weights in artifacts.
"""
function applyFullWeights!(m::HealpixMap{T,RingOrder}, wgt::Vector{T}) where T
    nside = m.resolution.nside
    @assert length(wgt) == n_fullweights(nside)
    pix, vpix = 0, 0
    npix = nside2npix(nside)
    for i in 0:(2nside - 1)
        shifted = (i < nside - 1) || Bool((i + nside) & 1)
        qpix = min(nside, i + 1)
        odd = Bool(qpix & 1)
        wpix = ((qpix + 1) >> 1) + ((odd || shifted) ? 0 : 1)
        psouth = npix - pix - (qpix << 2)
        for j in 0:((qpix << 2) - 1)
            j4 = j % qpix
            rpix = min(j4, qpix - (shifted ? 1 : 0) - j4)
            m.pixels[pix + j + 1] *= one(T) + wgt[vpix + rpix + 1]
            if (i != 2 * nside - 1)
                m.pixels[psouth + j + 1] *= one(T) + wgt[vpix + rpix + 1]
            end
        end
        pix += qpix << 2;
        vpix += wpix;
    end
end


function applyFullWeights!(m::HealpixMap{T,RingOrder}) where T
    nside = m.resolution.nside
    
    if nside ∈ (32, 64, 128, 256, 512, 1024, 2048)
        path = artifact"fullpixelweights_32_2048"
    elseif nside == 4096
        path = artifact"fullpixelweights_4096"
    elseif nside == 8192
        path = artifact"fullpixelweights_8192"
    else
        throw(ArgumentError("Unsupported nside $(nside)"))
    end

    nside_str = lpad(nside, 4, '0')
    wgt = readFullWeights(joinpath(path, "healpix_full_weights_nside_$(nside_str).fits"))
    applyFullWeights!(m, wgt)
end


"""
    applyFullWeights!(m::PolarizedHealpixMap{T, RingOrder}, [wgt::Vector{T}]) where T

Apply a pixel weighting to a polarized map for more accurate SHTs.

# Arguments:
- `m::PolarizedHealpixMap{T, RingOrder}`: map to modify
- `wgt::Vector{T}` (optional): compressed pixel weights. If not specified, an artifact 
        will be sought.
"""
function applyFullWeights!(m::PolarizedHealpixMap{T,RingOrder}, wgt::Vector{T}) where T
    applyFullWeights!(m.i, wgt)
    applyFullWeights!(m.q, wgt)
    applyFullWeights!(m.u, wgt)
end

function applyFullWeights!(m::PolarizedHealpixMap{T,RingOrder}) where T
    applyFullWeights!(m.i)
    applyFullWeights!(m.q)
    applyFullWeights!(m.u)
end
