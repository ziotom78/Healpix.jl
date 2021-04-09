################################################################################

function ringWeightPath(datadir, nside)
    joinpath(datadir, @sprintf("weight_ring_n%05d.fits", nside))
end

function readWeightRing(fileName, nside)
    f = FITSIO.fits_open_table(fileName)
    try
        weights = Array{Float64,1}(undef, 2 * nside)
        FITSIO.fits_read_col(f, 1, 1, 1, weights)

        return weights
    finally
        FITSIO.fits_close_file(f)
    end
end


# compressed size of pixel weights
n_fullweights(nside) = ((3*nside+1)*(nside+1))÷4

"""
    readfullweights(filename::String)

The easiest way to the pixel weight files is to run

```
git clone --depth 1 https://github.com/healpy/healpy-data
```

# Arguments:
- `filename::String`: filename of the full pixel weights

# Returns:
- `Vector{Float64}`: contains the compressed pixel weights
"""
function readfullweights(filename::String)
    f = FITS(filename)
    return read(f[2], "COMPRESSED PIXEL WEIGHTS")
end

"""
    applyweights!(m::Map{T, RingOrder}, wgt::Vector{T}) where T

Apply a pixel weighting to a map for more accurate SHTs.

# Arguments:
- `m::Map{T, RingOrder}`: map to modify
- `wgt::Vector{T}`: compressed pixel weights
"""
function applyweights!(m::Map{T, RingOrder}, wgt::Vector{T}) where T
    nside = m.resolution.nside
    @assert length(wgt) == n_fullweights(nside)
    pix, vpix = 0, 0
    npix = nside2npix(nside)
    for i in 0:(2nside-1)
        shifted = (i<nside-1) || Bool((i+nside)&1)
        qpix=min(nside,i+1)
        odd=Bool(qpix&1)
        wpix=((qpix+1)>>1) + ((odd||shifted) ? 0 : 1)
        psouth=npix-pix-(qpix<<2)
        for j in 0:((qpix<<2)-1)
            j4=j%qpix
            rpix=min(j4,qpix - (shifted ? 1 : 0) - j4)
            m.pixels[pix+j+1] *= one(T) + wgt[vpix+rpix+1]
            if (i!=2*nside-1)
                m.pixels[psouth+j+1] *= one(T) + wgt[vpix+rpix+1]
            end
        end
        pix+=qpix<<2;
        vpix+=wpix;
    end
end


"""
    applyweights!(m::PolarizedMap{T, RingOrder}, wgt::Vector{T}) where T

Apply a pixel weighting to a polarized map for more accurate SHTs.

# Arguments:
- `m::PolarizedMap{T, RingOrder}`: map to modify
- `wgt::Vector{T}`: compressed pixel weights
"""
function applyweights!(m::PolarizedMap{T, RingOrder}, wgt::Vector{T}) where T
    applyweights!(m.i, wgt)
    applyweights!(m.q, wgt)
    applyweights!(m.u, wgt)
end
