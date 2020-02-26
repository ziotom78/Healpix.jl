module Healpix
###
export nsideok, nside2pixarea, nside2resol
export Resolution, nside2npix, npix2nside
export ang2pixNest, ang2pixRing, pix2angNest, pix2angRing
export vec2pixNest, vec2pixRing, pix2vecNest, pix2vecRing
export pix2ringpos
export Order, RingOrder, NestedOrder, Map, PolarizedMap
export ang2vec, vec2ang, ang2pix, pix2ang
export readMapFromFITS, savePixelsToFITS, saveToFITS, conformables
export ringWeightPath, readWeightRing
export pixelWindowPath, readPixelWindowT, readPixelWindowP
export Alm, numberOfAlms, almIndexL0, almIndex, readAlmFromFITS
export getringinfo!, getringinfo, getinterpolRing
export pix2xyfRing, xyf2pixRing, pix2xyfNest, xyf2pixNest
export pix2zphiNest, pix2zphiRing, ringAbove
export ring2nest, nest2ring

import FITSIO
import Base: getindex, setindex!

include("nside.jl")
include("math.jl")
include("datatables.jl")
include("resolution.jl")
include("pixelfunc.jl")
include("interp.jl")
include("xyf.jl")
include("map.jl")
include("polarizedmap.jl")
include("map_io.jl")
include("conformables.jl")
include("ringweights.jl")
include("pixelwindow.jl")
include("map_pixelfunc.jl")
include("projections.jl")
include("alm.jl")
include("mapmaking.jl")

end
