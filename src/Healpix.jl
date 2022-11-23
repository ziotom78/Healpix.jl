module Healpix
###
export nsideok, nside2pixarea, nside2resol, nside2order, order2nside
export Resolution, nside2npix, npix2nside, numOfRings
export ang2pixNest, zphi2pixRing, ang2pixRing, pix2angNest, pix2angRing
export vec2pixNest, vec2pixRing, pix2vecNest, pix2vecRing
export pix2ringpos, ring2z
export Order, RingOrder, NestedOrder, AbstractHealpixMap, HealpixMap
export PolarizedHealpixMap, AbstractPolarizedHealpixMap
export ang2vec, vec2ang, ang2pix, pix2ang
export readMapFromFITS, readPolarizedMapFromFITS
export savePixelsToFITS, saveToFITS, conformables
export ringWeightPath, readWeightRing, readFullWeights, applyFullWeights!
export pixelWindowPath, readPixelWindowT, readPixelWindowP
export Alm, numberOfAlms, almIndexL0, almIndex, almExplicitIndex, readAlmFromFITS, writeAlmToFITS
export readClFromFITS, writeClToFITS, cl2dl, dl2cl, synalm!, synalm, synfast!, synfast, anafast
export map2alm, alm2map, map2alm!, alm2map!, adjoint_map2alm!, adjoint_alm2map!, alm2cl, pixwin, gaussbeam, almxfl, almxfl!
export getringinfo!, getringinfo, getinterpolRing
export pix2xyfRing, xyf2pixRing, pix2xyfNest, xyf2pixNest, xyf2loc
export pix2zphiNest, pix2zphiRing, ringAbove, max_pixrad
export interpolate
export ring2nest, nest2ring, ring2nest!, nest2ring!, udgrade
export queryDiscRing, queryStripRing
export boundariesRing, boundariesRing!
export ORDER_MAX, NSIDE_MAX, UNSEEN

using LinearAlgebra
using StaticArrays

using LazyArtifacts
using Random
import CFITSIO

import Libsharp
import Base: getindex, setindex!

"""
A constant commonly used by Healpix libraries to mark «missing» pixels.

This constant is useful if you need compatibility with other Healpix
libraries.

"""
const UNSEEN = -1.6375e+30

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
include("weights.jl")
include("pixelwindow.jl")
include("map_pixelfunc.jl")
include("projections.jl")
include("alm.jl")
include("cl.jl")
include("sphtfunc.jl")
include("mapmaking.jl")
include("query.jl")

end
