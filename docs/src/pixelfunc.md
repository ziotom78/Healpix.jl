```@meta
DocTestSetup = quote
    using Healpix
end
```

# Pixel functions

In this section we document the functions that convert from a
direction in the sky into a pixel index, and vice versa.

First of all, Healpix.jl implements the most basic functions to
convert between spherical and Cartesian coordinates. Note that Healpix
uses **co-latitude** instead of **latitude**:

```@repl pixelexample1
using Healpix # hide
ang2vec(0.0, 0.0)
vec2ang(0.0, 0.0, 1.0)
```

More interesting functions return the index of the pixel on a
Healpix-tessellated sphere. For these functions to work, you have to
provide a [`Resolution`](@ref) object:

```@repl pixelexample1
res = Resolution(16)
ang2pixRing(res, π/2, 0)
ang2pixNest(res, π/2, 0)
```

# Ring functions

The Healpix projection has the advantage of storing pixels along
iso-latitude rings; this allows to implement efficient
spherical-transform functions. Healpix.jl provides a number of
functions that manage rings. Many of them use the `RingInfo`
structure, which encodes details about a ring.

```@docs
RingInfo
getringinfo
getringinfo!
getinterpolRing
```

# Reference

```@docs
ang2vec(theta, phi)
vec2ang(x, y, z)
ang2pixNest(resol::Resolution, theta, phi)
ang2pixRing(resol::Resolution, theta, phi)
pix2angNest(resol::Resolution, pixel)
pix2angRing(resol::Resolution, pixel)
ring2nest(resol::Resolution, ipix)
nest2ring(resol::Resolution, ipix)
pix2ringpos(resol::Resolution, pixel)
pix2xyfNest(resol::Resolution, ipix)
pix2xyfRing(resol::Resolution, ipix)
xyf2pixNest(resol::Resolution, ix, iy, facenum)
xyf2pixRing(resol::Resolution, ix, iy, facenum)
pix2zphiRing(res::Resolution, pix)
pix2zphiNest(res::Resolution, pix)
ringAbove
```
