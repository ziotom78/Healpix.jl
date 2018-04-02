# Pixel functions

In this section we document the functions that convert from a direction in the sky into a pixel index, and vice versa.

First of all, Healpix.jl implements the most basic functions to convert between spherical and Cartesian coordinates. Note that Healpix uses **co-latitude** instead of **latitude**:

```@repl
using Healpix # hide
ang2vec(0.0, 0.0)
vec2ang(0.0, 0.0, 1.0)
```

More interesting functions return the index of the pixel on a Healpix-tessellated sphere. For these functions to work, you have to provide a [`Resolution`](@ref) object:

```@example
using Healpix # hide
res = Resolution(16)
ang2pixRing(res, π/2, 0)
ang2pixNest(res, π/2, 0)
```

```@docs
ang2vec(theta, phi)
vec2ang(x, y, z)
ang2pixNest(resol::Resolution, theta, phi)
ang2pixRing(resol::Resolution, theta, phi)
pix2angNest(resol::Resolution, pixel)
pix2angRing(resol::Resolution, pixel)
```
