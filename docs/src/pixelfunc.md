# Pixel functions

In this section we document the functions that convert from a direction in the sky into a pixel index, and vice versa.

```@docs
ang2vec(theta, phi)
vec2ang(x, y, z)
ang2pixNest(resol::Resolution, theta, phi)
ang2pixRing(resol::Resolution, theta, phi)
pix2angNest(resol::Resolution, pixel)
pix2angRing(resol::Resolution, pixel)
```
