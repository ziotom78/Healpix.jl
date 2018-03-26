# Working with resolutions

A Healpix tessellation is parametrized by a number, called `NSIDE`, which must be a positive power of 2. It is related to the number of pixels $N$ in the maps by the simple equation $N = 12 \mathrm{NSIDE}^2$, and it is therefore related to the resolution of the pixelization. Any function working on a Healpix tessellation needs to receive the value of `NSIDE`. Healpix.jl provides a wrapper around this parameter, the [`Resolution`](@ref) type, which internally keeps a number of precomputed coefficients to accelerate calculations.

```@docs
Resolution
Resolution(nside)
nside2npix(nside)
npix2nside(npix)
```
