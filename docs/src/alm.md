# Spherical harmonics

The support for spherical harmonics in Healpix.jl is still woefully inadequate. Only a few functions to load and store harmonic coefficients are available. Everything revolves around the `Alm` type:

```@docs
Alm
```

The number of coefficients in a spherical harmonic expansion is infinite. For obvious reasons, Healpix.jl only allows to store band-limited expansions. The function `numberOfAlms` returns the number of floating-point numbers used to store the expansion, as a function of the maximum value for $\ell$ and $m$.

```@docs
numberOfAlms
```

## Loading and saving harmonic coefficients

```@docs
readAlmFromFITS
```