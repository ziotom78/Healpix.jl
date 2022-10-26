```@meta
DocTestSetup = quote
    using Healpix
end
```

# Power Spectrum
Power spectrum components $C_{\ell}$ are encoded as Vector{T}.
You can use the function [`alm2cl`](@ref) to convert a set of $a_{\ell m}$
coefficients into the components $C_\ell$ of the power spectrum.

## Loading and saving power spectrum components
Healpix.jl implements functions to read/write the components $C_{\ell}$
from/to a FITS files.

```@docs
readClFromFITS
writeClToFITS
```

## Converting different power spectrum representation
It's often useful to represent, especially for plotting it, the power spectrum
in the form of $D_{\ell}$. Healpix.jl implements a couple of functions to convert
a power spectrum from/to such a representation.

```@docs
cl2dl
dl2cl
```

## Synthetizing harmonic coefficients from a given power spectrum
Random generate a set of [`Alm`](@ref) from a given power spectrum $C_{\ell}$.
Each harmonic coefficient $a_{\ell m}$ is a realization of a gaussian distribution
with zero mean and $C_{\ell}$ variance.

```@docs
synalm!
synalm
```

## Computing the power spectrum from a map
Compute the (cross-) power spectrum of one (or two) ['HealpixMap'](@ref).

```@docs
anafast
```
