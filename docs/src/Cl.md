```@meta
DocTestSetup = quote
    using Healpix
end
```

# Power Spectrum
Power spectrum components $C_{\ell}$ are encoded as Vector{T}.
Healpix.jl implements functions to perform sht operations on power spectra,
e.g. to obtain a map or a set of ['Alm'](@ref), as well as writing/reading a
power spectrum from a FITS file.
The functions ending with `!` are *mutating* functions, which means that they
assume that the result must be saved in a preallocated variable; they are
space- and time-efficient and should be used when you want your code to be performant,
or when you plan to apply the same operation several times (e.g., in a Monte Carlo
simulation).

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

## Synthesizing harmonic coefficients from power spectrum
Random generate a set of [`Alm`](@ref) from a given power spectrum $C_{\ell}$.
Each harmonic coefficient $a_{\ell m}$ is a realization of a gaussian distribution
with zero mean and $C_{\ell}$ variance.

```@docs
synalm!
synalm
```

## Generating a map from power spectrum
Synthesize a set of [`Alm`](@ref) through ['synalm'](@ref) and generates a map
from it through ['alm2map'](@ref).

```@docs
synfast!
synfast
```

## Computing the power spectrum from a map
Compute the (cross-) power spectrum of one (or two) ['HealpixMap'](@ref).

```@docs
anafast
```
