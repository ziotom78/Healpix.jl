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
Healpix.jl implements a functions to read/write the components $C_{\ell}$
from/to a FITS files.

```@docs
readClFromFITS
writeClToFITS
```
