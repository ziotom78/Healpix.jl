```@meta
DocTestSetup = quote
    using Healpix
end
```

# Working with resolutions

A Healpix tessellation is parametrized by a number, called `NSIDE`,
which must be a positive power of 2. It is related to the number of
pixels $N$ in the maps by the simple equation $N = 12
\mathrm{NSIDE}^2$, and it is therefore related to the resolution of
the pixelization. Any function working on a Healpix tessellation needs
to receive the value of `NSIDE`. Healpix.jl provides a wrapper around
this parameter, the [`Resolution`](@ref) type, which internally keeps
a number of precomputed coefficients to accelerate calculations.

The following example prints a table containing details about a few
Healpix resolutions:

```@example
using Healpix # hide
using Printf

@printf("%-6s\t%-12s\t%-12s\t%-12s\n",
        "NSIDE",
        "#pix",
        "#pix per face",
        "solid angle")
for poweroftwo in [0, 1, 2, 3, 4, 5]
    res = Resolution(2 ^ poweroftwo)
    @printf("%6d\t%12d\t%12d\t%12.4f\n",
            res.nside,
            res.numOfPixels,
            res.pixelsPerFace,
            4π / res.numOfPixels)
end
```

There is an upper limit to the value of the `NSIDE` parameter, which
is encoded in the constant [`NSIDE_MAX`](@ref). The value is determined
at runtime according to the size of the `Int` type; on 32-bit machines
it is 8192 (``2^{13}``), while on 64-bit machines it is
536870912 (``2^{29}``).

```@docs
Resolution
Resolution(nside::Integer)
nsideok(nside::Integer)
nside2npix(nside::Integer)
npix2nside(npix::Integer)
nside2pixarea(nside::Integer)
nside2resol(nside::Integer)
nside2order(nside::Integer)
order2nside(order::Integer)
ORDER_MAX
NSIDE_MAX
```
