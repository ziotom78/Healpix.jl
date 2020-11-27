```@meta
DocTestSetup = quote
    using Healpix
end
```

# Spherical harmonics

Starting from version 2.4, Healpix.jl implements generalized Fourier
transformations through the `libsharp` library, to convert a map from
its pixel-space representation to its decomposition in spherical
harmonics. This has multiple applications, the most relevant being the
analysis of Cosmic Microwave Background maps and the efficient
computation of convolution operators.

Everything revolves around the [`Alm`](@ref) type, which encodes a set of
spherical harmonics and is thus conceptually equivalent to the concept
of a [`Map`](@ref), only living in the harmonic space:

```@docs
Alm
```

In the general case, the number of coefficients in a spherical
harmonic expansion is infinite. For obvious reasons, Healpix.jl only
allows to store band-limited expansions. The function
[`numberOfAlms`](@ref) returns the number of floating-point numbers
used to store the expansion, as a function of the maximum value for
$\ell$ and $m$.

```@docs
numberOfAlms
```

## Converting between pixel space and harmonic space

Healpix.jl implements the four functions [`alm2map`](@ref),
[`map2alm`](@ref), [`alm2map!`](@ref), and [`map2alm!`](@ref) to
convert a map from a pixel-space representation to the harmonic space
and vice-versa. The functions ending with `!` are *mutating*
functions, which means that they assume that the result must be saved
in a preallocated variable; they are space- and time-efficient and
should be used when you want your code to be performant, or when you
plan to apply the same operation several times (e.g., in a Monte Carlo
simulation).

Here is an example:

```@example map2alm
using Healpix
using Random

# Ensure reproducibility by using a fixed seed
Random.seed!(1234)

nside = 8
m = Map{Float32,RingOrder}(nside)

# Initialize the pixels to random values in the 0…1 range
for i in 1:length(m)
    m[i] = rand(Float32)
end

alm = map2alm(m)

# Go back to pixel space
newm = alm2map(alm)
```

The variable `newm` is a map that is close enough to `m`, yet it is
not exactly the same because of the approximations done by both
`map2alm` and `alm2map`.

## From harmonic coefficients to the power spectrum

You can use the function `alm2cl` to convert a set of $a_{\ell m}$
coefficients into the components $C_\ell$ of the power spectrum.

```@docs
alm2cl
```

## Loading and saving harmonic coefficients

```@docs
readAlmFromFITS
```
