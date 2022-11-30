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
of a [`HealpixMap`](@ref), only living in the harmonic space:

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
The synthesis operation ([`alm2map`](@ref)) is generally referred to
with the matrix operator $\mathrm{Y}$, while his inverse ([`map2alm`](@ref))
with $\mathrm{Y}^{-1}$.

Healpix.jl also implements the two adjoint functions
[`adjoint_alm2map!`](@ref) and [`adjoint_map2alm!`](@ref), represented by
$\mathrm{Y}^{\mathrm{T}}$ and $(\mathrm{Y}^{-1})^\mathrm{T}$ respectively.
While the synthesis operator on a general scalar field $f(\theta, \phi)$
can be defined through an exact summation as $f(\theta, \phi) = \mathrm{Y} \, a_{\ell m} \quad \text{where} \quad f(\theta, \phi) = \sum_{\ell=0}^{\infty} \sum_{m=-\ell}^{\ell} a_{\ell m} Y_{\ell m} (\theta, \phi)$.
The analysis operator is defined through an integral operator as $a_{\ell m} = \mathrm{Y}^{-1} f(\theta, \phi) \quad \text{where} \quad a_{\ell m} = \int_0^{2\pi} \int_0^\pi Y^*_{\ell m}(\theta, \phi)\, f(\theta, \phi) \sin\theta \, d\theta \,d\phi$.
Though, in the real case wherein maps are pixelized, the latter ends
up being approximated through a summation over the pixels.
Here is where the adjoint of the synthesis operator, $\mathrm{Y}^{\mathrm{T}}$,
comes into play. It is defined through: $ \mathrm{Y}^{\mathrm{T}} f(\theta, \phi) \equiv \sum_{i = 1}^{N_{\mathrm{pix}}} Y^*_{\ell m,\, i} \, f_i,$
which is an exact operation. Note that the latter does not give directly the $a_{\ell m}$ coefficients, since $\mathrm{Y}^{-1} \simeq \mathrm{W}\, \mathrm{Y}^{\mathrm{T}}$,
where $\mathrm{W}$ is a diagonal matrix whose non-zero elements are approximately
constant and equal to $4 \pi / N_{\mathrm{pix}}$, depending on the map pixelization.
The latter realtion is also useful to obtain the adjoint of the analysis operator: $(\mathrm{Y}^{-1})^\mathrm{T} = \mathrm{W}^{\mathrm{T}}\,\mathrm{Y} =  \mathrm{W}\,\mathrm{Y}$.


Here is an example:

```@example map2alm
using Healpix # hide
using Random

# Ensure reproducibility by using a fixed seed
Random.seed!(1234)

nside = 8
m = HealpixMap{Float32,RingOrder}(nside)

# Initialize the pixels to random values in the 0…1 range
for i in 1:length(m)
    m[i] = rand(Float32)
end

alm = map2alm(m)

# Go back to pixel space
newm = alm2map(alm, nside)
```

The variable `newm` is a map that is close enough to `m`, yet it is
not exactly the same because of the approximations done by both
`map2alm` and `alm2map`.

```@docs
map2alm!
map2alm
alm2map!
alm2map
adjoint_alm2map!
adjoint_map2alm!
```

## From harmonic coefficients to the power spectrum

You can use the function [`alm2cl`](@ref) to convert a set of $a_{\ell m}$
coefficients into the components $C_\ell$ of the power spectrum.
The pixelization also induces a transfer function, which can be obtained from
[`pixwin`](@ref). A simple Gaussian beam window function in the asymptotic small-beam
limit can be computed with [`gaussbeam`](@ref).

```@docs
alm2cl
pixwin
gaussbeam
```

## Multiplying a set of Alm by a generic function of $\ell$

You can use the function [`almxfl`](@ref) (or [`almxfl!`](@ref)) to multiply (in-place)
 a set of $a_{\ell m}$ coefficients by an $\ell$-dependent generic function $f_\ell$.

```@docs
almxfl
almxfl!
```

## Loading and saving harmonic coefficients

```@docs
readAlmFromFITS
writeAlmToFITS
```

## Alm Indexing

You can use [`almExplicitIndex`](@ref) to compute the so-called explicit indexing.
It is exploited for instance in [`readAlmFromFITS`](@ref) and [`writeAlmToFITS`](@ref).

```@docs
almExplicitIndex
```

The following functions can be used, in an analogous way as [`eachindex`](@ref),
in the case of arrays, to obtain sets of indexes or values of $\ell$ and $m$.

On the same line as `eachindex`, these can be very useful when implementing for-cycles
over `Alm` objects.
Here is an example of how to exploit [`each_ell_m`](@ref) to print explicitly
the major-m ordering of a set of complex-stored `Alm`:

```@example each_ell_m

using Healpix # hide
using Random

# Initialize a random set of alm
alm = Alm(4,4, randn(ComplexF64, numberOfAlms(4,4)))

# Print it's values along with each corresponding l and m values
i=1
for (l,m) in each_ell_m(alm)
    a_lm = alm.alm[i]
    print("ℓ = $l, |m| = $m: a_ℓm = $a_lm \n")
    i+=1
end
```

```@docs
each_ell
each_ell_idx
each_m
each_m_idx
each_ell_m
```

## Full Pixel Weights

The default [`map2alm`](@ref) uses iteration to obtain an accurate transform.
One can instead apply a pixel weight to compute an accurate transform in a single
pass, like quadrature. The easiest way to the pixel weight files is to run

```
git clone --depth 1 https://github.com/healpy/healpy-data
```

These weights are in a compressed format that is read with [`readFullWeights`](@ref)
and multiplied into a map with [`applyFullWeights!`](@ref).

```julia
nside = 32
compressed_weights = Healpix.readFullWeights(
    "healpix_full_weights_nside_$(lpad(nside,4,'0')).fits")
m = Healpix.HealpixMap{Float64,Healpix.RingOrder}(ones(Healpix.nside2npix(nside)))
Healpix.applyFullWeights!(m, compressed_weights)
alm = Healpix.map2alm(m; niter=0)
```

The subsequent [`map2alm`](@ref) only needs `niter=0`.

```@docs
readFullWeights
applyFullWeights!
```
