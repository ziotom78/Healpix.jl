# Healpix

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ziotom78.github.io/Healpix.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ziotom78.github.io/Healpix.jl/dev)
[![Build Status](https://github.com/ziotom78/Healpix.jl/workflows/Unit%20tests/badge.svg)](https://github.com/ziotom78/Healpix.jl/actions?query=workflow%3A%22Unit+tests%22)
[![Codecov](https://codecov.io/gh/ziotom78/Healpix.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ziotom78/Healpix.jl)
<a href="https://ascl.net/2109.028"><img src="https://img.shields.io/badge/ascl-2109.028-blue.svg?colorB=262255" alt="ascl:2109.028" /></a>

Healpix.jl is a set of Julia functions that implement the
[HEALPix](https://en.wikipedia.org/wiki/HEALPix) algorithms to
pixelate a 2-sphere. The HEALPix pixelisation has a number of features
that make it efficient for the following tasks:

-   All pixels have the same area (*equal-area pixelisation*);
-   It is efficient to compute the Fourier decomposition in spherical
    harmonics of a discretized signal on the sphere (using the `RING`
    ordering scheme);
-   It is efficient to find neighbour pixels (using the `NESTED`
    ordering scheme);
-   The base pixelisation divides the 2-sphere in 12 pixels, but each
    of these pixels can be divided into four sub-pixels, and the
    process can be continued; in principle, the size of a pixel can be
    reduced to arbitrarily small values.

HEALPix is widely used in cosmology to store maps of the CMB
temperature and polarization anisotropies, which is the field studied
by the authors of this package.


# Supported platforms

The purpose of Healpix.jl, is to implement a Julia-only library,
instead of providing bindings to the original implementation of the
[C/C++/Fortran/Python Healpix library](http://healpix.jpl.nasa.gov/).
This means that, unlike the original implementation, Healpix.jl is
**fully supported on Windows systems**, among with Linux and Mac OS X.


## Installation

From the Julia REPL, run

````julia
import Pkg
Pkg.add("Healpix.jl")
````


## Usage examples

Here are some code snippets that show how to use `Healpix.jl`. It is
interesting to have a look at
[test/runtests.jl](https://github.com/ziotom78/Healpix.jl/blob/master/test/runtests.jl)
as well.

Refer to the
[documentation](https://ziotom78.github.io/Healpix.jl/stable) for more
examples.


### Dealing with resolutions

The resolution of a HEALPix map is uniquely determined by the `NSIDE`
parameter. Healpix.jl precalculates a number of values from `NSIDE` to
save time during computations; such values are saved in a
`Healpix.Resolution` object:

`````julia
import Healpix
res = Healpix.Resolution(256)
print("The pixel index is $(Healpix.ang2pixNest(res, 0.1, 0.2))\n")
`````


### Reading a map from a FITS file

This snippet loads a map named `planck_70GHz.fits` into an array of
64-bit floating-point numbers:

`````julia
import Healpix

m = Healpix.readMapFromFITS("planck_70GHz.fits", 1, Float64)
print("average: $(mean(m.pixels))\n")
`````


## How to contribute

See the document [CONTRIBUTING.md](https://github.com/ziotom78/Healpix.jl/blob/master/CONTRIBUTING.md).

## Citing Healpix.jl

If you have used Healpix.jl in your research and want to acknowledge the library in your academic publications, you can reference its  [ASCL entry](http://ascl.net/code/v/3025) through the code `ascl:2109.028`:

> Tomasi M., Li Z. 2021 Healpix.jl: Julia-only port of the HEALPix library, 0.30, Astrophysics Source Code Library ascl:2109.028

Here is a BibTeX entry ready to be used, generated by [ADS](https://ui.adsabs.harvard.edu/abs/2021ascl.soft09028T/abstract):

```
@MISC{2021ascl.soft09028T,
       author = {{Tomasi}, Maurizio and {Li}, Zack},
        title = "{Healpix.jl: Julia-only port of the HEALPix library}",
     keywords = {Software},
      version = {3.0},     
         year = 2021,
        month = sep,
          eid = {ascl:2109.028},
        pages = {ascl:2109.028},
archivePrefix = {ascl},
       eprint = {2109.028},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2021ascl.soft09028T},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

Please update the version number (`version = {…}`) accordingly.

See the [ASCL guidelines](http://ascl.net/home/getwp/351) for more information.

## License

Healpix.jl is released under the GPL license. Versions before 2.3.0
were released under a MIT license, but this was considered
incompatible with the way the code has been written
([#15](https://github.com/ziotom78/Healpix.jl/issues/15)).
